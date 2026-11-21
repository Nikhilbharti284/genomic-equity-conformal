import csv, gzip, json, math, os, sqlite3, time, urllib.request
from resource_logger import ResourceLogger

AM_URL = "https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz?download=1"

def _num(x):
    if x is None: return 0.0
    if isinstance(x, list):
        vals = [_num(v) for v in x]
        return max(vals) if vals else 0.0
    try:
        v = float(x)
        return 0.0 if math.isnan(v) else v
    except Exception:
        return 0.0

def stream_alphamissense(target_keys):
    """Streams 643 MB CC-BY-4.0 AlphaMissense_hg38.tsv.gz in memory without unpacking to disk."""
    print(f"1/2 Streaming CC-BY-4.0 AlphaMissense_hg38.tsv.gz for {len(target_keys)} target variants...")
    found = {}
    req = urllib.request.Request(AM_URL, headers={"User-Agent": "genomic-equity-conformal/1.0"})
    t0 = time.time()
    scanned = 0
    with urllib.request.urlopen(req, timeout=300) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            for raw in gz:
                if raw.startswith(b"#"): continue
                scanned += 1
                parts = raw.decode("utf-8", errors="ignore").rstrip("\n").split("\t")
                if len(parts) < 9: continue
                # #CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class
                key = f"{parts[0]}:g.{parts[1]}{parts[2]}>{parts[3]}"
                if key in target_keys:
                    try:
                        found[key] = float(parts[8])
                    except Exception:
                        pass
                if scanned % 10_000_000 == 0:
                    print(f"    Scanned {scanned//1_000_000}M rows ({round(time.time()-t0,1)}s) | Matched: {len(found)}...")
    print(f"    AlphaMissense matched: {len(found)}/{len(target_keys)} variants in {round(time.time()-t0,1)}s")
    return found

def fetch_gnomad_hg38(hgvs_ids):
    print(f"2/2 Fetching hg38 gnomAD ancestry AFs for {len(hgvs_ids)} variants via MyVariant.info...")
    out = {}
    url = "https://myvariant.info/v1/variant?assembly=hg38&fields=gnomad_exome.af,gnomad_genome.af,dbnsfp.gnomad_exomes"
    for i in range(0, len(hgvs_ids), 250):
        payload = json.dumps({"ids": hgvs_ids[i:i+250]}).encode("utf-8")
        req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json", "User-Agent": "genomic-equity-conformal/1.0"})
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                for r in json.loads(resp.read().decode("utf-8")):
                    q = r.get("query")
                    if not q or r.get("notfound"): continue
                    for k in ("gnomad_exome", "gnomad_genome"):
                        sub = r.get(k)
                        if isinstance(sub, list) and sub: sub = sub[0]
                        if isinstance(sub, dict) and isinstance(sub.get("af"), dict):
                            a = sub["af"]
                            vals = {"af_nfe": _num(a.get("af_nfe")), "af_sas": _num(a.get("af_sas")), "af_afr": _num(a.get("af_afr")), "af_all": _num(a.get("af"))}
                            if max(vals.values()) > 0:
                                out[q] = vals
                                break
        except Exception as e:
            print(f"    [warn] batch {i}: {e}")
        if (i + 250) % 2500 == 0:
            print(f"    gnomAD progress: {min(i+250, len(hgvs_ids))}/{len(hgvs_ids)}...")
    return out

def main():
    logger = ResourceLogger()
    with open("data/v1_clinvar_120genes.csv") as f:
        rows = list(csv.DictReader(f))
    target_keys = {r["hgvs_g"] for r in rows}

    am_map = stream_alphamissense(target_keys)
    af_map = fetch_gnomad_hg38([r["hgvs_g"] for r in rows])

    os.makedirs("artifacts", exist_ok=True)
    db_path = "artifacts/v1_commercial_features.sqlite"
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("""
        CREATE TABLE IF NOT EXISTS var_features (
            variant_id TEXT PRIMARY KEY,
            hgvs_g TEXT,
            am_score REAL,
            af_nfe REAL, af_sas REAL, af_afr REAL, af_all REAL,
            ancestry_group TEXT
        )
    """)
    counts = {"EUR": 0, "SAS": 0, "AFR": 0, "ULTRA_RARE": 0}
    for r in rows:
        af = af_map.get(r["hgvs_g"], {"af_nfe": 0.0, "af_sas": 0.0, "af_afr": 0.0, "af_all": 0.0})
        if af["af_sas"] > 0 and af["af_sas"] >= max(af["af_nfe"], af["af_afr"]): grp = "SAS"
        elif af["af_afr"] > 0 and af["af_afr"] >= max(af["af_nfe"], af["af_sas"]): grp = "AFR"
        elif af["af_nfe"] > 0: grp = "EUR"
        else: grp = "ULTRA_RARE"
        counts[grp] += 1
        conn.execute(
            "INSERT OR REPLACE INTO var_features VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            (r["variant_id"], r["hgvs_g"], am_map.get(r["hgvs_g"], -1.0), af["af_nfe"], af["af_sas"], af["af_afr"], af["af_all"], grp)
        )
    conn.commit()
    conn.close()
    stats = logger.log("v1_ingest_am_gnomad", {"n_rows": len(rows), "am_matched": len(am_map), "ancestry_counts": counts})
    print("Ancestry distribution (120 genes):", counts)
    print("Run summary:", stats)

if __name__ == "__main__":
    main()
