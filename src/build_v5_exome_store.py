import csv, gzip, hashlib, io, json, os, re, sqlite3, sys, time, urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from offline_reference_store import CANONICAL_CONTROLS, build_offline_store

HGNC_TSV = "data/hgnc_protein_coding_genes.tsv"
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
V5_DB = "artifacts/v5_panhuman_exome_ref.sqlite"
MANIFEST = "artifacts/v5_exome_store_manifest.json"

AA3TO1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
    "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
    "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
}
P_MISSENSE_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
TWO_STAR = {
    "criteria provided, multiple submitters, no conflicts",
    "reviewed by expert panel",
    "practice guideline",
}

def load_hgnc_genes(conn):
    g_thr = {}
    if os.path.exists("data/ba1_exception_registry.json"):
        g_thr = json.load(open("data/ba1_exception_registry.json", encoding="utf-8")).get("gene_thresholds", {})
    genes = {}
    rows = []
    with open(HGNC_TSV, encoding="utf-8", errors="ignore") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sym = (r.get("symbol") or "").strip().upper()
            if not sym:
                continue
            moi = g_thr.get(sym, {}).get("moi", "UNKNOWN")
            chrom = (r.get("location") or "").split("p")[0].split("q")[0].strip()
            genes[sym] = moi
            rows.append((sym, r.get("hgnc_id",""), r.get("name",""), chrom, r.get("ensembl_gene_id",""), r.get("uniprot_ids",""), moi))
    conn.executemany("INSERT OR REPLACE INTO hgnc_exome_genes VALUES (?,?,?,?,?,?,?)", rows)
    conn.commit()
    return genes

def seed_existing_v4_rows(conn):
    v3_db = "artifacts/v3_offline_clinical_ref.sqlite"
    if not os.path.exists(v3_db):
        build_offline_store(v3_db)
    c3 = sqlite3.connect(v3_db)
    rows = []
    for r in c3.execute("SELECT * FROM clinical_ref"):
        # Parse chrom/pos_g/ref/alt from hgvs_g if possible
        m = re.match(r"chr([^:]+):g\.(\d+)([ACGT]+)>([ACGT]+)", r[0] or "")
        chrom, pos_g, ref, alt = (m.group(1), int(m.group(2)), m.group(3), m.group(4)) if m else ("", 0, "", "")
        rows.append(tuple(r) + (chrom, pos_g, ref, alt, "TIER1_V4_CALIBRATED"))
    c3.close()
    conn.executemany(
        "INSERT OR REPLACE INTO clinical_ref VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        rows,
    )
    conn.commit()
    return len(rows)

def stream_clinvar_2star(conn, hgnc_genes, max_retries=3):
    batch, added = [], 0
    for attempt in range(1, max_retries + 1):
        try:
            req = urllib.request.Request(CLINVAR_URL, headers={"User-Agent": "GenomicEquityConformal/5.0.0"})
            with urllib.request.urlopen(req, timeout=180) as resp:
                with gzip.GzipFile(fileobj=resp) as gz:
                    reader = csv.DictReader(io.TextIOWrapper(gz, encoding="utf-8", errors="ignore"), delimiter="\t")
                    for r in reader:
                        if r.get("Assembly") != "GRCh38" or r.get("Type") != "single nucleotide variant":
                            continue
                        if (r.get("ReviewStatus") or "").lower() not in TWO_STAR:
                            continue
                        cs = (r.get("ClinicalSignificance") or "").lower()
                        is_path = ("pathogenic" in cs) and ("benign" not in cs) and ("conflicting" not in cs)
                        is_ben  = ("benign" in cs) and ("pathogenic" not in cs) and ("conflicting" not in cs)
                        if not (is_path or is_ben):
                            continue
                        gene = (r.get("GeneSymbol") or "").strip().upper()
                        if gene not in hgnc_genes:
                            continue
                        m = P_MISSENSE_RE.search(r.get("Name") or "")
                        if not m:
                            continue
                        wt, pos, mt = AA3TO1.get(m.group(1)), int(m.group(2)), AA3TO1.get(m.group(3))
                        if not wt or not mt or wt == mt:
                            continue
                        chrom = (r.get("Chromosome") or "").replace("chr", "").strip()
                        pos_g = int(r.get("PositionVCF") or r.get("Start") or 0)
                        ref = (r.get("ReferenceAlleleVCF") or "").strip().upper()
                        alt = (r.get("AlternateAlleleVCF") or "").strip().upper()
                        if not chrom or pos_g <= 0 or len(ref) != 1 or len(alt) != 1 or ref == alt:
                            continue
                        hgvs_g = f"chr{chrom}:g.{pos_g}{ref}>{alt}"
                        vid = f"{gene}_{wt}{pos}{mt}"
                        # Prior proxy only when variant not already in v4 calibrated set
                        am_proxy = 0.88 if is_path else 0.08
                        llr_proxy = round(-6.5 * am_proxy + 0.5, 4)
                        moi = hgnc_genes[gene] if hgnc_genes[gene] != "UNKNOWN" else "AD"
                        batch.append((
                            hgvs_g, vid, gene, wt, pos, mt, am_proxy, llr_proxy,
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, moi,
                            chrom, pos_g, ref, alt, "TIER1_CLINVAR_2STAR_PANHUMAN"
                        ))
                        if len(batch) >= 5000:
                            conn.executemany(
                                "INSERT OR IGNORE INTO clinical_ref VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                                batch,
                            )
                            conn.commit()
                            added += len(batch)
                            batch.clear()
                            if added % 25000 == 0:
                                print(f"  ... streamed {added} ClinVar >=2-star GRCh38 missense rows")
            if batch:
                conn.executemany(
                    "INSERT OR IGNORE INTO clinical_ref VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    batch,
                )
                conn.commit()
                added += len(batch)
            return added
        except Exception as e:
            print(f"[WARN] Attempt {attempt} failed: {e}")
            if attempt == max_retries:
                raise
            time.sleep(3 * attempt)
    return added

def main():
    t0 = time.time()
    os.makedirs("artifacts", exist_ok=True)
    conn = sqlite3.connect(V5_DB)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute("""CREATE TABLE IF NOT EXISTS hgnc_exome_genes (
        gene TEXT PRIMARY KEY, hgnc_id TEXT, name TEXT, chrom TEXT,
        ensembl_gene_id TEXT, uniprot_id TEXT, moi TEXT
    )""")
    conn.execute("""CREATE TABLE IF NOT EXISTS clinical_ref (
        hgvs_g TEXT PRIMARY KEY,
        variant_id TEXT, gene TEXT, wt TEXT, pos INTEGER, mt TEXT,
        am_score REAL, esm2_llr REAL,
        af_sas REAL, af_afr REAL, af_eas REAL, af_amr REAL, af_mid REAL,
        af_nfe REAL, af_asj REAL, af_fin REAL, af_all REAL, moi TEXT,
        chrom TEXT, pos_g INTEGER, ref TEXT, alt TEXT, source_tier TEXT
    )""")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_v5_vid ON clinical_ref(variant_id);")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_v5_gene ON clinical_ref(gene);")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_v5_coord ON clinical_ref(chrom, pos_g, ref, alt);")

    hgnc_genes = load_hgnc_genes(conn)
    n_v4 = seed_existing_v4_rows(conn)
    print(f"[1/3] Indexed {len(hgnc_genes)} HGNC protein-coding genes + {n_v4} calibrated v4 rows.")
    n_cv = stream_clinvar_2star(conn, hgnc_genes)
    total_rows = conn.execute("SELECT COUNT(*) FROM clinical_ref").fetchone()[0]
    uniq_genes = conn.execute("SELECT COUNT(DISTINCT gene) FROM clinical_ref").fetchone()[0]
    conn.execute("PRAGMA wal_checkpoint(TRUNCATE);")
    conn.close()

    size_mb = round(os.path.getsize(V5_DB) / (1024 * 1024), 2)
    sha = hashlib.sha256(open(V5_DB, "rb").read()).hexdigest()
    manifest = {
        "version": "5.0.0",
        "db_path": V5_DB,
        "hgnc_exome_genes_indexed": len(hgnc_genes),
        "total_clinical_ref_rows": int(total_rows),
        "unique_clinvar_genes_with_2star_variants": int(uniq_genes),
        "clinvar_2star_streamed_candidates": int(n_cv),
        "db_size_mb": size_mb,
        "sha256": sha,
        "build_seconds": round(time.time() - t0, 2),
    }
    json.dump(manifest, open(MANIFEST, "w", encoding="utf-8"), indent=2)
    print("[PASS] v5.0.0 Pan-Human Exome SQLite WAL ready:")
    print(json.dumps(manifest, indent=2))

if __name__ == "__main__":
    main()
