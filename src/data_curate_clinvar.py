import csv, gzip, json, os, re, urllib.request
from collections import defaultdict
from resource_logger import ResourceLogger

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
TARGET_GENES = {
    "BRCA1", "BRCA2", "TP53", "MLH1", "MSH2", "PTEN", "LDLR", "MYH7",
    "KCNQ1", "SCN5A", "CFTR", "PAH", "G6PD", "HBB", "GLA", "GAA",
    "APC", "ATM", "CHEK2", "PALB2"
}
AA3TO1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
    "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
    "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"
}
HGVS_P_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
GOOD_REVIEW = {
    "reviewed by expert panel",
    "practice guideline",
    "criteria provided, multiple submitters, no conflicts",
    "criteria provided, single submitter",
}

def fetch_uniprot_sequences(genes):
    seqs = {}
    for g in sorted(genes):
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{g}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1"
        req = urllib.request.Request(url, headers={"User-Agent": "genomic-equity-conformal/0.1"})
        try:
            with urllib.request.urlopen(req, timeout=20) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                results = data.get("results", [])
                if results:
                    seqs[g] = results[0]["sequence"]["value"]
        except Exception as e:
            print(f"[warn] UniProt fetch failed for {g}: {e}")
    return seqs

def fetch_gnomad_batch(hgvs_ids):
    out = {}
    url = "https://myvariant.info/v1/variant?fields=gnomad_exome.af,gnomad_genome.af"
    for i in range(0, len(hgvs_ids), 250):
        chunk = hgvs_ids[i:i+250]
        payload = json.dumps({"ids": chunk}).encode("utf-8")
        req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json", "User-Agent": "genomic-equity-conformal/0.1"})
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                rows = json.loads(resp.read().decode("utf-8"))
                for r in rows:
                    q = r.get("query")
                    if not q or r.get("notfound"):
                        continue
                    af_src = r.get("gnomad_exome", {}).get("af") or r.get("gnomad_genome", {}).get("af") or {}
                    if isinstance(af_src, dict):
                        out[q] = {
                            "af_nfe": float(af_src.get("af_nfe") or 0.0),
                            "af_sas": float(af_src.get("af_sas") or 0.0),
                            "af_afr": float(af_src.get("af_afr") or 0.0),
                            "af_all": float(af_src.get("af") or 0.0),
                        }
        except Exception as e:
            print(f"[warn] MyVariant batch {i} failed: {e}")
    return out

def main(max_per_gene=80):
    logger = ResourceLogger()
    print("1/3 Fetching canonical human protein sequences from UniProt...")
    seqs = fetch_uniprot_sequences(TARGET_GENES)
    print(f"    Loaded {len(seqs)}/{len(TARGET_GENES)} UniProt sequences.")

    print("2/3 Streaming ClinVar GRCh38 variant_summary.txt.gz (no raw file saved)...")
    req = urllib.request.Request(CLINVAR_URL, headers={"User-Agent": "genomic-equity-conformal/0.1"})
    counts = defaultdict(lambda: {0: 0, 1: 0})
    records = []
    seen = set()

    with urllib.request.urlopen(req, timeout=120) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            header = gz.readline().decode("utf-8", errors="replace").rstrip("\n").split("\t")
            idx = {col: i for i, col in enumerate(header)}
            for raw_line in gz:
                parts = raw_line.decode("utf-8", errors="replace").rstrip("\n").split("\t")
                if len(parts) <= idx["ReviewStatus"]:
                    continue
                if parts[idx["Assembly"]] != "GRCh38" or parts[idx["Type"]] != "single nucleotide variant":
                    continue
                gene = parts[idx["GeneSymbol"]]
                if gene not in seqs:
                    continue
                if parts[idx["ReviewStatus"]].lower() not in GOOD_REVIEW:
                    continue
                clinsig = parts[idx["ClinicalSignificance"]].lower()
                if "pathogenic" in clinsig and "benign" not in clinsig and "conflicting" not in clinsig:
                    label = 1
                elif "benign" in clinsig and "pathogenic" not in clinsig and "conflicting" not in clinsig:
                    label = 0
                else:
                    continue
                if counts[gene][label] >= max_per_gene // 2:
                    continue
                m = HGVS_P_RE.search(parts[idx["Name"]])
                if not m:
                    continue
                wt3, pos_s, mt3 = m.groups()
                if wt3 not in AA3TO1 or mt3 not in AA3TO1 or wt3 == mt3:
                    continue
                wt, mt, pos = AA3TO1[wt3], AA3TO1[mt3], int(pos_s)
                seq = seqs[gene]
                if pos < 1 or pos > len(seq) or seq[pos - 1] != wt:
                    continue
                vid = f"{gene}_{wt}{pos}{mt}"
                if vid in seen:
                    continue
                chrom = parts[idx["Chromosome"]]
                gpos = parts[idx["Start"]]
                ref = parts[idx["ReferenceAlleleVCF"]]
                alt = parts[idx["AlternateAlleleVCF"]]
                if ref not in "ACGT" or alt not in "ACGT":
                    continue
                hgvs_g = f"chr{chrom}:g.{gpos}{ref}>{alt}"
                seen.add(vid)
                counts[gene][label] += 1
                records.append({
                    "variant_id": vid, "gene": gene, "pos": pos, "wt": wt, "mt": mt,
                    "label": label, "hgvs_g": hgvs_g, "seq": seq
                })

    print(f"    Collected {len(records)} verified missense variants across {len(counts)} genes.")
    print("3/3 Fetching gnomAD ancestry AFs (EUR/SAS/AFR) via MyVariant.info...")
    afs = fetch_gnomad_batch([r["hgvs_g"] for r in records])

    os.makedirs("data", exist_ok=True)
    out_csv = "data/pilot_variants.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "variant_id","gene","pos","wt","mt","label","hgvs_g",
            "af_nfe","af_sas","af_afr","af_all","ancestry_group","seq"
        ])
        w.writeheader()
        for r in records:
            af = afs.get(r["hgvs_g"], {"af_nfe": 0.0, "af_sas": 0.0, "af_afr": 0.0, "af_all": 0.0})
            # Assign dominant observed ancestry cohort for Mondrian subgroup evaluation
            if af["af_sas"] > 0 and af["af_sas"] >= max(af["af_nfe"], af["af_afr"]):
                grp = "SAS"
            elif af["af_afr"] > 0 and af["af_afr"] >= max(af["af_nfe"], af["af_sas"]):
                grp = "AFR"
            elif af["af_nfe"] > 0:
                grp = "EUR"
            else:
                grp = "ULTRA_RARE"
            r.update(af)
            r["ancestry_group"] = grp
            w.writerow(r)

    sha = logger.sha256_file(out_csv)
    stats = logger.log("curate_clinvar", {"n_variants": len(records), "sha256": sha})
    print(f"Saved {out_csv} (SHA-256: {sha[:16]}...)")
    print("Run summary:", stats)

if __name__ == "__main__":
    main()
