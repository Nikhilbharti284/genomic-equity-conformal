import csv, gzip, json, os, re, urllib.request
from collections import Counter, defaultdict
from resource_logger import ResourceLogger

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
AA3TO1 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
HGVS_P_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
TWO_STAR = {"reviewed by expert panel", "practice guideline", "criteria provided, multiple submitters, no conflicts"}

def fetch_uniprot_batch(genes):
    seqs, accs = {}, {}
    for g in sorted(genes):
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{g}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1"
        req = urllib.request.Request(url, headers={"User-Agent": "genomic-equity-conformal/1.0"})
        try:
            with urllib.request.urlopen(req, timeout=15) as resp:
                res = json.loads(resp.read().decode("utf-8")).get("results", [])
                if res:
                    seqs[g] = res[0]["sequence"]["value"]
                    accs[g] = res[0]["primaryAccession"]
        except Exception:
            pass
    return seqs, accs

def main(top_n_genes=120, max_total=10000):
    logger = ResourceLogger()
    print("1/3 Streaming ClinVar GRCh38 (>=2-star only, natural class prior, post-2023 flagged)...")
    req = urllib.request.Request(CLINVAR_URL, headers={"User-Agent": "genomic-equity-conformal/1.0"})
    by_gene = defaultdict(list)
    seen = set()
    with urllib.request.urlopen(req, timeout=180) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            header = gz.readline().decode("utf-8", errors="replace").rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(header)}
            for raw in gz:
                p = raw.decode("utf-8", errors="replace").rstrip("\n").split("\t")
                if len(p) <= idx["ReviewStatus"]: continue
                if p[idx["Assembly"]] != "GRCh38" or p[idx["Type"]] != "single nucleotide variant": continue
                if p[idx["ReviewStatus"]].lower() not in TWO_STAR: continue
                cs = p[idx["ClinicalSignificance"]].lower()
                if "pathogenic" in cs and "benign" not in cs and "conflicting" not in cs: label = 1
                elif "benign" in cs and "pathogenic" not in cs and "conflicting" not in cs: label = 0
                else: continue
                m = HGVS_P_RE.search(p[idx["Name"]])
                if not m: continue
                wt3, pos_s, mt3 = m.groups()
                if wt3 not in AA3TO1 or mt3 not in AA3TO1 or wt3 == mt3: continue
                gene = p[idx["GeneSymbol"]]
                if not gene or ";" in gene or "-" == gene: continue
                chrom, gpos, ref, alt = p[idx["Chromosome"]], p[idx["Start"]], p[idx["ReferenceAlleleVCF"]], p[idx["AlternateAlleleVCF"]]
                if ref not in "ACGT" or alt not in "ACGT": continue
                vid = f"{gene}_{AA3TO1[wt3]}{pos_s}{AA3TO1[mt3]}"
                if vid in seen: continue
                seen.add(vid)
                last_eval = p[idx["LastEvaluated"]] if "LastEvaluated" in idx and len(p) > idx["LastEvaluated"] else ""
                is_post_2023 = int(any(yr in last_eval for yr in ("2024", "2025", "2026")))
                by_gene[gene].append({
                    "variant_id": vid, "gene": gene, "pos": int(pos_s), "wt": AA3TO1[wt3], "mt": AA3TO1[mt3],
                    "label": label, "hgvs_g": f"chr{chrom}:g.{gpos}{ref}>{alt}", "is_post_2023": is_post_2023
                })

    # Pick top 120 genes with both pathogenic and benign variants (or highest total >=2-star count)
    ranked = sorted(by_gene.keys(), key=lambda g: (len({r["label"] for r in by_gene[g]}), len(by_gene[g])), reverse=True)[:top_n_genes]
    print(f"2/3 Fetching canonical UniProt sequences for top {len(ranked)} genes...")
    seqs, accs = fetch_uniprot_batch(ranked)

    final_rows = []
    for g in ranked:
        if g not in seqs: continue
        seq, acc = seqs[g], accs[g]
        gene_added = 0
        for r in by_gene[g]:
            if gene_added >= 85: break
            if 1 <= r["pos"] <= len(seq) and seq[r["pos"] - 1] == r["wt"]:
                r["uniprot_acc"] = acc
                r["seq"] = seq
                final_rows.append(r)
                gene_added += 1
                if len(final_rows) >= max_total:
                    break
        if len(final_rows) >= max_total:
            break

    os.makedirs("data", exist_ok=True)
    out_csv = "data/v1_clinvar_120genes.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["variant_id","gene","uniprot_acc","pos","wt","mt","label","hgvs_g","is_post_2023","seq"])
        w.writeheader()
        w.writerows(final_rows)
    sha = logger.sha256_file(out_csv)
    post23 = sum(r["is_post_2023"] for r in final_rows)
    print(f"3/3 Saved {len(final_rows)} >=2-star variants across {len({r['gene'] for r in final_rows})} genes (post-2023 holdout: {post23}) -> {out_csv}")
    print(f"    SHA-256: {sha}")
    logger.log("v1_curate_120genes", {"n_variants": len(final_rows), "n_genes": len({r["gene"] for r in final_rows}), "post_2023": post23, "sha256": sha})

if __name__ == "__main__":
    main()
