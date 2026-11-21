import argparse, csv, io, json, math, os, re, urllib.request
import numpy as np
import torch
from scipy.stats import spearmanr
from transformers import AutoTokenizer, AutoModelForMaskedLM
from alphafold_v6_features import fetch_uniprot_acc, fetch_afdb_v6_plddt
from resource_logger import ResourceLogger

AA3TO1 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
HGVS_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")

SUITE = [
    ("TP53",  "urn:mavedb:00000068-0-1"),
    ("BRCA1", "urn:mavedb:00001222-b-2"),
    ("PTEN",  "urn:mavedb:00000054-a-1"),
    ("MSH2",  "urn:mavedb:00000050-a-1"),
    ("LDLR",  "urn:mavedb:00001269-a-1"),
]

def fetch_seq(gene):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1"
    req = urllib.request.Request(url, headers={"User-Agent": "genomic-equity-conformal/1.0"})
    with urllib.request.urlopen(req, timeout=20) as resp:
        return json.loads(resp.read().decode("utf-8"))["results"][0]["sequence"]["value"]

def crop_window(seq, pos_1idx, max_len=512):
    if len(seq) <= max_len: return seq, pos_1idx
    start = max(0, (pos_1idx - 1) - max_len // 2)
    end = min(len(seq), start + max_len)
    return seq[max(0, end - max_len):end], pos_1idx - max(0, end - max_len)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-per-gene", type=int, default=150, help="Use 150 for fast pilot, 1000 for full benchmark")
    args = ap.parse_args()

    logger = ResourceLogger()
    tok = AutoTokenizer.from_pretrained("facebook/esm2_t12_35M_UR50D")
    mod = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t12_35M_UR50D")
    mod.eval()

    results = {}
    wins = 0
    for gene, urn in SUITE:
        print(f"\n--- Evaluating {gene} ({urn}) ---")
        try:
            url = f"https://api.mavedb.org/api/v1/score-sets/{urn}/scores"
            req = urllib.request.Request(url, headers={"User-Agent": "genomic-equity-conformal/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                text = resp.read().decode("utf-8", errors="replace")
            seq = fetch_seq(gene)
            acc = fetch_uniprot_acc(gene)
            plddt_map, _ = fetch_afdb_v6_plddt(acc)

            parsed = []
            for row in csv.DictReader(io.StringIO(text)):
                m = HGVS_RE.match((row.get("hgvs_pro") or "").strip())
                if not m: continue
                wt3, pos_s, mt3 = m.groups()
                if wt3 not in AA3TO1 or mt3 not in AA3TO1 or wt3 == mt3: continue
                wt, mt, pos = AA3TO1[wt3], AA3TO1[mt3], int(pos_s)
                if pos < 1 or pos > len(seq) or seq[pos - 1] != wt: continue
                try: sc = float(row["score"])
                except Exception: continue
                if math.isnan(sc): continue
                parsed.append((wt, pos, mt, sc))

            rng = np.random.default_rng(42)
            idx = rng.choice(len(parsed), size=min(args.max_per_gene, len(parsed)), replace=False)
            sample = [parsed[i] for i in idx]

            pos_cache = {}
            raw_s, struct_s, lab_s = [], [], []
            for wt, pos, mt, sc in sample:
                if pos not in pos_cache:
                    wseq, wpos = crop_window(seq, pos, max_len=512)
                    inp = tok(wseq, return_tensors="pt")
                    inp["input_ids"][0, wpos] = tok.mask_token_id
                    with torch.no_grad():
                        out = mod(**inp)
                        pos_cache[pos] = torch.log_softmax(out.logits[0, wpos], dim=-1)
                lp = pos_cache[pos]
                llr = float(lp[tok.convert_tokens_to_ids(mt)] - lp[tok.convert_tokens_to_ids(wt)])
                plddt = plddt_map.get(pos, 70.0) / 100.0
                # Zero-circularity AlphaFold DB v6 structural attenuation (ordered core amplified, disordered loop damped)
                s_struct = llr * (0.35 + 0.65 * plddt)
                raw_s.append(llr); struct_s.append(s_struct); lab_s.append(sc)

            r_raw, p_raw = spearmanr(raw_s, lab_s)
            r_str, p_str = spearmanr(struct_s, lab_s)
            beat = abs(float(r_str)) > abs(float(r_raw))
            if beat: wins += 1
            results[gene] = {
                "urn": urn, "n": len(sample), "uniprot": acc,
                "esm2_35m_raw_rho": round(abs(float(r_raw)), 4),
                "esm2_35m_afdb_v6_rho": round(abs(float(r_str)), 4),
                "delta_rho": round(abs(float(r_str)) - abs(float(r_raw)), 4),
                "p_val": float(f"{p_str:.2e}"), "beats_raw_35m": beat
            }
            print("   ", results[gene])
        except Exception as e:
            print(f"    [warn] {gene} failed: {e}")

    summary = {"version": "1.0.0", "model": "esm2_t12_35M_UR50D + AlphaFold DB v6 pLDDT", "genes_beating_raw_35m": f"{wins}/{len(results)}", "per_gene": results}
    out_json = "artifacts/v1_mavedb_5gene_report.json"
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)
    summary["sha256"] = logger.sha256_file(out_json)
    logger.log("v1_mavedb_5gene_suite", {"wins": summary["genes_beating_raw_35m"], "sha256": summary["sha256"]})
    print("\nFINAL 5-GENE MAVEDB SUMMARY:\n", json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
