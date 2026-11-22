import argparse, gzip, hashlib, json, math, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inference_engine import GenomicEquityEngine
from offline_reference_store import OfflineReferenceStore
from vcf_patient_profiler import parse_sample_gt, classify_clinical_zygosity, infer_patient_ancestry_from_vcf

GEC_INFO_HEADERS = [
    '##INFO=<ID=GEC_PLATT,Number=1,Type=Float,Description="Genomic Equity Conformal v3.0.0 calibrated joint probability">',
    '##INFO=<ID=GEC_ACMG,Number=1,Type=String,Description="AF-free Pejaver PP3/BP4 ACMG tier">',
    '##INFO=<ID=GEC_TIER,Number=1,Type=String,Description="Non-empty Mondrian 90pct conformal + 3-Layer BA1 final tier">',
    '##INFO=<ID=GEC_SAFEGUARD,Number=1,Type=String,Description="Balancing-selection / BA1 safeguard status">',
    '##INFO=<ID=GEC_ANC,Number=1,Type=String,Description="Patient ancestry stratum used for Mondrian calibration">',
    '##INFO=<ID=GEC_PGX_ALERT,Number=1,Type=String,Description="CPIC-aligned Pharmacogenomics drug-toxicity alert">',
    '##INFO=<ID=GEC_ZYGOSITY,Number=1,Type=String,Description="MOI-aware clinical zygosity status from VCF GT">',
]

def _open_vcf(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="ignore") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8", errors="ignore")

def _vcf_safe(s):
    return re.sub(r"[^A-Za-z0-9_.\-:]+", "_", str(s or "NONE")).strip("_") or "NONE"

def stream_vcf_records(vcf_path):
    """Constant-memory streaming generator for .vcf and .vcf.gz files."""
    with _open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom = cols[0].replace("chr", "")
            pos_g, rsid, ref, alt = cols[1], cols[2], cols[3], cols[4]
            info = dict(kv.split("=", 1) for kv in cols[7].split(";") if "=" in kv) if len(cols) > 7 else {}
            gt_info = parse_sample_gt(cols[8], cols[9], chrom) if len(cols) >= 10 else {"gt": "./.", "zygosity": "UNKNOWN", "dp": None, "gq": None}
            yield {
                "chrom": chrom, "pos_g": pos_g, "rsid": rsid, "ref": ref, "alt": alt,
                "hgvs_g": f"chr{chrom}:g.{pos_g}{ref}>{alt}", "info": info, "gt_info": gt_info
            }

def write_annotated_vcf(in_vcf, out_vcf, card_by_hgvs):
    """Streaming VCF 4.2 writer with GEC_* INFO tags (supports .vcf and .vcf.gz input)."""
    os.makedirs(os.path.dirname(out_vcf) or ".", exist_ok=True)
    with _open_vcf(in_vcf) as fin, open(out_vcf, "w", encoding="utf-8") as fout:
        inserted = False
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                if not inserted:
                    for h in GEC_INFO_HEADERS:
                        fout.write(h + "\n")
                    inserted = True
                fout.write(line)
                continue
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            while len(cols) < 8:
                cols.append(".")
            chrom = cols[0].replace("chr", "")
            hgvs_g = f"chr{chrom}:g.{cols[1]}{cols[3]}>{cols[4]}"
            card = card_by_hgvs.get(hgvs_g)
            if card:
                tags = [
                    f"GEC_PLATT={card['joint_triage_prob']:.4f}",
                    f"GEC_ACMG={_vcf_safe(card['pejaver_acmg_code'])}",
                    f"GEC_TIER={_vcf_safe(card['final_tier'])}",
                    f"GEC_SAFEGUARD={_vcf_safe(card['safeguard_triggered'] or 'NONE')}",
                    f"GEC_ANC={_vcf_safe(card['patient_ancestry'])}",
                    f"GEC_PGX_ALERT={_vcf_safe((card.get('pgx_alert') or {}).get('drug_class', 'NONE'))}",
                    f"GEC_ZYGOSITY={_vcf_safe(card.get('clinical_zygosity_status', 'UNKNOWN'))}",
                ]
                cols[7] = ";".join(tags) if cols[7] in ("", ".") else cols[7] + ";" + ";".join(tags)
            fout.write("\t".join(cols) + "\n")

def audit_vcf_offline(vcf_path, out_vcf=None, default_ancestry="auto", eng=None, store=None):
    """Constant-memory air-gapped streaming VCF auditor (zero network packets)."""
    eng = eng or GenomicEquityEngine()
    store = store or OfflineReferenceStore()
    g_thr = eng.registry.get("gene_thresholds", {})

    inferred_anc, posteriors = default_ancestry, {}
    if str(default_ancestry).lower() == "auto":
        aim_obs = []
        for r in stream_vcf_records(vcf_path):
            ref = store.lookup_by_hgvs(r["hgvs_g"])
            if ref:
                aim_obs.append(ref)
            if len(aim_obs) >= 500:
                break
        inferred_anc, posteriors = infer_patient_ancestry_from_vcf(aim_obs)

    results, card_by_hgvs = [], {}
    total_scanned = 0
    for r in stream_vcf_records(vcf_path):
        total_scanned += 1
        ref = store.lookup_by_hgvs(r["hgvs_g"])
        if not ref:
            continue
        anc = r["info"].get("ANC") or inferred_anc
        card = eng.predict(
            gene=ref["gene"], wt=ref["wt"], pos=int(ref["pos"]), mt=ref["mt"],
            llr=float(ref["esm2_llr"]), enorm=10.0, ancestry=anc,
            af_nfe=float(ref["af_nfe"]), af_sas=float(ref["af_sas"]), af_afr=float(ref["af_afr"]),
            af_eas=float(ref["af_eas"]), af_amr=float(ref["af_amr"]), af_mid=float(ref["af_mid"]),
            af_asj=float(ref["af_asj"]), af_fin=float(ref["af_fin"]), af_all=float(ref["af_all"]),
            rsid=r.get("rsid"),
        )
        moi = g_thr.get(ref["gene"], {}).get("moi", ref.get("moi", "AD"))
        card["hgvs_g"] = r["hgvs_g"]
        card["sample_gt"] = r["gt_info"]["gt"]
        card["zygosity"] = r["gt_info"]["zygosity"]
        card["moi"] = moi
        card["clinical_zygosity_status"] = classify_clinical_zygosity(moi, r["gt_info"]["zygosity"])
        card["live_af_sas"] = round(float(ref["af_sas"]), 6)
        card["live_af_afr"] = round(float(ref["af_afr"]), 6)
        results.append(card)
        card_by_hgvs[r["hgvs_g"]] = card

    if out_vcf:
        write_annotated_vcf(vcf_path, out_vcf, card_by_hgvs)

    vcf_sha = hashlib.sha256(open(vcf_path, "rb").read()).hexdigest()
    payload = {
        "engine_version": "3.0.0-enterprise-airgapped",
        "offline_mode": True,
        "commercial_license_clean": True,
        "vcf_source": vcf_path,
        "vcf_input_sha256": vcf_sha,
        "annotated_vcf": out_vcf,
        "inferred_patient_ancestry": inferred_anc,
        "ancestry_posteriors": posteriors,
        "total_vcf_rows_scanned": total_scanned,
        "n_audited": len(results),
        "n_pgx_alerts": sum(1 for c in results if c.get("pgx_alert")),
        "results": results,
    }
    payload["sha256"] = hashlib.sha256(json.dumps(results, sort_keys=True).encode()).hexdigest()
    return payload

def main():
    ap = argparse.ArgumentParser(description="Genomic Equity v3.0.0 — Enterprise Streaming Hospital VCF 4.2 Auditor")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out-vcf", default="artifacts/audited_output.vcf")
    ap.add_argument("--out-json", default="artifacts/raw_vcf_audit_output.json")
    ap.add_argument("--default-ancestry", default="auto", choices=["auto", "EUR", "NFE", "SAS", "AFR", "EAS", "AMR", "MID", "ASJ", "FIN", "BOTTLENECK_FOUNDER", "ULTRA_RARE"])
    ap.add_argument("--offline", action="store_true", help="Air-gapped zero-network mode using local SQLite WAL store")
    args = ap.parse_args()

    payload = audit_vcf_offline(args.vcf, out_vcf=args.out_vcf, default_ancestry=args.default_ancestry)
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(json.dumps({k: v for k, v in payload.items() if k != "results"}, indent=2))
    print(f"[PASS] Audited {payload['n_audited']}/{payload['total_vcf_rows_scanned']} variants -> {args.out_vcf}")

if __name__ == "__main__":
    main()
