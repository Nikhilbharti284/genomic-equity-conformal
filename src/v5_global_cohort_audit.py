"""
Genomic Equity Conformal v5.0.0 — Multi-Continental Global Cohort Benchmark
Evaluates SAS, AFR, EAS, AMR, and EUR cohorts across Tier-1 SQLite + Tier-2 Pan-Human Exome
with BA1 balancing-selection protection, CPIC Level A PGx alerts, throughput & peak RSS RAM.
"""
import hashlib, json, os, resource, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inference_engine import GenomicEquityEngine
from offline_reference_store import OfflineReferenceStore
from trio_phasing import analyze_trio_and_compound_het
from v5_hybrid_fallback import evaluate_variant_v5

REPORT_PATH = "artifacts/v50_global_realworld_cohort_report.json"

COHORT_PANELS = {
    "SAS": [
        ("HBB", "E", 7, "V", 0.0006, 0.047, 0.0, 0.002, 0.0001, "chr11:g.5227002T>A"),
        ("NUDT15", "R", 139, "C", 0.068, 0.001, 0.104, 0.042, 0.002, "chr13:g.48045719C>T"),
        ("TP53", "R", 175, "H", 0.0, 0.0, 0.0, 0.0, 0.0, "chr17:g.7675088C>T"),
        ("TP53", "P", 72, "R", 0.50, 0.38, 0.58, 0.65, 0.73, "chr17:g.7676154G>C"),
        ("EGFR", "L", 858, "R", 0.0, 0.0, 0.0, 0.0, 0.0, None),
    ],
    "AFR": [
        ("HBB", "E", 7, "V", 0.0006, 0.047, 0.0, 0.002, 0.0001, "chr11:g.5227002T>A"),
        ("TP53", "R", 175, "H", 0.0, 0.0, 0.0, 0.0, 0.0, "chr17:g.7675088C>T"),
        ("TP53", "P", 72, "R", 0.50, 0.38, 0.58, 0.65, 0.73, "chr17:g.7676154G>C"),
        ("G6PD", "V", 68, "M", 0.0, 0.11, 0.0, 0.01, 0.0001, None),
    ],
    "EAS": [
        ("NUDT15", "R", 139, "C", 0.068, 0.001, 0.104, 0.042, 0.002, "chr13:g.48045719C>T"),
        ("TP53", "R", 175, "H", 0.0, 0.0, 0.0, 0.0, 0.0, "chr17:g.7675088C>T"),
        ("TP53", "P", 72, "R", 0.50, 0.38, 0.58, 0.65, 0.73, "chr17:g.7676154G>C"),
        ("GJB2", "V", 37, "I", 0.005, 0.0005, 0.082, 0.003, 0.0008, None),
    ],
    "AMR": [
        ("HBB", "E", 7, "V", 0.0006, 0.047, 0.0, 0.002, 0.0001, "chr11:g.5227002T>A"),
        ("NUDT15", "R", 139, "C", 0.068, 0.001, 0.104, 0.042, 0.002, "chr13:g.48045719C>T"),
        ("TP53", "R", 175, "H", 0.0, 0.0, 0.0, 0.0, 0.0, "chr17:g.7675088C>T"),
        ("NOTCH1", "G", 200, "R", 0.0, 0.0, 0.0, 0.0, 0.0, None),
    ],
    "EUR": [
        ("TP53", "R", 175, "H", 0.0, 0.0, 0.0, 0.0, 0.0, "chr17:g.7675088C>T"),
        ("TP53", "P", 72, "R", 0.50, 0.38, 0.58, 0.65, 0.73, "chr17:g.7676154G>C"),
        ("INS", "C", 96, "Y", 0.0, 0.0, 0.0, 0.0, 0.0, None),
    ],
}

def run_global_cohort_audit():
    t0 = time.perf_counter()
    eng = GenomicEquityEngine()
    store = OfflineReferenceStore()
    st = store.stats()
    cohorts = {}
    total_variants = 0

    for pop, variants in COHORT_PANELS.items():
        anc = "NFE" if pop == "EUR" else pop
        tier1, tier2, resolved, ba1_protected, pgx_alerts = 0, 0, 0, 0, 0
        for g, w, pos, m, af_sas, af_afr, af_eas, af_amr, af_nfe, hgvs in variants:
            card = evaluate_variant_v5(
                gene=g, wt=w, pos=pos, mt=m, ancestry=anc,
                af_sas=af_sas, af_afr=af_afr, af_eas=af_eas, af_amr=af_amr, af_nfe=af_nfe,
                hgvs_g=hgvs, eng=eng, store=store,
            )
            total_variants += 1
            if card["v5_tier_status"] == "TIER1_SQLITE_INDEXED":
                tier1 += 1
            elif card["v5_tier_status"] == "TIER2_CALIBRATED_NOVEL_MISSENSE":
                tier2 += 1
            if not str(card.get("final_tier", "")).startswith("VUS"):
                resolved += 1
            if card.get("safeguard_triggered"):
                ba1_protected += 1
            if card.get("pgx_alert"):
                pgx_alerts += 1
        cohorts[pop] = {
            "variants_evaluated": len(variants),
            "tier1_sqlite_hits": tier1,
            "tier2_panhuman_novel_hits": tier2,
            "definitive_non_vus_rate": round(resolved / max(len(variants), 1), 4),
            "ba1_safeguard_triggers": ba1_protected,
            "cpic_pgx_alerts": pgx_alerts,
        }

    # Verify Trio Phasing integration on compound-het Wilson disease control
    trio_res = analyze_trio_and_compound_het([
        {"variant_id": "ATP7B_H1069Q", "gene": "ATP7B", "moi": "AR", "proband_gt": "0/1", "mother_gt": "0/1", "father_gt": "0/0", "final_tier": "Pathogenic_90pct"},
        {"variant_id": "ATP7B_R778L",  "gene": "ATP7B", "moi": "AR", "proband_gt": "0/1", "mother_gt": "0/0", "father_gt": "0/1", "final_tier": "Pathogenic_90pct"},
    ])
    elapsed = max(time.perf_counter() - t0, 1e-6)
    peak_rss_mb = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 2)

    report = {
        "version": "5.0.0",
        "exome_store": st,
        "super_populations": cohorts,
        "trio_compound_het_verified": len(trio_res.get("compound_het_alerts", [])) >= 1,
        "total_variants_evaluated": total_variants,
        "elapsed_seconds": round(elapsed, 4),
        "throughput_variants_per_sec": round(total_variants / elapsed, 1),
        "peak_rss_mb": peak_rss_mb,
        "ram_under_250mb": peak_rss_mb < 250.0,
    }
    report["sha256"] = hashlib.sha256(json.dumps(report, sort_keys=True).encode()).hexdigest()
    os.makedirs("artifacts", exist_ok=True)
    json.dump(report, open(REPORT_PATH, "w", encoding="utf-8"), indent=2)
    print(json.dumps(report, indent=2))
    return report

if __name__ == "__main__":
    run_global_cohort_audit()
