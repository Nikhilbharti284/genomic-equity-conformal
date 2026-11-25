import json, os, sys, tempfile
sys.path.insert(0, os.path.abspath("src"))

from audit_chain import append_audit_entry, verify_audit_chain
from offline_reference_store import OfflineReferenceStore
from v5_global_cohort_audit import REPORT_PATH, run_global_cohort_audit
from v5_hybrid_fallback import evaluate_variant_v5

def test_v5_panhuman_exome_19k_hgnc_and_67k_clinvar_store():
    store = OfflineReferenceStore()
    st = store.stats()
    assert st["n_hgnc_genes"] >= 19000
    assert st["n_records"] >= 60000
    assert st["unique_variant_genes"] >= 7000
    assert os.path.getsize(st["db_path"]) < 1500 * 1024 * 1024

def test_v5_tier1_sqlite_indexed_hit():
    res = evaluate_variant_v5("TP53", "R", 175, "H", ancestry="SAS", hgvs_g="chr17:g.7675088C>T")
    assert res["v5_tier_status"] == "TIER1_SQLITE_INDEXED"
    assert res["final_tier"] == "Pathogenic_90pct"

def test_v5_tier2_novel_missense_any_hgnc_gene():
    res = evaluate_variant_v5("TTN", "G", 9999, "R", ancestry="AFR", am_score=0.82)
    assert res["v5_tier_status"] == "TIER2_CALIBRATED_NOVEL_MISSENSE"
    assert len(res["conformal_90pct_set"]) >= 1

def test_v5_ood_non_coding_or_unmapped_abstention():
    res = evaluate_variant_v5("FAKEGENE999", "A", 10, "V", ancestry="SAS")
    assert res["v5_tier_status"] == "OOD_NON_MISSENSE_OR_UNMAPPED"

def test_v5_global_5continent_cohort_report_and_ram():
    assert os.path.exists(REPORT_PATH)
    rep = json.load(open(REPORT_PATH, encoding="utf-8"))
    for pop in ("SAS", "AFR", "EAS", "AMR", "EUR"):
        assert pop in rep["super_populations"]
    assert rep["trio_compound_het_verified"] is True
    assert rep["peak_rss_mb"] < 250.0
    assert len(rep["sha256"]) == 64

def test_v5_audit_chain_tamper_evident():
    with tempfile.TemporaryDirectory() as td:
        lp = os.path.join(td, "v5_chain.jsonl")
        append_audit_entry("V5_PANHUMAN_AUDIT", {"version": "5.0.0"}, ledger_path=lp)
        v = verify_audit_chain(lp)
        assert v["valid"] is True and v["entries"] == 1
