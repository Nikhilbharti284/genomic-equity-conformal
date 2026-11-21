import os, sys
import numpy as np
sys.path.insert(0, os.path.abspath("src"))
from mondrian_conformal import AncestryMondrianConformal
from inference_engine import GenomicEquityEngine, pejaver_acmg_tier

def test_zero_empty_prediction_sets():
    """Verifies v1.0.0 fix for mondrian_conformal.py empty-set bug when q_g < 0.5."""
    mc = AncestryMondrianConformal(alpha=0.10)
    mc.q_by_group = {"SAS": 0.20}
    mc.q_global = 0.20
    sets = mc.predict_sets(np.array([0.45, 0.55, 0.50]), np.array(["SAS", "SAS", "SAS"]))
    assert all(len(s) >= 1 for s in sets), f"Empty set returned: {sets}"

def test_pejaver_acmg_thresholds():
    assert pejaver_acmg_tier(0.98) == ("PP3_Strong", 4)
    assert pejaver_acmg_tier(0.88) == ("PP3_Moderate", 2)
    assert pejaver_acmg_tier(0.50) == ("Indeterminate", 0)
    assert pejaver_acmg_tier(0.01) == ("BP4_Strong", -4)

def test_ba1_exception_registry_blocks_false_benign():
    if not os.path.exists("artifacts/equity_model_v10_clean.joblib"):
        return
    eng = GenomicEquityEngine()
    res = eng.predict(gene="HBB", wt="E", pos=7, mt="V", llr=-2.277, ancestry="AFR", af_afr=0.047)
    assert "Benign" not in res["final_tier"], f"Sickle Cell HBB_E7V falsely called Benign: {res}"
    assert res["safeguard_triggered"] is not None
    assert res["commercial_license_clean"] is True
