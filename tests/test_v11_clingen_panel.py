import json, os, re, sys
sys.path.insert(0, os.path.abspath("src"))
from inference_engine import GenomicEquityEngine

def test_all_registered_clingen_ba1_exceptions_blocked_from_benign():
    eng = GenomicEquityEngine()
    reg = json.load(open("data/ba1_exception_registry.json", encoding="utf-8"))
    assert len(reg["exceptions"]) >= 8
    for vid, meta in reg["exceptions"].items():
        if not meta.get("block_benign"):
            continue
        m = re.match(r"^([A-Z0-9]+)_([A-Z])(\d+)([A-Z])$", vid)
        assert m, f"Invalid variant_id format in registry: {vid}"
        gene, wt, pos, mt = m.group(1), m.group(2), int(m.group(3)), m.group(4)
        card = eng.predict(gene=gene, wt=wt, pos=pos, mt=mt, llr=-1.0, enorm=10.0,
                           ancestry="SAS", af_sas=0.06, af_afr=0.06, af_all=0.06)
        assert card["final_tier"] != "Benign_90pct", f"BA1 exception {vid} leaked as Benign!"
        assert card["safeguard_triggered"] and "Layer-1 BA1 Exception" in card["safeguard_triggered"]

def test_layer2_standalone_ba1_vs_layer3_grayzone_conflict():
    eng = GenomicEquityEngine()
    # Layer-2: Non-exception variant with AF >= 0.05 -> Stand-Alone Benign even with harsh LLR
    l2 = eng.predict(gene="TP53", wt="P", pos=72, mt="R", llr=-3.3, enorm=9.5,
                     ancestry="SAS", af_sas=0.50, af_afr=0.38, af_all=0.73)
    assert l2["final_tier"] == "Benign_90pct"
    assert "Layer-2 ACMG BA1" in (l2["safeguard_triggered"] or "")

    # Layer-3: Non-exception variant in BS1..0.05 gray zone (AF=0.02) with strong deleterious LLR (-5.5)
    l3 = eng.predict(gene="BRCA1", wt="M", pos=1775, mt="R", llr=-5.5, enorm=10.5,
                     ancestry="AFR", af_afr=0.02, af_all=0.02)
    assert l3["final_tier"] != "Benign_90pct"
