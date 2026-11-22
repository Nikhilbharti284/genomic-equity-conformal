import os, sys, tempfile
sys.path.insert(0, os.path.abspath("src"))
from inference_engine import GenomicEquityEngine
from audit_vcf import write_annotated_vcf

def test_engine_commercial_clean_and_hbb_ba1_protection():
    eng = GenomicEquityEngine()
    card = eng.predict(gene="HBB", wt="E", pos=7, mt="V", llr=-4.5, enorm=10.0,
                       ancestry="AFR", af_nfe=0.0, af_sas=0.000621, af_afr=0.046998, af_all=0.004)
    assert card["commercial_license_clean"] is True
    assert card["final_tier"] != "Benign_90pct", "Sickle Cell HBB_E7V must NEVER be called Benign!"
    assert card["safeguard_triggered"] is not None and "Layer-1 BA1 Exception" in card["safeguard_triggered"]
    assert card["pejaver_points"] >= 0

def test_tp53_pathogenic_and_benign_controls():
    eng = GenomicEquityEngine()
    path = eng.predict(gene="TP53", wt="R", pos=175, mt="H", llr=-6.2, enorm=11.0, ancestry="SAS")
    ben = eng.predict(gene="TP53", wt="P", pos=72, mt="R", llr=0.4, enorm=9.5,
                      ancestry="SAS", af_nfe=0.25, af_sas=0.500196, af_afr=0.379906, af_all=0.35)
    assert path["final_tier"] == "Pathogenic_90pct"
    assert ben["final_tier"] == "Benign_90pct"
    assert 0.0 <= path["joint_triage_prob"] <= 1.0
    assert 0.0 <= ben["joint_triage_prob"] <= 1.0

def test_zero_empty_prediction_sets_all_ancestries():
    eng = GenomicEquityEngine()
    for anc in ("EUR", "SAS", "AFR", "ULTRA_RARE"):
        for llr in (-8.0, -3.0, -0.5, 1.5):
            card = eng.predict(gene="BRCA1", wt="C", pos=61, mt="G", llr=llr, enorm=10.0, ancestry=anc)
            assert len(card["conformal_90pct_set"]) >= 1, f"Empty set detected for {anc}, llr={llr}"

def test_annotated_vcf_v42_output_format():
    eng = GenomicEquityEngine()
    card = eng.predict(gene="HBB", wt="E", pos=7, mt="V", llr=-4.5, enorm=10.0,
                       ancestry="AFR", af_sas=0.000621, af_afr=0.046998)
    sample_vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr11\t5227002\trs334\tT\tA\t.\tPASS\tGENE=HBB\n"
    with tempfile.TemporaryDirectory() as td:
        in_p, out_p = os.path.join(td, "in.vcf"), os.path.join(td, "out.vcf")
        open(in_p, "w").write(sample_vcf)
        write_annotated_vcf(in_p, out_p, {"chr11:g.5227002T>A": card})
        txt = open(out_p).read()
        assert "##INFO=<ID=GEC_PLATT," in txt
        assert "##INFO=<ID=GEC_TIER," in txt
        assert "GEC_SAFEGUARD=Layer-1_BA1_Exception" in txt
        assert "GEC_ANC=AFR" in txt
