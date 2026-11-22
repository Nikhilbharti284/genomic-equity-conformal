import os, sys, tempfile
sys.path.insert(0, os.path.abspath("src"))
from inference_engine import GenomicEquityEngine
from audit_vcf import write_annotated_vcf
from pgx_drug_safety import PGxDrugSafetySidecar
from mondrian_conformal import GLOBAL_7_ANCESTRIES, normalize_global_ancestry

def test_pan_global_7_ancestry_non_empty_conformal():
    eng = GenomicEquityEngine()
    for anc in (*GLOBAL_7_ANCESTRIES, "EUR", "ASJ", "FIN", "ULTRA_RARE"):
        card = eng.predict(gene="BRCA1", wt="C", pos=61, mt="G", llr=-5.0, enorm=10.0, ancestry=anc)
        assert len(card["conformal_90pct_set"]) >= 1
        assert 0.0 <= card["joint_triage_prob"] <= 1.0
    assert normalize_global_ancestry("ASJ") == "BOTTLENECK_FOUNDER"
    assert normalize_global_ancestry("FIN") == "BOTTLENECK_FOUNDER"

def test_cpic_pgx_drug_toxicity_alerts_and_vcf_tag():
    eng = GenomicEquityEngine()
    pgx = PGxDrugSafetySidecar()
    # NUDT15_R139C (EAS/SAS thiopurine toxicity)
    c1 = eng.predict(gene="NUDT15", wt="R", pos=139, mt="C", llr=-4.0, enorm=10.0, ancestry="EAS", af_eas=0.10)
    assert c1["pgx_alert"] is not None and "Thiopurines" in c1["pgx_alert"]["drug_class"]
    # DPYD_D949V (5-FU toxicity)
    c2 = eng.predict(gene="DPYD", wt="D", pos=949, mt="V", llr=-3.5, enorm=10.0, ancestry="SAS")
    assert c2["pgx_alert"] is not None and "Fluoropyrimidines" in c2["pgx_alert"]["drug_class"]
    # HLA-B*57:01 tag via rsid
    tag = pgx.lookup(rsid="rs2395029")
    assert tag is not None and "Abacavir" in tag["drug_class"]
    # VCF 4.2 GEC_PGX_ALERT output
    sample_vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr13\t48045719\trs116855232\tC\tT\t.\tPASS\tGENE=NUDT15\n"
    with tempfile.TemporaryDirectory() as td:
        inp, outp = os.path.join(td, "in.vcf"), os.path.join(td, "out.vcf")
        open(inp, "w").write(sample_vcf)
        write_annotated_vcf(inp, outp, {"chr13:g.48045719C>T": c1})
        txt = open(outp).read()
        assert "##INFO=<ID=GEC_PGX_ALERT," in txt
        assert "GEC_PGX_ALERT=Thiopurines" in txt

def test_v20_discovery_and_benchmark_artifacts_valid():
    import json
    for fp in ("artifacts/v20_global_7ancestry_report.json", "artifacts/v20_novel_founder_discoveries.json", "artifacts/v20_mavedb_rsa_report.json"):
        if os.path.exists(fp):
            d = json.load(open(fp, encoding="utf-8"))
            assert len(d.get("sha256", "")) == 64
