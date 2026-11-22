import os
import sys

sys.path.insert(0, os.path.abspath("src"))
from audit_chain import append_audit_entry, verify_audit_chain
from fhir_exporter import build_fhir_r4_bundle
from lan_hospital_bridge import annotate_ood
from trio_phasing import analyze_trio_and_compound_het


def test_trio_trans_compound_het():
    recs = [
        {"variant_id": "ATP7B_H1069Q", "gene": "ATP7B", "moi": "AR", "proband_gt": "0/1", "mother_gt": "0/1", "father_gt": "0/0"},
        {"variant_id": "ATP7B_R778L",  "gene": "ATP7B", "moi": "AR", "proband_gt": "0/1", "mother_gt": "0/0", "father_gt": "0/1"},
    ]
    res = analyze_trio_and_compound_het(recs)
    assert res["per_variant"]["ATP7B_H1069Q"]["compound_het_status"] == "AR_Compound_Het_Trans_Affected"
    assert len(res["compound_het_alerts"]) == 1


def test_trio_cis_carrier_only():
    recs = [
        {"variant_id": "ABCA4_V1", "gene": "ABCA4", "moi": "AR", "proband_gt": "0/1", "mother_gt": "0/1", "father_gt": "0/0"},
        {"variant_id": "ABCA4_V2", "gene": "ABCA4", "moi": "AR", "proband_gt": "0/1", "mother_gt": "0/1", "father_gt": "0/0"},
    ]
    res = analyze_trio_and_compound_het(recs)
    assert res["per_variant"]["ABCA4_V1"]["compound_het_status"] == "AR_Cis_Carrier_Only"


def test_trio_phased_trans_without_parents():
    recs = [
        {"variant_id": "PAH_V1", "gene": "PAH", "moi": "AR", "proband_gt": "0|1"},
        {"variant_id": "PAH_V2", "gene": "PAH", "moi": "AR", "proband_gt": "1|0"},
    ]
    res = analyze_trio_and_compound_het(recs)
    assert res["per_variant"]["PAH_V1"]["compound_het_status"] == "AR_Compound_Het_Trans_Affected"


def test_trio_de_novo_detection():
    recs = [{"variant_id": "TP53_R175H", "gene": "TP53", "moi": "AD", "proband_gt": "0/1", "mother_gt": "0/0", "father_gt": "0/0"}]
    res = analyze_trio_and_compound_het(recs)
    assert res["per_variant"]["TP53_R175H"]["inheritance"] == "De_Novo_Candidate"


def test_fhir_r4_bundle_loinc_structure():
    bundle = build_fhir_r4_bundle({"results": [{"variant_id": "HBB_E7V", "gene": "HBB", "final_tier": "VUS_Review_Required", "platt_calibrated_prob": 0.0195}]})
    assert bundle["resourceType"] == "Bundle"
    assert bundle["entry"][0]["resource"]["resourceType"] == "DiagnosticReport"
    assert bundle["entry"][0]["resource"]["code"]["coding"][0]["code"] == "51969-4"
    assert bundle["entry"][1]["resource"]["code"]["coding"][0]["code"] == "69548-6"


def test_audit_chain_detects_single_byte_tamper(tmp_path):
    log = str(tmp_path / "chain.jsonl")
    append_audit_entry("VCF_AUDIT", {"sample": "P1", "n": 2}, log)
    append_audit_entry("FHIR_EXPORT", {"sample": "P1", "ok": True}, log)
    assert verify_audit_chain(log)["valid"] is True
    # Tamper 1 character in line 1
    txt = open(log, encoding="utf-8").read().replace('"P1"', '"P9"', 1)
    open(log, "w", encoding="utf-8").write(txt)
    ver = verify_audit_chain(log)
    assert ver["valid"] is False
    assert ver["broken_at_seq"] == 1


def test_ood_panel_abstention():
    res = annotate_ood([{"variant_id": "NOVELXYZ_A1V", "gene": "NOVELXYZ"}])
    assert res[0]["ood_status"] in ("IN_OFFLINE_PANEL", "OOD_NOT_IN_OFFLINE_PANEL")
