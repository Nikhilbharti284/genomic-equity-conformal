"""
Genomic Equity Conformal v4.0.0 — Global HL7 FHIR R4 Genomics Exporter (MIT License)
Converts v3.0.0/v4.0.0 clinical QA audit results into standard HL7 FHIR R4
Bundle (DiagnosticReport + Genomics Observation resources using LOINC codes).
Compatible worldwide with Epic, Cerner, Bahmni, OpenMRS, UK NHS, EU & ABDM LIMS.
"""
import hashlib
import json
import time


def build_fhir_r4_bundle(audit_payload: dict, patient_ref: str = "Patient/anonymous-local") -> dict:
    ts = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    results = audit_payload.get("results", [])
    bundle_id = "gec-fhir-" + hashlib.sha256(
        json.dumps(results, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()[:16]

    obs_refs = []
    entries = []

    for idx, r in enumerate(results, 1):
        obs_id = f"obs-var-{idx}"
        obs_refs.append({"reference": f"Observation/{obs_id}"})
        vid = r.get("variant_id") or r.get("variant") or f"VAR_{idx}"
        gene = r.get("gene") or (vid.split("_")[0] if "_" in vid else "UNKNOWN")
        tier = r.get("final_tier") or r.get("conformal_tier") or "VUS_Review_Required"
        zyg = r.get("clinical_zygosity_status") or r.get("zygosity") or "Unknown"
        trio = r.get("compound_het_status") or "None"
        pgx = r.get("gec_pgx_alert") or "None"
        ood = r.get("ood_status") or "IN_PANEL"
        prob = r.get("platt_calibrated_prob", r.get("platt_prob"))

        components = [
            {"code": {"coding": [{"system": "http://loinc.org", "code": "48018-6", "display": "Gene studied"}]}, "valueString": gene},
            {"code": {"coding": [{"system": "http://loinc.org", "code": "53034-5", "display": "Allelic state"}]}, "valueString": zyg},
            {"code": {"coding": [{"system": "http://loinc.org", "code": "53037-8", "display": "Genetic variation clinical significance"}]}, "valueString": tier},
            {"code": {"text": "GEC_Trio_Compound_Het"}, "valueString": trio},
            {"code": {"text": "GEC_CPIC_PGx_Alert"}, "valueString": pgx},
            {"code": {"text": "GEC_OOD_Panel_Status"}, "valueString": ood},
        ]
        if prob is not None:
            components.append({"code": {"text": "GEC_Platt_Calibrated_Prob"}, "valueQuantity": {"value": round(float(prob), 4)}})

        entries.append({
            "fullUrl": f"urn:uuid:{obs_id}",
            "resource": {
                "resourceType": "Observation",
                "id": obs_id,
                "status": "final",
                "category": [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category", "code": "laboratory"}]}],
                "code": {"coding": [{"system": "http://loinc.org", "code": "69548-6", "display": "Genetic variant assessment"}]},
                "subject": {"reference": patient_ref},
                "effectiveDateTime": ts,
                "valueCodeableConcept": {"text": vid},
                "component": components,
            },
        })

    diag_report = {
        "resourceType": "DiagnosticReport",
        "id": f"dr-{bundle_id}",
        "status": "final",
        "category": [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074", "code": "GE", "display": "Genetics"}]}],
        "code": {"coding": [{"system": "http://loinc.org", "code": "51969-4", "display": "Genetic analysis report"}]},
        "subject": {"reference": patient_ref},
        "issued": ts,
        "result": obs_refs,
        "conclusion": (
            f"Genomic Equity Conformal v4.0.0 Air-Gapped QA Audit ({len(results)} variants). "
            "Computational research & clinical-pipeline QA sidecar only — not standalone SaMD."
        ),
    }

    entries.insert(0, {"fullUrl": f"urn:uuid:dr-{bundle_id}", "resource": diag_report})
    return {
        "resourceType": "Bundle",
        "id": bundle_id,
        "type": "collection",
        "timestamp": ts,
        "total": len(entries),
        "entry": entries,
    }


if __name__ == "__main__":
    sample = {"results": [{"variant_id": "HBB_E7V", "gene": "HBB", "clinical_zygosity_status": "AR_Carrier_Only", "final_tier": "VUS_Review_Required (Layer-1 BA1 Exception)", "platt_calibrated_prob": 0.0195}]}
    out = build_fhir_r4_bundle(sample)
    print("FHIR Bundle OK:", out["id"], "| entries:", out["total"])
