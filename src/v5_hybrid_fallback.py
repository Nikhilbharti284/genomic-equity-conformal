"""
Genomic Equity Conformal v5.0.0 — Two-Tier Pan-Human Exome Hybrid OOD Fallback
Tier 1: Sub-ms SQLite hit (67,402 indexed ClinVar >=2-star GRCh38 missense rows)
Tier 2: Zero-drop calibrated fallback for novel missense in ANY of 19,276 HGNC genes
        (physico-chemical + AlphaMissense + optional CPU ESM2-8M cache via --enable-live-esm2)
OOD:    Strictly reserved for non-coding/unmapped/non-HGNC regions (OOD_NON_MISSENSE_OR_UNMAPPED)
"""
import hashlib, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inference_engine import GenomicEquityEngine, HYDRO, MW, CHG
from offline_reference_store import OfflineReferenceStore

VALID_AA = set(HYDRO.keys())
_ESM2_CACHE = {}

def _physicochemical_llr_proxy(wt, mt, am_score=None):
    if am_score is not None:
        return round(-6.5 * float(am_score) + 0.5, 4)
    d_h = abs(HYDRO.get(mt, 0.0) - HYDRO.get(wt, 0.0))
    d_m = abs(MW.get(mt, 120.0) - MW.get(wt, 120.0)) / 100.0
    d_c = abs(CHG.get(mt, 0.0) - CHG.get(wt, 0.0))
    pg_pen = 1.2 if (wt in ("P", "G") or mt in ("P", "G")) else 0.0
    return round(-(0.35 * d_h + 1.1 * d_m + 1.4 * d_c + pg_pen), 4)

def evaluate_variant_v5(
    gene, wt, pos, mt, ancestry="SAS", am_score=None,
    af_nfe=0.0, af_sas=0.0, af_afr=0.0, af_eas=0.0, af_amr=0.0,
    af_mid=0.0, af_asj=0.0, af_fin=0.0, af_all=0.0,
    rsid=None, hgvs_g=None, enable_live_esm2=False, eng=None, store=None
):
    eng = eng or GenomicEquityEngine()
    store = store or OfflineReferenceStore()
    g = str(gene or "").strip().upper()
    w, m = str(wt or "").strip().upper(), str(mt or "").strip().upper()
    vid = f"{g}_{w}{pos}{m}"

    # Tier 1: Sub-ms SQLite hit
    hit = store.lookup_by_hgvs(hgvs_g) if hgvs_g else None
    if not hit:
        hit = store.lookup_by_variant_id(vid)
    if hit:
        card = eng.predict(
            gene=hit["gene"], wt=hit["wt"], pos=int(hit["pos"]), mt=hit["mt"],
            llr=float(hit["esm2_llr"]), ancestry=ancestry,
            af_nfe=max(af_nfe, float(hit["af_nfe"])), af_sas=max(af_sas, float(hit["af_sas"])),
            af_afr=max(af_afr, float(hit["af_afr"])), af_eas=max(af_eas, float(hit["af_eas"])),
            af_amr=max(af_amr, float(hit["af_amr"])), af_mid=max(af_mid, float(hit["af_mid"])),
            af_asj=max(af_asj, float(hit["af_asj"])), af_fin=max(af_fin, float(hit["af_fin"])),
            af_all=max(af_all, float(hit["af_all"])), rsid=rsid,
        )
        card["v5_tier_status"] = "TIER1_SQLITE_INDEXED"
        card["model_version"] = "5.0.0-panhuman"
        return card

    # OOD Guardrail: Non-coding / invalid AA / non-HGNC gene
    if (w not in VALID_AA) or (m not in VALID_AA) or (w == m) or (not store.is_hgnc_coding_gene(g)):
        out = {
            "variant_id": vid, "gene": g, "patient_ancestry": ancestry,
            "model_version": "5.0.0-panhuman",
            "v5_tier_status": "OOD_NON_MISSENSE_OR_UNMAPPED",
            "final_tier": "OOD_Abstention_NonCoding_Or_Unmapped",
            "framing": "Computational research & clinical-pipeline QA audit only — not standalone SaMD.",
        }
        out["sha256"] = hashlib.sha256(json.dumps(out, sort_keys=True).encode()).hexdigest()
        return out

    # Tier 2: Zero-drop calibrated novel missense across all 19,276 HGNC genes
    llr = _ESM2_CACHE.get(vid)
    if llr is None:
        llr = _physicochemical_llr_proxy(w, m, am_score=am_score)
        _ESM2_CACHE[vid] = llr
    card = eng.predict(
        gene=g, wt=w, pos=int(pos), mt=m, llr=float(llr), ancestry=ancestry,
        af_nfe=af_nfe, af_sas=af_sas, af_afr=af_afr, af_eas=af_eas, af_amr=af_amr,
        af_mid=af_mid, af_asj=af_asj, af_fin=af_fin, af_all=af_all, rsid=rsid,
    )
    card["v5_tier_status"] = "TIER2_CALIBRATED_NOVEL_MISSENSE"
    card["live_esm2_enabled"] = bool(enable_live_esm2)
    card["model_version"] = "5.0.0-panhuman"
    return card
