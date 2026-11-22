import hashlib, json, math, os
import sys, joblib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pgx_drug_safety import PGxDrugSafetySidecar
from mondrian_conformal import normalize_global_ancestry
import numpy as np

HYDRO = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
MW    = {"A":89,"R":174,"N":132,"D":133,"C":121,"Q":146,"E":147,"G":75,"H":155,"I":131,"L":131,"K":146,"M":149,"F":165,"P":115,"S":105,"T":119,"W":204,"Y":181,"V":117}
CHG   = {"R":1,"K":1,"H":0.1,"D":-1,"E":-1}

def _logit(p):
    p = np.clip(p, 1e-6, 1.0 - 1e-6)
    return np.log(p / (1.0 - p))

def pejaver_acmg_tier(p_func_cal):
    """Maps AF-free calibrated probability (50/50 prior) to Pejaver et al. (2022) ACMG PP3/BP4 tiers."""
    if p_func_cal >= 0.97104: return ("PP3_Strong", +4)
    if p_func_cal >= 0.85273: return ("PP3_Moderate", +2)
    if p_func_cal >= 0.70642: return ("PP3_Supporting", +1)
    if p_func_cal <= 0.02896: return ("BP4_Strong", -4)
    if p_func_cal <= 0.14727: return ("BP4_Moderate", -2)
    if p_func_cal <= 0.29358: return ("BP4_Supporting", -1)
    return ("Indeterminate", 0)

def build_func_features(wt, mt, llr, enorm):
    d_h = HYDRO.get(mt, 0.0) - HYDRO.get(wt, 0.0)
    d_m = (MW.get(mt, 120.0) - MW.get(wt, 120.0)) / 100.0
    d_c = abs(CHG.get(mt, 0.0) - CHG.get(wt, 0.0))
    return [llr, enorm, d_h, abs(d_h), d_m, abs(d_m), d_c, float(wt in ("P","G") or mt in ("P","G"))]

def build_triage_features(wt, mt, llr, enorm, af_nfe=0.0, af_sas=0.0, af_afr=0.0, af_all=0.0):
    af_max = max(af_nfe, af_sas, af_afr, af_all)
    return build_func_features(wt, mt, llr, enorm) + [
        math.log10(af_nfe + 1e-7), math.log10(af_sas + 1e-7),
        math.log10(af_afr + 1e-7), math.log10(af_max + 1e-7),
    ]

class GenomicEquityEngine:
    """Centralized v1.0.0 inference, Pejaver PP3/BP4 & 3-Layer Balancing-Selection Engine."""
    def __init__(self, bundle_path="artifacts/equity_model_v10_clean.joblib", registry_path="data/ba1_exception_registry.json"):
        self.bundle = joblib.load(bundle_path)
        with open(registry_path) as f:
            self.registry = json.load(f)
        self.pgx = PGxDrugSafetySidecar()

    def predict(self, gene, wt, pos, mt, llr, enorm=10.0, ancestry="SAS", af_nfe=0.0, af_sas=0.0, af_afr=0.0, af_eas=0.0, af_amr=0.0, af_mid=0.0, af_asj=0.0, af_fin=0.0, af_all=0.0, rsid=None):
        gene = gene.upper()
        wt, mt = wt.upper(), mt.upper()
        vid = f"{gene}_{wt}{pos}{mt}"
        anc_norm = normalize_global_ancestry(ancestry)
        grpmax_af = max(af_nfe, af_sas, af_afr, af_eas, af_amr, af_mid, af_asj, af_fin, af_all)

        # 1. AF-Free Functional Head -> ClinGen Pejaver PP3/BP4 (zero AF double-counting!)
        xf = np.array([build_func_features(wt, mt, llr, enorm)], dtype=float)
        p_func_raw = float(self.bundle["func_clf"].predict_proba(self.bundle["func_scaler"].transform(xf))[0, 1])
        p_func_cal = float(self.bundle["func_platt"].predict_proba(_logit(np.array([p_func_raw])).reshape(-1, 1))[0, 1])
        acmg_code, acmg_pts = pejaver_acmg_tier(p_func_cal)

        # 2. Commercially Clean Joint Triage Head -> Nested OOF Platt + Non-Empty Mondrian Conformal
        xt = np.array([build_triage_features(wt, mt, llr, enorm, af_nfe, af_sas, af_afr, af_all)], dtype=float)
        xts = self.bundle["triage_scaler"].transform(xt)
        p_tr_raw = float(0.50 * self.bundle["triage_mlp"].predict_proba(xts)[0, 1] + 0.35 * self.bundle["triage_hgb"].predict_proba(xt)[0, 1] + 0.15 * self.bundle["triage_lr"].predict_proba(xts)[0, 1])
        p_tr_cal = float(self.bundle["triage_platt"].predict_proba(_logit(np.array([p_tr_raw])).reshape(-1, 1))[0, 1])
        ps = self.bundle["mondrian"].predict_sets(np.array([p_tr_cal]), np.array([anc_norm]))[0]
        tier = "Pathogenic_90pct" if ps == [1] else ("Benign_90pct" if ps == [0] else "VUS_Review_Required")

        # 3. 3-Layer Balancing-Selection & Recessive Safeguard
        safeguard_triggered = None
        exc = self.registry["exceptions"].get(vid)
        g_thr = self.registry["gene_thresholds"].get(gene, {"moi": "AD", "ba1": 0.05, "bs1": 0.01})
        if exc and exc.get("block_benign"):
            if tier == "Benign_90pct":
                tier = f"VUS_Review_Required (Layer-1 BA1 Exception: {exc['condition']})"
            safeguard_triggered = f"Layer-1 BA1 Exception ({exc['condition']})"
            if acmg_pts < 0:
                acmg_code, acmg_pts = ("Indeterminate (BA1 Exception Protected)", 0)
        elif grpmax_af >= max(g_thr["ba1"], 0.05):
            ps = [0]
            tier = "Benign_90pct"
            safeguard_triggered = f"Layer-2 ACMG BA1 Stand-Alone Benign (AF={grpmax_af:.4f} >= 0.05, non-exception)"
        elif tier == "Benign_90pct" and grpmax_af >= g_thr["bs1"] and (acmg_pts >= 2 or llr <= -2.0):
            tier = f"VUS_Review_Required (Layer-3 High-AF Functional Conflict; {g_thr['moi']} BS1={g_thr['bs1']})"
            safeguard_triggered = f"Layer-3 Functional Conflict (AF={grpmax_af:.4f} >= BS1={g_thr['bs1']}, LLR={llr:.2f})"

        out = {
            "variant_id": vid, "gene": gene, "patient_ancestry": ancestry,
            "model_version": self.bundle.get("version", "1.0.0-clean"),
            "commercial_license_clean": True,
            "esm2_llr": round(llr, 4),
            "af_free_functional_prob": round(p_func_cal, 4),
            "pejaver_acmg_code": acmg_code, "pejaver_points": acmg_pts,
            "joint_triage_prob": round(p_tr_cal, 4),
            "conformal_90pct_set": ["Benign" if c == 0 else "Pathogenic" for c in ps],
            "final_tier": tier,
            "ancestry_stratum": anc_norm,
            "safeguard_triggered": safeguard_triggered,
            "pgx_alert": self.pgx.lookup(variant_id=vid, rsid=rsid),
            "framing": "Computational research & clinical-pipeline QA audit only — not standalone SaMD."
        }
        out["sha256"] = hashlib.sha256(json.dumps(out, sort_keys=True).encode()).hexdigest()
        return out
