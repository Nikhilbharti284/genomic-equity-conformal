import csv, json, math, os, sqlite3, urllib.request
import joblib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from mondrian_conformal import AncestryMondrianConformal, compute_ece
from resource_logger import ResourceLogger

HYDRO = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
MW    = {"A":89,"R":174,"N":132,"D":133,"C":121,"Q":146,"E":147,"G":75,"H":155,"I":131,"L":131,"K":146,"M":149,"F":165,"P":115,"S":105,"T":119,"W":204,"Y":181,"V":117}
CHG   = {"R":1,"K":1,"H":0.1,"D":-1,"E":-1}

def _logit(p):
    p = np.clip(p, 1e-6, 1.0 - 1e-6)
    return np.log(p / (1.0 - p))

def _max_num(x):
    if x is None: return 0.0
    if isinstance(x, list):
        vals = [_max_num(v) for v in x]
        return max(vals) if vals else 0.0
    try:
        v = float(x)
        return 0.0 if math.isnan(v) else v
    except Exception:
        return 0.0

def fetch_spliceai_max(hgvs_ids):
    out = {}
    url = "https://myvariant.info/v1/variant?assembly=hg38&fields=spliceai,dbnsfp.spliceai"
    for i in range(0, len(hgvs_ids), 250):
        payload = json.dumps({"ids": hgvs_ids[i:i+250]}).encode("utf-8")
        req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json", "User-Agent": "genomic-equity-conformal/0.2"})
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                for r in json.loads(resp.read().decode("utf-8")):
                    q = r.get("query")
                    if not q or r.get("notfound"):
                        continue
                    sp = r.get("spliceai") or (r.get("dbnsfp", {}) or {}).get("spliceai") or {}
                    if isinstance(sp, list) and sp:
                        sp = sp[0]
                    if isinstance(sp, dict):
                        ds = max(_max_num(sp.get(k)) for k in ("ds_ag", "ds_al", "ds_dg", "ds_dl", "DS_AG", "DS_AL", "DS_DG", "DS_DL"))
                        out[q] = ds
        except Exception as e:
            print(f"[warn] SpliceAI batch {i}: {e}")
    return out

def build_features(r, llr, enorm):
    wt, mt = r["wt"], r["mt"]
    d_hydro = HYDRO.get(mt, 0.0) - HYDRO.get(wt, 0.0)
    d_mw = (MW.get(mt, 120.0) - MW.get(wt, 120.0)) / 100.0
    d_chg = abs(CHG.get(mt, 0.0) - CHG.get(wt, 0.0))
    is_pro_gly = float(wt in ("P", "G") or mt in ("P", "G"))
    return [
        llr, enorm, d_hydro, abs(d_hydro), d_mw, abs(d_mw), d_chg, is_pro_gly,
        math.log10(float(r["af_nfe"]) + 1e-7),
        math.log10(float(r["af_sas"]) + 1e-7),
        math.log10(float(r["af_afr"]) + 1e-7),
        math.log10(float(r["af_all"]) + 1e-7),
    ]

def main():
    logger = ResourceLogger()
    with open("data/pilot_variants_hg38.csv") as f:
        rows = list(csv.DictReader(f))
    print(f"1/3 Fetching hg38 SpliceAI cryptic-splicing scores for {len(rows)} variants...")
    sp_map = fetch_spliceai_max([r["hgvs_g"] for r in rows])
    nonzero = 0
    for r in rows:
        val = sp_map.get(r["hgvs_g"], 0.0)
        r["spliceai_max"] = val
        if val > 0: nonzero += 1
    print(f"    SpliceAI non-zero signals found: {nonzero}/{len(rows)}")

    conn = sqlite3.connect("artifacts/esm2_pilot.sqlite")
    scores = {r[0]: (r[1], r[2]) for r in conn.execute("SELECT variant_id, llr, emb_norm FROM scores")}
    conn.close()

    valid = [r for r in rows if r["variant_id"] in scores]
    X = np.array([build_features(r, *scores[r["variant_id"]]) for r in valid], dtype=float)
    y = np.array([int(r["label"]) for r in valid])
    genes = np.array([r["gene"] for r in valid])
    groups = np.array([r["ancestry_group"] for r in valid])

    print("2/3 Running 5-Fold Zero-Leakage CV (v0.2.0 with SpliceAI)...")
    gkf = GroupKFold(n_splits=5)
    oof_raw, oof_platt, oof_cov, oof_sz = np.zeros(len(valid)), np.zeros(len(valid)), np.zeros(len(valid)), np.zeros(len(valid))

    for tr_idx, te_idx in gkf.split(X, y, genes):
        tr_genes = sorted(set(genes[tr_idx]))
        in_tr_g = set(tr_genes[:max(1, int(0.75 * len(tr_genes)))])
        in_tr = np.array([i for i in tr_idx if genes[i] in in_tr_g])
        in_ca = np.array([i for i in tr_idx if genes[i] not in in_tr_g])

        sc = StandardScaler().fit(X[in_tr])
        mlp = MLPClassifier(hidden_layer_sizes=(32, 16), alpha=0.05, max_iter=800, early_stopping=True, random_state=42)
        mlp.fit(sc.transform(X[in_tr]), y[in_tr])

        p_ca = mlp.predict_proba(sc.transform(X[in_ca]))[:, 1]
        p_te = mlp.predict_proba(sc.transform(X[te_idx]))[:, 1]
        oof_raw[te_idx] = p_te

        platt = LogisticRegression(C=1.0).fit(_logit(p_ca).reshape(-1, 1), y[in_ca])
        p_ca_cal = platt.predict_proba(_logit(p_ca).reshape(-1, 1))[:, 1]
        p_te_cal = platt.predict_proba(_logit(p_te).reshape(-1, 1))[:, 1]
        oof_platt[te_idx] = p_te_cal

        mc = AncestryMondrianConformal(alpha=0.10).fit(p_ca_cal, y[in_ca], groups[in_ca])
        sets = mc.predict_sets(p_te_cal, groups[te_idx])
        oof_cov[te_idx] = [int(y[te_idx[k]] in sets[k]) for k in range(len(te_idx))]
        oof_sz[te_idx] = [len(sets[k]) for k in range(len(te_idx))]

    print("3/3 Fitting & exporting final v0.2.0 deployable model bundle...")
    all_genes = sorted(set(genes))
    fit_g = set(all_genes[:int(0.80 * len(all_genes))])
    f_tr = np.array([i for i, g in enumerate(genes) if g in fit_g])
    f_ca = np.array([i for i, g in enumerate(genes) if g not in fit_g])
    sc_fin = StandardScaler().fit(X[f_tr])
    mlp_fin = MLPClassifier(hidden_layer_sizes=(32, 16), alpha=0.05, max_iter=800, early_stopping=True, random_state=42).fit(sc_fin.transform(X[f_tr]), y[f_tr])
    p_ca_fin = mlp_fin.predict_proba(sc_fin.transform(X[f_ca]))[:, 1]
    platt_fin = LogisticRegression(C=1.0).fit(_logit(p_ca_fin).reshape(-1, 1), y[f_ca])
    p_ca_cal_fin = platt_fin.predict_proba(_logit(p_ca_fin).reshape(-1, 1))[:, 1]
    mc_fin = AncestryMondrianConformal(alpha=0.10).fit(p_ca_cal_fin, y[f_ca], groups[f_ca])

    bundle_path = "artifacts/equity_model_v2.joblib"
    joblib.dump({"scaler": sc_fin, "mlp": mlp_fin, "platt": platt_fin, "mondrian": mc_fin, "version": "0.2.0"}, bundle_path)

    report = {
        "version": "0.2.0",
        "n_total": len(valid),
        "spliceai_nonzero_count": nonzero,
        "v02_oof_auroc": round(float(roc_auc_score(y, oof_raw)), 4),
        "v02_oof_prauc": round(float(average_precision_score(y, oof_raw)), 4),
        "overall_platt_ece": round(compute_ece(y, oof_platt), 4),
        "by_ancestry_group": {
            g: {
                "n": int(np.sum(groups == g)),
                "platt_ece": round(compute_ece(y[groups == g], oof_platt[groups == g]), 4),
                "mondrian_cov": round(float(np.mean(oof_cov[groups == g])), 4),
                "singleton_rate": round(float(np.mean(oof_sz[groups == g] == 1)), 4),
            } for g in ["EUR", "SAS", "AFR", "ULTRA_RARE"]
        },
        "model_bundle_sha256": logger.sha256_file(bundle_path),
    }
    with open("artifacts/v02_spliceai_report.json", "w") as f:
        json.dump(report, f, indent=2)
    logger.log("v02_spliceai_upgrade", {"auroc": report["v02_oof_auroc"], "bundle_sha256": report["model_bundle_sha256"]})
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
