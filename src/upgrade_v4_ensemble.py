import csv, json, math, os, sqlite3, urllib.request
import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
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

def _num(x):
    if x is None: return None
    if isinstance(x, list):
        vals = [_num(v) for v in x]
        vals = [v for v in vals if v is not None]
        return max(vals) if vals else None
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except Exception:
        return None

def fetch_cadd_revel(hgvs_ids):
    out = {}
    url = "https://myvariant.info/v1/variant?assembly=hg38&fields=cadd.phred,dbnsfp.revel.score"
    for i in range(0, len(hgvs_ids), 250):
        payload = json.dumps({"ids": hgvs_ids[i:i+250]}).encode("utf-8")
        req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json", "User-Agent": "genomic-equity-conformal/0.4"})
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                for r in json.loads(resp.read().decode("utf-8")):
                    q = r.get("query")
                    if not q or r.get("notfound"): continue
                    cadd = r.get("cadd")
                    if isinstance(cadd, list) and cadd: cadd = cadd[0]
                    cadd_v = _num(cadd.get("phred")) if isinstance(cadd, dict) else None
                    db = r.get("dbnsfp")
                    if isinstance(db, list) and db: db = db[0]
                    rev_v = None
                    if isinstance(db, dict) and isinstance(db.get("revel"), dict):
                        rev_v = _num(db["revel"].get("score"))
                    out[q] = (cadd_v, rev_v)
        except Exception as e:
            print(f"[warn] CADD/REVEL batch {i}: {e}")
    return out

def bootstrap_ci(y, s1, s2=None, n_boot=1000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(y)
    vals1, diffs = [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y[idx])) < 2: continue
        a1 = roc_auc_score(y[idx], s1[idx])
        vals1.append(a1)
        if s2 is not None:
            diffs.append(a1 - roc_auc_score(y[idx], s2[idx]))
    ci = (round(float(np.percentile(vals1, 2.5)), 4), round(float(np.percentile(vals1, 97.5)), 4))
    pval = round(float(np.mean(np.array(diffs) <= 0.0)), 4) if diffs else None
    return ci, pval

def main():
    logger = ResourceLogger()
    with open("data/pilot_variants_hg38.csv") as f:
        rows = list(csv.DictReader(f))
    conn = sqlite3.connect("artifacts/esm2_pilot.sqlite")
    scores = {r[0]: (r[1], r[2]) for r in conn.execute("SELECT variant_id, llr, emb_norm FROM scores")}
    conn.close()

    print(f"1/3 Fetching hg38 CADD Phred & REVEL annotations for {len(rows)} variants...")
    ann = fetch_cadd_revel([r["hgvs_g"] for r in rows])
    cadd_vals = [ann.get(r["hgvs_g"], (None, None))[0] for r in rows]
    rev_vals  = [ann.get(r["hgvs_g"], (None, None))[1] for r in rows]
    n_cadd = sum(v is not None for v in cadd_vals)
    n_rev  = sum(v is not None for v in rev_vals)
    use_cadd = n_cadd >= 400
    use_rev  = n_rev >= 400
    med_cadd = float(np.median([v for v in cadd_vals if v is not None])) if use_cadd else 0.0
    med_rev  = float(np.median([v for v in rev_vals if v is not None])) if use_rev else 0.0
    print(f"    Coverage -> CADD: {n_cadd}/{len(rows)} (use={use_cadd}) | REVEL: {n_rev}/{len(rows)} (use={use_rev})")

    valid = [r for r in rows if r["variant_id"] in scores]
    X_list = []
    for r in valid:
        llr, enorm = scores[r["variant_id"]]
        wt, mt = r["wt"], r["mt"]
        d_hydro = HYDRO.get(mt, 0.0) - HYDRO.get(wt, 0.0)
        d_mw = (MW.get(mt, 120.0) - MW.get(wt, 120.0)) / 100.0
        d_chg = abs(CHG.get(mt, 0.0) - CHG.get(wt, 0.0))
        feat = [
            llr, enorm, d_hydro, abs(d_hydro), d_mw, abs(d_mw), d_chg,
            float(wt in ("P","G") or mt in ("P","G")),
            math.log10(float(r["af_nfe"]) + 1e-7), math.log10(float(r["af_sas"]) + 1e-7),
            math.log10(float(r["af_afr"]) + 1e-7), math.log10(float(r["af_all"]) + 1e-7),
        ]
        c_v, r_v = ann.get(r["hgvs_g"], (None, None))
        if use_cadd: feat.append(c_v if c_v is not None else med_cadd)
        if use_rev:  feat.append(r_v if r_v is not None else med_rev)
        X_list.append(feat)

    X = np.array(X_list, dtype=float)
    y = np.array([int(r["label"]) for r in valid])
    genes = np.array([r["gene"] for r in valid])
    groups = np.array([r["ancestry_group"] for r in valid])

    print("2/3 Running 5-Fold Zero-Leakage CV with 3-Model Ensemble (MLP + HistGBM + LogReg)...")
    gkf = GroupKFold(n_splits=5)
    oof_ens, oof_platt = np.zeros(len(valid)), np.zeros(len(valid))
    oof_cov, oof_sz = np.zeros(len(valid)), np.zeros(len(valid))

    def fit_ensemble(Xtr, ytr):
        sc = StandardScaler().fit(Xtr)
        Xtr_s = sc.transform(Xtr)
        m1 = MLPClassifier(hidden_layer_sizes=(32, 16), alpha=0.05, max_iter=800, early_stopping=True, random_state=42).fit(Xtr_s, ytr)
        m2 = HistGradientBoostingClassifier(max_depth=3, learning_rate=0.05, max_iter=150, l2_regularization=1.0, random_state=42).fit(Xtr, ytr)
        m3 = LogisticRegression(C=0.5, max_iter=500).fit(Xtr_s, ytr)
        return sc, m1, m2, m3

    def pred_ensemble(models, Xq):
        sc, m1, m2, m3 = models
        Xq_s = sc.transform(Xq)
        return 0.45 * m1.predict_proba(Xq_s)[:, 1] + 0.35 * m2.predict_proba(Xq)[:, 1] + 0.20 * m3.predict_proba(Xq_s)[:, 1]

    for tr_idx, te_idx in gkf.split(X, y, genes):
        tr_genes = sorted(set(genes[tr_idx]))
        in_tr_g = set(tr_genes[:max(1, int(0.75 * len(tr_genes)))])
        in_tr = np.array([i for i in tr_idx if genes[i] in in_tr_g])
        in_ca = np.array([i for i in tr_idx if genes[i] not in in_tr_g])

        models = fit_ensemble(X[in_tr], y[in_tr])
        p_ca = pred_ensemble(models, X[in_ca])
        p_te = pred_ensemble(models, X[te_idx])
        oof_ens[te_idx] = p_te

        platt = LogisticRegression(C=1.0).fit(_logit(p_ca).reshape(-1, 1), y[in_ca])
        p_ca_cal = platt.predict_proba(_logit(p_ca).reshape(-1, 1))[:, 1]
        p_te_cal = platt.predict_proba(_logit(p_te).reshape(-1, 1))[:, 1]
        oof_platt[te_idx] = p_te_cal

        mc = AncestryMondrianConformal(alpha=0.10).fit(p_ca_cal, y[in_ca], groups[in_ca])
        sets = mc.predict_sets(p_te_cal, groups[te_idx])
        oof_cov[te_idx] = [int(y[te_idx[k]] in sets[k]) for k in range(len(te_idx))]
        oof_sz[te_idx]  = [len(sets[k]) for k in range(len(te_idx))]

    print("3/3 Computing 1,000-Bootstrap 95% CIs, Paired p-value & Exporting v0.4.0 Ensemble...")
    raw_esm2 = -X[:, 0]
    esm2_ci, _ = bootstrap_ci(y, raw_esm2)
    ens_ci, p_vs_esm2 = bootstrap_ci(y, oof_ens, raw_esm2)

    all_genes = sorted(set(genes))
    fit_g = set(all_genes[:int(0.80 * len(all_genes))])
    f_tr = np.array([i for i, g in enumerate(genes) if g in fit_g])
    f_ca = np.array([i for i, g in enumerate(genes) if g not in fit_g])
    fin_models = fit_ensemble(X[f_tr], y[f_tr])
    p_ca_fin = pred_ensemble(fin_models, X[f_ca])
    platt_fin = LogisticRegression(C=1.0).fit(_logit(p_ca_fin).reshape(-1, 1), y[f_ca])
    mc_fin = AncestryMondrianConformal(alpha=0.10).fit(platt_fin.predict_proba(_logit(p_ca_fin).reshape(-1, 1))[:, 1], y[f_ca], groups[f_ca])

    bundle_path = "artifacts/equity_model_v4.joblib"
    joblib.dump({
        "scaler": fin_models[0], "mlp": fin_models[1], "hgb": fin_models[2], "logreg": fin_models[3],
        "platt": platt_fin, "mondrian": mc_fin,
        "use_cadd": use_cadd, "use_rev": use_rev, "med_cadd": med_cadd, "med_rev": med_rev,
        "version": "0.4.0"
    }, bundle_path)

    report = {
        "version": "0.4.0",
        "n_total": len(valid),
        "features_used": {"n_features": int(X.shape[1]), "cadd_included": use_cadd, "revel_included": use_rev},
        "baseline_raw_esm2_auroc": round(float(roc_auc_score(y, raw_esm2)), 4),
        "baseline_raw_esm2_95ci": list(esm2_ci),
        "v4_ensemble_oof_auroc": round(float(roc_auc_score(y, oof_ens)), 4),
        "v4_ensemble_95ci": list(ens_ci),
        "paired_bootstrap_p_vs_esm2": p_vs_esm2,
        "v4_ensemble_oof_prauc": round(float(average_precision_score(y, oof_ens)), 4),
        "overall_platt_ece": round(compute_ece(y, oof_platt), 4),
        "by_ancestry_group": {
            g: {
                "n": int(np.sum(groups == g)),
                "platt_ece": round(compute_ece(y[groups == g], oof_platt[groups == g]), 4),
                "mondrian_cov": round(float(np.mean(oof_cov[groups == g])), 4),
                "singleton_rate": round(float(np.mean(oof_sz[groups == g] == 1)), 4),
            } for g in ["EUR", "SAS", "AFR", "ULTRA_RARE"]
        },
        "bundle_sha256": logger.sha256_file(bundle_path),
    }
    out_json = "artifacts/v04_scientific_proof.json"
    with open(out_json, "w") as f:
        json.dump(report, f, indent=2)
    logger.log("v04_ensemble_proof", {"auroc": report["v4_ensemble_oof_auroc"], "p_val": p_vs_esm2})
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
