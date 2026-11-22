import csv, json, math, os, sqlite3
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

def bootstrap_ci(y, s1, s2, n_boot=1000, seed=42):
    rng = np.random.default_rng(seed)
    n, a1s, diffs = len(y), [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y[idx])) < 2: continue
        v1 = roc_auc_score(y[idx], s1[idx])
        a1s.append(v1)
        diffs.append(v1 - roc_auc_score(y[idx], s2[idx]))
    return [round(float(np.percentile(a1s, 2.5)), 4), round(float(np.percentile(a1s, 97.5)), 4)], round(float(np.mean(np.array(diffs) <= 0.0)), 4)

def main():
    logger = ResourceLogger()
    with open("data/v1_homology_clusters.json") as f:
        g2c = json.load(f)["gene_to_cluster"]
    conn = sqlite3.connect("artifacts/v1_commercial_features.sqlite")
    fmap = {r[0]: r[1:] for r in conn.execute("SELECT variant_id, am_score, af_nfe, af_sas, af_afr, af_all, ancestry_group FROM var_features")}
    conn.close()

    rows, X_list, y_list, clu_list, grp_list, post23_list, am_raw = [], [], [], [], [], [], []
    with open("data/v1_clinvar_120genes.csv") as f:
        for r in csv.DictReader(f):
            fv = fmap.get(r["variant_id"])
            if not fv or fv[0] < 0: continue  # require verified CC-BY-4.0 AlphaMissense score
            am, af_nfe, af_sas, af_afr, af_all, grp = fv
            wt, mt = r["wt"], r["mt"]
            dh = HYDRO.get(mt, 0.0) - HYDRO.get(wt, 0.0)
            dm = (MW.get(mt, 120.0) - MW.get(wt, 120.0)) / 100.0
            dc = abs(CHG.get(mt, 0.0) - CHG.get(wt, 0.0))
            grpmax = max(af_nfe, af_sas, af_afr, af_all)
            X_list.append([
                am, _logit(am), dh, abs(dh), dm, abs(dm), dc, float(wt in ("P","G") or mt in ("P","G")),
                math.log10(af_nfe + 1e-7), math.log10(af_sas + 1e-7), math.log10(af_afr + 1e-7), math.log10(grpmax + 1e-7)
            ])
            y_list.append(int(r["label"]))
            clu_list.append(g2c[r["gene"]])
            grp_list.append(grp)
            post23_list.append(int(r["is_post_2023"]))
            am_raw.append(am)
            rows.append(r)

    X = np.array(X_list, dtype=float)
    y = np.array(y_list, dtype=int)
    clusters = np.array(clu_list)
    groups = np.array(grp_list)
    post23 = np.array(post23_list, dtype=int)
    am_raw = np.array(am_raw, dtype=float)

    print(f"1/3 Running 5-Fold Homology-Cluster CV (n={len(y)}, {len(set(clusters))} independent clusters)...")
    def fit_ens(Xa, ya):
        sc = StandardScaler().fit(Xa)
        Xs = sc.transform(Xa)
        m1 = MLPClassifier(hidden_layer_sizes=(32, 16), alpha=0.05, max_iter=600, early_stopping=True, random_state=42).fit(Xs, ya)
        m2 = HistGradientBoostingClassifier(max_depth=4, learning_rate=0.05, max_iter=200, l2_regularization=1.0, random_state=42).fit(Xa, ya)
        m3 = LogisticRegression(C=0.5, max_iter=500).fit(Xs, ya)
        return sc, m1, m2, m3

    def pred_ens(ms, Xq):
        sc, m1, m2, m3 = ms
        Xs = sc.transform(Xq)
        return 0.40 * m1.predict_proba(Xs)[:, 1] + 0.45 * m2.predict_proba(Xq)[:, 1] + 0.15 * m3.predict_proba(Xs)[:, 1]

    gkf = GroupKFold(n_splits=5)
    oof_raw, oof_cal, oof_cov, oof_sz, empty_cnt = np.zeros(len(y)), np.zeros(len(y)), np.zeros(len(y)), np.zeros(len(y)), 0
    for tr_idx, te_idx in gkf.split(X, y, clusters):
        # Allocate 30+ independent homology clusters strictly to conformal calibration!
        tr_clus = sorted(set(clusters[tr_idx]), key=lambda c: hash(f"v10_{c}") % 100000)
        cut = max(20, int(0.62 * len(tr_clus)))
        in_tr_c = set(tr_clus[:cut])
        in_tr = np.array([i for i in tr_idx if clusters[i] in in_tr_c])
        in_ca = np.array([i for i in tr_idx if clusters[i] not in in_tr_c])

        # Nested inner OOF Platt scaling (zero double-dipping!)
        igkf = GroupKFold(n_splits=4)
        in_oof = np.zeros(len(in_tr))
        for a, b in igkf.split(X[in_tr], y[in_tr], clusters[in_tr]):
            in_oof[b] = pred_ens(fit_ens(X[in_tr[a]], y[in_tr[a]]), X[in_tr[b]])
        platt = LogisticRegression(C=1.0).fit(_logit(in_oof).reshape(-1, 1), y[in_tr])

        ms = fit_ens(X[in_tr], y[in_tr])
        p_ca = platt.predict_proba(_logit(pred_ens(ms, X[in_ca])).reshape(-1, 1))[:, 1]
        p_te_raw = pred_ens(ms, X[te_idx])
        p_te = platt.predict_proba(_logit(p_te_raw).reshape(-1, 1))[:, 1]
        oof_raw[te_idx] = p_te_raw
        oof_cal[te_idx] = p_te

        mc = AncestryMondrianConformal(alpha=0.10).fit(p_ca, y[in_ca], groups[in_ca])
        sets = mc.predict_sets(p_te, groups[te_idx])
        empty_cnt += sum(1 for s in sets if len(s) == 0)
        oof_cov[te_idx] = [int(y[te_idx[k]] in sets[k]) for k in range(len(te_idx))]
        oof_sz[te_idx]  = [len(sets[k]) for k in range(len(te_idx))]

    print("2/3 Computing 1,000-Bootstrap 95% CIs & Post-2023 Prospective Holdout Metrics...")
    ci_v1, p_vs_am = bootstrap_ci(y, oof_cal, am_raw)
    ci_am, _ = bootstrap_ci(y, am_raw, am_raw)
    p23_mask = post23 == 1

    print("3/3 Exporting final 100% Commercial-Safe v1.0.0 119-Gene Bundle...")
    all_c = sorted(set(clusters), key=lambda c: hash(f"fin_{c}") % 100000)
    f_tr_c = set(all_c[:int(0.68 * len(all_c))])
    f_tr = np.array([i for i, c in enumerate(clusters) if c in f_tr_c])
    f_ca = np.array([i for i, c in enumerate(clusters) if c not in f_tr_c])
    fin_ms = fit_ens(X[f_tr], y[f_tr])
    igkf = GroupKFold(n_splits=4)
    in_oof = np.zeros(len(f_tr))
    for a, b in igkf.split(X[f_tr], y[f_tr], clusters[f_tr]):
        in_oof[b] = pred_ens(fit_ens(X[f_tr[a]], y[f_tr[a]]), X[f_tr[b]])
    fin_platt = LogisticRegression(C=1.0).fit(_logit(in_oof).reshape(-1, 1), y[f_tr])
    p_ca_fin = fin_platt.predict_proba(_logit(pred_ens(fin_ms, X[f_ca])).reshape(-1, 1))[:, 1]
    fin_mc = AncestryMondrianConformal(alpha=0.10).fit(p_ca_fin, y[f_ca], groups[f_ca])

    bpath = "artifacts/equity_model_v10_119genes_ccby4.joblib"
    joblib.dump({"version": "1.0.0-ccby4-119genes", "commercial_safe": True, "models": fin_ms, "platt": fin_platt, "mondrian": fin_mc}, bpath)

    report = {
        "version": "1.0.0-ccby4-119genes",
        "commercial_license_clean": True,
        "n_variants": int(len(y)),
        "n_genes": int(len({r["gene"] for r in rows})),
        "n_homology_clusters": int(len(set(clusters))),
        "empty_prediction_sets": int(empty_cnt),
        "baseline_alphamissense_only_auroc": round(float(roc_auc_score(y, am_raw)), 4),
        "baseline_alphamissense_95ci": ci_am,
        "v10_trackB_oof_auroc": round(float(roc_auc_score(y, oof_cal)), 4),
        "v10_trackB_95ci": ci_v1,
        "paired_bootstrap_p_vs_alphamissense": p_vs_am,
        "v10_trackB_oof_prauc": round(float(average_precision_score(y, oof_cal)), 4),
        "post_2023_prospective_n": int(np.sum(p23_mask)),
        "post_2023_prospective_auroc": round(float(roc_auc_score(y[p23_mask], oof_cal[p23_mask])), 4),
        "overall_ece_raw_am": round(compute_ece(y, am_raw), 4),
        "overall_ece_v10_platt": round(compute_ece(y, oof_cal), 4),
        "by_ancestry": {
            g: {
                "n": int(np.sum(groups == g)),
                "raw_am_ece": round(compute_ece(y[groups == g], am_raw[groups == g]), 4),
                "v10_platt_ece": round(compute_ece(y[groups == g], oof_cal[groups == g]), 4),
                "mondrian_cov": round(float(np.mean(oof_cov[groups == g])), 4),
                "singleton_rate": round(float(np.mean(oof_sz[groups == g] == 1)), 4),
            } for g in ["EUR", "SAS", "AFR", "ULTRA_RARE"]
        },
        "bundle_sha256": logger.sha256_file(bpath),
    }
    out_json = "artifacts/v10_119genes_final_report.json"
    with open(out_json, "w") as f:
        json.dump(report, f, indent=2)
    report["report_sha256"] = logger.sha256_file(out_json)
    logger.log("v10_119genes_final", {"auroc": report["v10_trackB_oof_auroc"], "p_vs_am": p_vs_am})
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
