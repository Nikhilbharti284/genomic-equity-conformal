import csv, json, math, sqlite3
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from mondrian_conformal import AncestryMondrianConformal, compute_ece
from resource_logger import ResourceLogger

# Kyte-Doolittle hydrophobicity, molecular weight, formal charge (pH 7.4)
HYDRO = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
MW    = {"A":89,"R":174,"N":132,"D":133,"C":121,"Q":146,"E":147,"G":75,"H":155,"I":131,"L":131,"K":146,"M":149,"F":165,"P":115,"S":105,"T":119,"W":204,"Y":181,"V":117}
CHG   = {"R":1,"K":1,"H":0.1,"D":-1,"E":-1}

def _logit(p):
    p = np.clip(p, 1e-6, 1.0 - 1e-6)
    return np.log(p / (1.0 - p))

def main():
    logger = ResourceLogger()
    with open("data/pilot_variants_hg38.csv") as f:
        rows = list(csv.DictReader(f))
    conn = sqlite3.connect("artifacts/esm2_pilot.sqlite")
    scores = {r[0]: (r[1], r[2]) for r in conn.execute("SELECT variant_id, llr, emb_norm FROM scores")}
    conn.close()

    valid = [r for r in rows if r["variant_id"] in scores]
    X = []
    for r in valid:
        llr, enorm = scores[r["variant_id"]]
        wt, mt = r["wt"], r["mt"]
        d_hydro = HYDRO.get(mt, 0.0) - HYDRO.get(wt, 0.0)
        d_mw = (MW.get(mt, 120.0) - MW.get(wt, 120.0)) / 100.0
        d_chg = abs(CHG.get(mt, 0.0) - CHG.get(wt, 0.0))
        is_pro_gly = float(wt in ("P", "G") or mt in ("P", "G"))
        X.append([
            llr, enorm, d_hydro, abs(d_hydro), d_mw, abs(d_mw), d_chg, is_pro_gly,
            math.log10(float(r["af_nfe"]) + 1e-7),
            math.log10(float(r["af_sas"]) + 1e-7),
            math.log10(float(r["af_afr"]) + 1e-7),
            math.log10(float(r["af_all"]) + 1e-7),
        ])
    X = np.array(X, dtype=float)
    y = np.array([int(r["label"]) for r in valid])
    genes = np.array([r["gene"] for r in valid])
    groups = np.array([r["ancestry_group"] for r in valid])

    gkf = GroupKFold(n_splits=5)
    oof_raw = np.zeros(len(valid))
    oof_platt = np.zeros(len(valid))
    oof_global_cov = np.zeros(len(valid))
    oof_mondrian_cov = np.zeros(len(valid))
    oof_mondrian_size = np.zeros(len(valid))

    for tr_idx, te_idx in gkf.split(X, y, genes):
        tr_genes = sorted(set(genes[tr_idx]))
        cut = max(1, int(0.75 * len(tr_genes)))
        in_tr_g = set(tr_genes[:cut])
        in_tr = np.array([i for i in tr_idx if genes[i] in in_tr_g])
        in_ca = np.array([i for i in tr_idx if genes[i] not in in_tr_g])

        sc = StandardScaler().fit(X[in_tr])
        mlp = MLPClassifier(hidden_layer_sizes=(32, 16), alpha=0.05, max_iter=800, early_stopping=True, random_state=42)
        mlp.fit(sc.transform(X[in_tr]), y[in_tr])

        p_ca = mlp.predict_proba(sc.transform(X[in_ca]))[:, 1]
        p_te = mlp.predict_proba(sc.transform(X[te_idx]))[:, 1]
        oof_raw[te_idx] = p_te

        # Platt scaling (2-parameter logistic calibration on log-odds — low variance, no overfit)
        platt = LogisticRegression(C=1.0).fit(_logit(p_ca).reshape(-1, 1), y[in_ca])
        p_ca_cal = platt.predict_proba(_logit(p_ca).reshape(-1, 1))[:, 1]
        p_te_cal = platt.predict_proba(_logit(p_te).reshape(-1, 1))[:, 1]
        oof_platt[te_idx] = p_te_cal

        mc_g = AncestryMondrianConformal(alpha=0.10).fit(p_ca_cal, y[in_ca], np.array(["ALL"] * len(in_ca)))
        g_sets = mc_g.predict_sets(p_te_cal, np.array(["ALL"] * len(te_idx)))
        oof_global_cov[te_idx] = [int(y[te_idx[k]] in g_sets[k]) for k in range(len(te_idx))]

        mc_a = AncestryMondrianConformal(alpha=0.10).fit(p_ca_cal, y[in_ca], groups[in_ca])
        m_sets = mc_a.predict_sets(p_te_cal, groups[te_idx])
        oof_mondrian_cov[te_idx] = [int(y[te_idx[k]] in m_sets[k]) for k in range(len(te_idx))]
        oof_mondrian_size[te_idx] = [len(m_sets[k]) for k in range(len(te_idx))]

    by_group = {}
    for g in ["EUR", "SAS", "AFR", "ULTRA_RARE"]:
        idx = np.where(groups == g)[0]
        raw_ece = compute_ece(y[idx], oof_raw[idx])
        cal_ece = compute_ece(y[idx], oof_platt[idx])
        by_group[g] = {
            "n": int(len(idx)),
            "raw_ece": round(raw_ece, 4),
            "platt_ece": round(cal_ece, 4),
            "ece_reduction_pct": round((raw_ece - cal_ece) / max(raw_ece, 1e-9) * 100.0, 1),
            "global_cov": round(float(np.mean(oof_global_cov[idx])), 4),
            "mondrian_cov": round(float(np.mean(oof_mondrian_cov[idx])), 4),
            "avg_set_size": round(float(np.mean(oof_mondrian_size[idx])), 4),
            "singleton_rate": round(float(np.mean(oof_mondrian_size[idx] == 1)), 4),
        }

    report = {
        "n_total": len(valid),
        "baseline_raw_esm2_auroc": round(float(roc_auc_score(y, -X[:, 0])), 4),
        "v2_mlp_oof_auroc": round(float(roc_auc_score(y, oof_raw)), 4),
        "v2_mlp_oof_prauc": round(float(average_precision_score(y, oof_raw)), 4),
        "overall_raw_ece": round(compute_ece(y, oof_raw), 4),
        "overall_platt_ece": round(compute_ece(y, oof_platt), 4),
        "by_ancestry_group": by_group,
    }
    out_json = "artifacts/hg38_5fold_v2_report.json"
    with open(out_json, "w") as f:
        json.dump(report, f, indent=2)
    sha = logger.sha256_file(out_json)
    logger.log("hg38_5fold_v2", {"report_sha256": sha, "auroc": report["v2_mlp_oof_auroc"]})
    print(json.dumps(report, indent=2))
    print(f"Verified V2 Report SHA-256: {sha}")

if __name__ == "__main__":
    main()
