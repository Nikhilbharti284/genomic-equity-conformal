import csv, json, math, os, sqlite3, urllib.request
import numpy as np
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from mondrian_conformal import AncestryMondrianConformal, compute_ece
from resource_logger import ResourceLogger

def _num(x):
    if x is None:
        return 0.0
    if isinstance(x, list):
        vals = [_num(v) for v in x]
        return max(vals) if vals else 0.0
    try:
        v = float(x)
        return 0.0 if math.isnan(v) else v
    except Exception:
        return 0.0

def extract_af(rec):
    # Try gnomad_exome, gnomad_genome, then dbnsfp fallback
    for key in ("gnomad_exome", "gnomad_genome"):
        sub = rec.get(key)
        if isinstance(sub, list) and sub:
            sub = sub[0]
        if isinstance(sub, dict) and isinstance(sub.get("af"), dict):
            af = sub["af"]
            vals = {
                "af_nfe": _num(af.get("af_nfe")),
                "af_sas": _num(af.get("af_sas")),
                "af_afr": _num(af.get("af_afr")),
                "af_all": _num(af.get("af")),
            }
            if max(vals.values()) > 0:
                return vals
    db = rec.get("dbnsfp")
    if isinstance(db, list) and db:
        db = db[0]
    if isinstance(db, dict):
        for key in ("gnomad_exomes", "gnomad_genomes"):
            sub = db.get(key)
            if isinstance(sub, dict):
                af = sub.get("af") if isinstance(sub.get("af"), dict) else sub
                vals = {
                    "af_nfe": _num(af.get("af_nfe") or sub.get("nfe_af")),
                    "af_sas": _num(af.get("af_sas") or sub.get("sas_af")),
                    "af_afr": _num(af.get("af_afr") or sub.get("afr_af")),
                    "af_all": _num(af.get("af") or sub.get("af")),
                }
                if max(vals.values()) > 0:
                    return vals
    return {"af_nfe": 0.0, "af_sas": 0.0, "af_afr": 0.0, "af_all": 0.0}

def fetch_hg38_gnomad(hgvs_ids):
    out = {}
    url = "https://myvariant.info/v1/variant?assembly=hg38&fields=gnomad_exome.af,gnomad_genome.af,dbnsfp.gnomad_exomes,dbnsfp.gnomad_genomes"
    for i in range(0, len(hgvs_ids), 250):
        chunk = hgvs_ids[i:i+250]
        payload = json.dumps({"ids": chunk}).encode("utf-8")
        req = urllib.request.Request(
            url, data=payload,
            headers={"Content-Type": "application/json", "User-Agent": "genomic-equity-conformal/0.2"}
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                rows = json.loads(resp.read().decode("utf-8"))
                for r in rows:
                    q = r.get("query")
                    if q and not r.get("notfound"):
                        out[q] = extract_af(r)
        except Exception as e:
            print(f"[warn] batch {i} error: {e}")
    return out

def main():
    logger = ResourceLogger()
    with open("data/pilot_variants.csv") as f:
        rows = list(csv.DictReader(f))

    print(f"1/3 Fetching hg38 gnomAD ancestry AFs for {len(rows)} variants...")
    afs = fetch_hg38_gnomad([r["hgvs_g"] for r in rows])

    grp_counts = {"EUR": 0, "SAS": 0, "AFR": 0, "ULTRA_RARE": 0}
    for r in rows:
        af = afs.get(r["hgvs_g"], {"af_nfe": 0.0, "af_sas": 0.0, "af_afr": 0.0, "af_all": 0.0})
        r.update(af)
        if af["af_sas"] > 0 and af["af_sas"] >= max(af["af_nfe"], af["af_afr"]):
            grp = "SAS"
        elif af["af_afr"] > 0 and af["af_afr"] >= max(af["af_nfe"], af["af_sas"]):
            grp = "AFR"
        elif af["af_nfe"] > 0:
            grp = "EUR"
        else:
            grp = "ULTRA_RARE"
        r["ancestry_group"] = grp
        grp_counts[grp] += 1

    print("    Ancestry counts (hg38):", grp_counts)
    out_csv = "data/pilot_variants_hg38.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print("2/3 Loading cached ESM-2 scores from artifacts/esm2_pilot.sqlite (0s compute)...")
    conn = sqlite3.connect("artifacts/esm2_pilot.sqlite")
    scores = {r[0]: (r[1], r[2]) for r in conn.execute("SELECT variant_id, llr, emb_norm FROM scores")}
    conn.close()

    valid = [r for r in rows if r["variant_id"] in scores]
    X = np.array([[
        scores[r["variant_id"]][0],
        scores[r["variant_id"]][1],
        math.log10(float(r["af_nfe"]) + 1e-7),
        math.log10(float(r["af_sas"]) + 1e-7),
        math.log10(float(r["af_afr"]) + 1e-7),
        math.log10(float(r["af_all"]) + 1e-7),
    ] for r in valid])
    y = np.array([int(r["label"]) for r in valid])
    genes = np.array([r["gene"] for r in valid])
    groups = np.array([r["ancestry_group"] for r in valid])

    print("3/3 Running 5-Fold Zero-Leakage Gene-Holdout CV + Ancestry Mondrian Conformal...")
    gkf = GroupKFold(n_splits=5)
    oof_raw_p = np.zeros(len(valid))
    oof_cal_p = np.zeros(len(valid))
    oof_global_cov = np.zeros(len(valid))
    oof_mondrian_cov = np.zeros(len(valid))
    oof_mondrian_size = np.zeros(len(valid))

    for tr_idx, te_idx in gkf.split(X, y, genes):
        # Split train genes into inner-train (75%) and calibration (25%) by gene
        tr_genes = sorted(set(genes[tr_idx]))
        cut = max(1, int(0.75 * len(tr_genes)))
        in_tr_g = set(tr_genes[:cut])
        in_tr = np.array([i for i in tr_idx if genes[i] in in_tr_g])
        in_ca = np.array([i for i in tr_idx if genes[i] not in in_tr_g])

        sc = StandardScaler().fit(X[in_tr])
        mlp = MLPClassifier(hidden_layer_sizes=(32, 16), max_iter=300, random_state=42)
        mlp.fit(sc.transform(X[in_tr]), y[in_tr])

        p_ca = mlp.predict_proba(sc.transform(X[in_ca]))[:, 1]
        p_te = mlp.predict_proba(sc.transform(X[te_idx]))[:, 1]
        oof_raw_p[te_idx] = p_te

        # Isotonic probability calibration on calibration genes
        iso = IsotonicRegression(out_of_bounds="clip").fit(p_ca, y[in_ca])
        p_ca_cal = iso.predict(p_ca)
        p_te_cal = iso.predict(p_te)
        oof_cal_p[te_idx] = p_te_cal

        # Global conformal vs Ancestry-Conditional Mondrian Conformal
        mc_global = AncestryMondrianConformal(alpha=0.10).fit(p_ca_cal, y[in_ca], np.array(["ALL"] * len(in_ca)))
        g_sets = mc_global.predict_sets(p_te_cal, np.array(["ALL"] * len(te_idx)))
        oof_global_cov[te_idx] = [int(y[te_idx[k]] in g_sets[k]) for k in range(len(te_idx))]

        mc_anc = AncestryMondrianConformal(alpha=0.10).fit(p_ca_cal, y[in_ca], groups[in_ca])
        m_sets = mc_anc.predict_sets(p_te_cal, groups[te_idx])
        oof_mondrian_cov[te_idx] = [int(y[te_idx[k]] in m_sets[k]) for k in range(len(te_idx))]
        oof_mondrian_size[te_idx] = [len(m_sets[k]) for k in range(len(te_idx))]

    by_group = {}
    for g in ["EUR", "SAS", "AFR", "ULTRA_RARE"]:
        idx = np.where(groups == g)[0]
        if len(idx) == 0:
            continue
        by_group[g] = {
            "n": int(len(idx)),
            "raw_ece": round(compute_ece(y[idx], oof_raw_p[idx]), 4),
            "calibrated_ece": round(compute_ece(y[idx], oof_cal_p[idx]), 4),
            "global_conformal_coverage": round(float(np.mean(oof_global_cov[idx])), 4),
            "mondrian_conformal_coverage": round(float(np.mean(oof_mondrian_cov[idx])), 4),
            "mondrian_avg_set_size": round(float(np.mean(oof_mondrian_size[idx])), 4),
            "mondrian_singleton_rate": round(float(np.mean(oof_mondrian_size[idx] == 1)), 4),
        }

    report = {
        "n_total": len(valid),
        "cv_scheme": "5-Fold GroupKFold (Zero Gene Leakage)",
        "baseline_raw_esm2_auroc": round(float(roc_auc_score(y, -X[:, 0])), 4),
        "mlp_oof_auroc": round(float(roc_auc_score(y, oof_raw_p)), 4),
        "mlp_oof_prauc": round(float(average_precision_score(y, oof_raw_p)), 4),
        "overall_raw_ece": round(compute_ece(y, oof_raw_p), 4),
        "overall_calibrated_ece": round(compute_ece(y, oof_cal_p), 4),
        "by_ancestry_group": by_group,
    }
    out_json = "artifacts/hg38_5fold_report.json"
    with open(out_json, "w") as f:
        json.dump(report, f, indent=2)
    sha = logger.sha256_file(out_json)
    logger.log("hg38_5fold_benchmark", {"report_sha256": sha, "mlp_oof_auroc": report["mlp_oof_auroc"]})
    print(json.dumps(report, indent=2))
    print(f"Verified Report SHA-256: {sha}")

if __name__ == "__main__":
    main()
