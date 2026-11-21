import csv, json, math, sqlite3
import numpy as np
from resource_logger import ResourceLogger

def empirical_lr_bootstrap(y, scores, threshold, direction="ge", n_boot=2000, seed=42):
    """Computes empirical LR+ and 1-sided 95% conservative bootstrap bound (Pejaver 2022)."""
    y = np.asarray(y, dtype=int)
    scores = np.asarray(scores, dtype=float)
    rng = np.random.default_rng(seed)
    n = len(y)
    lrs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb, sb = y[idx], scores[idx]
        pos, neg = (yb == 1), (yb == 0)
        if np.sum(pos) == 0 or np.sum(neg) == 0: continue
        hit = (sb >= threshold) if direction == "ge" else (sb <= threshold)
        # Laplace smoothing (+0.5) prevents divide-by-zero on extreme tails
        tpr = (np.sum(hit & pos) + 0.5) / (np.sum(pos) + 1.0)
        fpr = (np.sum(hit & neg) + 0.5) / (np.sum(neg) + 1.0)
        lrs.append(tpr / fpr)
    lrs = np.array(lrs)
    if direction == "ge":
        # Pathogenic PP3: conservative = 5th percentile lower bound
        return round(float(np.median(lrs)), 3), round(float(np.percentile(lrs, 5)), 3)
    else:
        # Benign BP4: conservative = 95th percentile upper bound (must stay <= target)
        return round(float(np.median(lrs)), 3), round(float(np.percentile(lrs, 95)), 3)

def main():
    logger = ResourceLogger()
    conn = sqlite3.connect("artifacts/v1_commercial_features.sqlite")
    fmap = {r[0]: (r[1], r[2]) for r in conn.execute("SELECT variant_id, am_score, ancestry_group FROM var_features")}
    conn.close()

    y_all, am_all, grp_all = [], [], []
    with open("data/v1_clinvar_120genes.csv") as f:
        for r in csv.DictReader(f):
            fv = fmap.get(r["variant_id"])
            if not fv or fv[0] < 0: continue
            y_all.append(int(r["label"]))
            am_all.append(float(fv[0]))
            grp_all.append(fv[1])

    y_all, am_all, grp_all = np.array(y_all), np.array(am_all), np.array(grp_all)

    # Pejaver (2022) / Bergquist (2025) AF-free AlphaMissense empirical LR+ verification per ancestry!
    tiers = [
        ("PP3_Strong (+4, LR+>=33.53)",     0.97, "ge", 33.526),
        ("PP3_Moderate (+2, LR+>=5.79)",    0.85, "ge", 5.790),
        ("PP3_Supporting (+1, LR+>=2.41)",  0.56, "ge", 2.406),
        ("BP4_Supporting (-1, LR+<=0.416)", 0.34, "le", 0.4156),
        ("BP4_Moderate (-2, LR+<=0.173)",   0.15, "le", 0.1727),
        ("BP4_Strong (-4, LR+<=0.030)",     0.07, "le", 0.0298),
    ]

    report = {"n_total": int(len(y_all)), "signal": "AF-Free CC-BY-4.0 AlphaMissense", "bootstrap_replicates": 2000, "subgroups": {}}
    for g in ["ALL", "EUR", "SAS", "AFR"]:
        mask = np.ones(len(y_all), dtype=bool) if g == "ALL" else (grp_all == g)
        sub_res = {}
        for name, thr, d, target in tiers:
            med_lr, bound_lr = empirical_lr_bootstrap(y_all[mask], am_all[mask], thr, direction=d, n_boot=2000)
            passed = (bound_lr >= target) if d == "ge" else (bound_lr <= target)
            sub_res[name] = {"threshold": thr, "median_lr": med_lr, "conservative_95_bound": bound_lr, "pejaver_target": target, "meets_clingen_95_bound": bool(passed)}
        report["subgroups"][g] = {"n": int(np.sum(mask)), "tiers": sub_res}

    out_json = "artifacts/v10_pejaver_ancestry_lr_audit.json"
    with open(out_json, "w") as f:
        json.dump(report, f, indent=2)
    sha = logger.sha256_file(out_json)
    print(json.dumps(report, indent=2))
    print(f"Pejaver Empirical LR+ Audit SHA-256: {sha}")

if __name__ == "__main__":
    main()
