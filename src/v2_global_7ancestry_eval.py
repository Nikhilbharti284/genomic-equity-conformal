import csv, hashlib, json, math, os, sqlite3, sys, urllib.request
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mondrian_conformal import AncestryMondrianConformal, compute_ece

HYDRO = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
MW    = {"A":89,"R":174,"N":132,"D":133,"C":121,"Q":146,"E":147,"G":75,"H":155,"I":131,"L":131,"K":146,"M":149,"F":165,"P":115,"S":105,"T":119,"W":204,"Y":181,"V":117}
CHG   = {"R":1,"K":1,"H":0.1,"D":-1,"E":-1}

def _num(x):
    if x is None: return 0.0
    if isinstance(x, list):
        vals = [_num(v) for v in x]
        return max(vals) if vals else 0.0
    try:
        v = float(x)
        return 0.0 if math.isnan(v) else v
    except Exception:
        return 0.0

def _logit(p):
    p = np.clip(p, 1e-6, 1.0 - 1e-6)
    return np.log(p / (1.0 - p))

def assign_7ancestry_stratum(afs):
    # afs dict: sas, afr, eas, amr, mid, nfe, asj, fin, all
    founder = max(afs["asj"], afs["fin"])
    cand = {
        "SAS": afs["sas"], "AFR": afs["afr"], "EAS": afs["eas"],
        "AMR": afs["amr"], "MID": afs["mid"], "NFE": afs["nfe"],
        "BOTTLENECK_FOUNDER": founder
    }
    best_g, best_v = max(cand.items(), key=lambda kv: kv[1])
    return best_g if best_v > 0.0 else "ULTRA_RARE"

def fetch_7ancestry_afs(hgvs_list, cache_db="artifacts/v2_global_7ancestry_features.sqlite"):
    os.makedirs(os.path.dirname(cache_db), exist_ok=True)
    conn = sqlite3.connect(cache_db)
    conn.execute("""CREATE TABLE IF NOT EXISTS gnomad7 (
        hgvs_g TEXT PRIMARY KEY,
        af_sas REAL, af_afr REAL, af_eas REAL, af_amr REAL, af_mid REAL,
        af_nfe REAL, af_asj REAL, af_fin REAL, af_all REAL
    )""")
    existing = {r[0]: r[1:] for r in conn.execute("SELECT * FROM gnomad7")}
    todo = [h for h in hgvs_list if h and h not in existing]
    if todo:
        print(f"Fetching 7-ancestry gnomAD AFs for {len(todo)} variants via MyVariant.info...")
        url = "https://myvariant.info/v1/variant?assembly=hg38&fields=gnomad_exome.af,gnomad_genome.af"
        for i in range(0, len(todo), 250):
            batch = todo[i:i+250]
            payload = json.dumps({"ids": batch}).encode("utf-8")
            req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json", "User-Agent": "genomic-equity-conformal/2.0"})
            rows_ins = []
            try:
                with urllib.request.urlopen(req, timeout=30) as resp:
                    for r in json.loads(resp.read().decode("utf-8")):
                        q = r.get("query")
                        if not q: continue
                        vals = [0.0]*9
                        for k in ("gnomad_exome", "gnomad_genome"):
                            sub = r.get(k)
                            if isinstance(sub, list) and sub: sub = sub[0]
                            if isinstance(sub, dict) and isinstance(sub.get("af"), dict):
                                a = sub["af"]
                                cur = [
                                    _num(a.get("af_sas")), _num(a.get("af_afr")), _num(a.get("af_eas")),
                                    _num(a.get("af_amr")), _num(a.get("af_mid")), _num(a.get("af_nfe")),
                                    _num(a.get("af_asj")), _num(a.get("af_fin")), _num(a.get("af"))
                                ]
                                vals = [max(v1, v2) for v1, v2 in zip(vals, cur)]
                        rows_ins.append((q, *vals))
                        existing[q] = tuple(vals)
            except Exception as e:
                print(f"  [WARN] batch {i}: {e}")
            if rows_ins:
                conn.executemany("INSERT OR REPLACE INTO gnomad7 VALUES (?,?,?,?,?,?,?,?,?,?)", rows_ins)
                conn.commit()
    conn.close()
    return existing

def main():
    with open("data/v1_homology_clusters.json") as f:
        g2c = json.load(f)["gene_to_cluster"]
    conn = sqlite3.connect("artifacts/v1_commercial_features.sqlite")
    fmap = {r[0]: r[1:] for r in conn.execute("SELECT variant_id, am_score, af_nfe, af_sas, af_afr, af_all, ancestry_group FROM var_features")}
    conn.close()

    clin_rows = []
    with open("data/v1_clinvar_120genes.csv") as f:
        for r in csv.DictReader(f):
            fv = fmap.get(r["variant_id"])
            if not fv or fv[0] < 0: continue
            hgvs = r.get("hgvs_g") or f"chr{str(r.get('chrom','')).replace('chr','')}:g.{r.get('pos_g') or r.get('pos','')}{r.get('ref','')}>{r.get('alt','')}"
            r["_hgvs"] = hgvs
            r["_am"] = float(fv[0])
            r["_fv"] = fv
            clin_rows.append(r)

    g7 = fetch_7ancestry_afs(list({r["_hgvs"] for r in clin_rows}))

    X_list, y_list, clu_list, grp_list, am_raw = [], [], [], [], []
    for r in clin_rows:
        am, old_nfe, old_sas, old_afr, old_all, _ = r["_fv"]
        g = g7.get(r["_hgvs"], (old_sas, old_afr, 0.0, 0.0, 0.0, old_nfe, 0.0, 0.0, old_all))
        afs = {
            "sas": max(g[0], old_sas), "afr": max(g[1], old_afr), "eas": g[2],
            "amr": g[3], "mid": g[4], "nfe": max(g[5], old_nfe),
            "asj": g[6], "fin": g[7], "all": max(g[8], old_all)
        }
        grp = assign_7ancestry_stratum(afs)
        wt, mt = r["wt"], r["mt"]
        dh = HYDRO.get(mt, 0.0) - HYDRO.get(wt, 0.0)
        dm = (MW.get(mt, 120.0) - MW.get(wt, 120.0)) / 100.0
        dc = abs(CHG.get(mt, 0.0) - CHG.get(wt, 0.0))
        grpmax = max(afs.values())
        X_list.append([
            am, dh, abs(dh), dm, abs(dm), dc, float(wt in ("P","G") or mt in ("P","G")),
            math.log10(afs["sas"]+1e-7), math.log10(afs["afr"]+1e-7), math.log10(afs["eas"]+1e-7),
            math.log10(afs["amr"]+1e-7), math.log10(afs["mid"]+1e-7), math.log10(afs["nfe"]+1e-7),
            math.log10(max(afs["asj"], afs["fin"])+1e-7), math.log10(grpmax+1e-7)
        ])
        y_list.append(int(r.get("label", r.get("y", 0))))
        clu_list.append(g2c.get(r["gene"], r["gene"]))
        grp_list.append(grp)
        am_raw.append(am)

    X, y, clu, grp, am_raw = np.array(X_list), np.array(y_list), np.array(clu_list), np.array(grp_list), np.array(am_raw)
    gkf = GroupKFold(n_splits=5)
    oof_raw = np.zeros(len(y))
    oof_cal = np.zeros(len(y))
    oof_sets = [None] * len(y)

    for tr, te in gkf.split(X, y, clu):
        # True Nested Inner-GroupKFold(3) OOF predictions for zero-overfit Platt + Mondrian calibration
        inner_oof = np.zeros(len(tr))
        for itr, iva in GroupKFold(n_splits=3).split(X[tr], y[tr], clu[tr]):
            isc = StandardScaler().fit(X[tr][itr])
            Xit_s, Xiv_s = isc.transform(X[tr][itr]), isc.transform(X[tr][iva])
            imlp = MLPClassifier(hidden_layer_sizes=(64, 32), max_iter=300, random_state=42).fit(Xit_s, y[tr][itr])
            ihgb = HistGradientBoostingClassifier(max_depth=4, max_iter=150, random_state=42).fit(X[tr][itr], y[tr][itr])
            ilr  = LogisticRegression(max_iter=500).fit(Xit_s, y[tr][itr])
            inner_oof[iva] = 0.50*imlp.predict_proba(Xiv_s)[:,1] + 0.35*ihgb.predict_proba(X[tr][iva])[:,1] + 0.15*ilr.predict_proba(Xiv_s)[:,1]

        sc = StandardScaler().fit(X[tr])
        Xtr_s, Xte_s = sc.transform(X[tr]), sc.transform(X[te])
        mlp = MLPClassifier(hidden_layer_sizes=(64, 32), max_iter=300, random_state=42).fit(Xtr_s, y[tr])
        hgb = HistGradientBoostingClassifier(max_depth=4, max_iter=150, random_state=42).fit(X[tr], y[tr])
        lr  = LogisticRegression(max_iter=500).fit(Xtr_s, y[tr])
        p_te = 0.50*mlp.predict_proba(Xte_s)[:,1] + 0.35*hgb.predict_proba(X[te])[:,1] + 0.15*lr.predict_proba(Xte_s)[:,1]

        platt = LogisticRegression().fit(_logit(inner_oof).reshape(-1, 1), y[tr])
        c_cal = platt.predict_proba(_logit(inner_oof).reshape(-1, 1))[:, 1]
        c_te  = platt.predict_proba(_logit(p_te).reshape(-1, 1))[:, 1]
        oof_raw[te] = p_te
        oof_cal[te] = c_te
        mc = AncestryMondrianConformal(alpha=0.10).fit(c_cal, y[tr], grp[tr])
        ps = mc.predict_sets(c_te, grp[te])
        for idx, s in zip(te, ps):
            oof_sets[idx] = s

    overall = {
        "n_variants": int(len(y)),
        "oof_auroc": round(float(roc_auc_score(y, oof_cal)), 4),
        "oof_pr_auc": round(float(average_precision_score(y, oof_cal)), 4),
        "raw_am_ece": round(float(compute_ece(y, am_raw)), 4),
        "v20_platt_ece": round(float(compute_ece(y, oof_cal)), 4),
        "empty_prediction_sets": int(sum(1 for s in oof_sets if len(s) == 0))
    }
    overall["ece_reduction_pct"] = round(100.0 * (overall["raw_am_ece"] - overall["v20_platt_ece"]) / max(overall["raw_am_ece"], 1e-9), 1)

    subgroups = {}
    for gname in ("SAS", "AFR", "EAS", "AMR", "MID", "NFE", "BOTTLENECK_FOUNDER", "ULTRA_RARE"):
        m = (grp == gname)
        n = int(np.sum(m))
        if n == 0:
            continue
        r_ece = compute_ece(y[m], am_raw[m])
        c_ece = compute_ece(y[m], oof_cal[m])
        idxs = np.where(m)[0]
        cov = float(np.mean([y[i] in oof_sets[i] for i in idxs]))
        sing = float(np.mean([len(oof_sets[i]) == 1 for i in idxs]))
        subgroups[gname] = {
            "n": n, "raw_am_ece": round(r_ece, 4), "v20_ece": round(c_ece, 4),
            "ece_cut_pct": round(100.0*(r_ece-c_ece)/max(r_ece, 1e-9), 1),
            "mondrian_cov": round(cov, 4), "singleton_pct": round(100.0*sing, 1),
            "empty_sets": int(sum(1 for i in idxs if len(oof_sets[i]) == 0))
        }

    report = {"version": "2.0.0-pan-global", "overall": overall, "subgroups": subgroups}
    report["sha256"] = hashlib.sha256(json.dumps(report, sort_keys=True).encode()).hexdigest()
    open("artifacts/v20_global_7ancestry_report.json", "w").write(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
