import csv, json, math, os, signal, sqlite3, time
import numpy as np
import torch
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from transformers import AutoTokenizer, AutoModelForMaskedLM
from mondrian_conformal import AncestryMondrianConformal, compute_ece
from resource_logger import ResourceLogger

STOP = False
def _sig(s, f):
    global STOP
    STOP = True
signal.signal(signal.SIGTERM, _sig)
signal.signal(signal.SIGINT, _sig)

def crop_window(seq, pos_1idx, max_len=1020):
    if len(seq) <= max_len:
        return seq, pos_1idx
    half = max_len // 2
    start = max(0, (pos_1idx - 1) - half)
    end = min(len(seq), start + max_len)
    start = max(0, end - max_len)
    return seq[start:end], pos_1idx - start

def run_esm2_cache(rows, db_path="artifacts/esm2_pilot.sqlite", max_wall_min=40):
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("CREATE TABLE IF NOT EXISTS scores (variant_id TEXT PRIMARY KEY, llr REAL, emb_norm REAL)")
    done = {r[0]: (r[1], r[2]) for r in conn.execute("SELECT variant_id, llr, emb_norm FROM scores")}
    todo = [r for r in rows if r["variant_id"] not in done]

    if todo:
        print(f"Scoring {len(todo)} new variants with Frozen ESM-2 (8M)...")
        tok = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
        mod = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D")
        mod.eval()
        t0 = time.time()
        for idx, r in enumerate(todo, 1):
            if STOP or (time.time() - t0) > max_wall_min * 60:
                print(f"Checkpointed cleanly at {idx-1}/{len(todo)} variants. Re-run to resume.")
                break
            wseq, wpos = crop_window(r["seq"], int(r["pos"]))
            inp = tok(wseq, return_tensors="pt")
            inp["input_ids"][0, wpos] = tok.mask_token_id
            with torch.no_grad():
                out = mod(**inp, output_hidden_states=True)
                lp = torch.log_softmax(out.logits[0, wpos], dim=-1)
                llr = float(lp[tok.convert_tokens_to_ids(r["mt"])] - lp[tok.convert_tokens_to_ids(r["wt"])])
                enorm = float(out.hidden_states[-1][0, wpos].norm().item())
            conn.execute("INSERT OR REPLACE INTO scores VALUES (?, ?, ?)", (r["variant_id"], llr, enorm))
            done[r["variant_id"]] = (llr, enorm)
            if idx % 25 == 0:
                conn.commit()
                print(f"  Progress: {idx}/{len(todo)} variants scored...")
        conn.commit()
    conn.close()
    return done

def main():
    logger = ResourceLogger()
    with open("data/pilot_variants.csv") as f:
        rows = list(csv.DictReader(f))
    scores = run_esm2_cache(rows)

    valid = [r for r in rows if r["variant_id"] in scores]
    if len(valid) < 100:
        print(f"Only {len(valid)} variants scored so far; re-run script to finish scoring.")
        return

    # Strict gene-holdout split (zero gene leakage!)
    genes = sorted({r["gene"] for r in valid})
    n_train = max(1, int(0.60 * len(genes)))
    n_cal = max(1, int(0.20 * len(genes)))
    train_genes = set(genes[:n_train])
    cal_genes = set(genes[n_train:n_train+n_cal])
    test_genes = set(genes[n_train+n_cal:])
    assert not (train_genes & cal_genes) and not (train_genes & test_genes) and not (cal_genes & test_genes)

    def feat(r):
        llr, enorm = scores[r["variant_id"]]
        return [
            llr, enorm,
            math.log10(float(r["af_nfe"]) + 1e-7),
            math.log10(float(r["af_sas"]) + 1e-7),
            math.log10(float(r["af_afr"]) + 1e-7),
            math.log10(float(r["af_all"]) + 1e-7),
        ]

    split = {"train": [], "cal": [], "test": []}
    for r in valid:
        if r["gene"] in train_genes: split["train"].append(r)
        elif r["gene"] in cal_genes: split["cal"].append(r)
        else: split["test"].append(r)

    X_tr = np.array([feat(r) for r in split["train"]])
    y_tr = np.array([int(r["label"]) for r in split["train"]])
    X_ca = np.array([feat(r) for r in split["cal"]])
    y_ca = np.array([int(r["label"]) for r in split["cal"]])
    g_ca = np.array([r["ancestry_group"] for r in split["cal"]])
    X_te = np.array([feat(r) for r in split["test"]])
    y_te = np.array([int(r["label"]) for r in split["test"]])
    g_te = np.array([r["ancestry_group"] for r in split["test"]])

    scaler = StandardScaler().fit(X_tr)
    mlp = MLPClassifier(hidden_layer_sizes=(32, 16), max_iter=300, random_state=42)
    mlp.fit(scaler.transform(X_tr), y_tr)

    p_ca = mlp.predict_proba(scaler.transform(X_ca))[:, 1]
    p_te = mlp.predict_proba(scaler.transform(X_te))[:, 1]
    raw_llr_te = -X_te[:, 0]  # more negative LLR = more pathogenic

    mc = AncestryMondrianConformal(alpha=0.10).fit(p_ca, y_ca, g_ca)
    group_metrics = mc.evaluate_by_group(p_te, y_te, g_te)

    report = {
        "n_total_scored": len(valid),
        "genes_split": {"train": sorted(train_genes), "cal": sorted(cal_genes), "test": sorted(test_genes)},
        "baseline_raw_esm2_test_auroc": round(float(roc_auc_score(y_te, raw_llr_te)), 4),
        "mlp_test_auroc": round(float(roc_auc_score(y_te, p_te)), 4),
        "mlp_test_prauc": round(float(average_precision_score(y_te, p_te)), 4),
        "mlp_test_ece_overall": round(compute_ece(y_te, p_te), 4),
        "mondrian_conformal_by_ancestry": group_metrics,
    }
    out_json = "artifacts/pilot_report.json"
    with open(out_json, "w") as f:
        json.dump(report, f, indent=2)
    sha = logger.sha256_file(out_json)
    logger.log("pilot_benchmark", {"report_sha256": sha, "mlp_auroc": report["mlp_test_auroc"]})
    print(json.dumps(report, indent=2))
    print(f"Report SHA-256: {sha}")

if __name__ == "__main__":
    main()
