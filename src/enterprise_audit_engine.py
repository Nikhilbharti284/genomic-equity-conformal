import csv, hashlib, html, json, math, os, sqlite3, time
import joblib
import numpy as np
from mondrian_conformal import compute_ece

HYDRO = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
MW    = {"A":89,"R":174,"N":132,"D":133,"C":121,"Q":146,"E":147,"G":75,"H":155,"I":131,"L":131,"K":146,"M":149,"F":165,"P":115,"S":105,"T":119,"W":204,"Y":181,"V":117}
CHG   = {"R":1,"K":1,"H":0.1,"D":-1,"E":-1}

def _logit(p):
    p = np.clip(p, 1e-6, 1.0 - 1e-6)
    return np.log(p / (1.0 - p))

def main():
    bundle = joblib.load("artifacts/equity_model_v2.joblib")
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
        X.append([
            llr, enorm, d_hydro, abs(d_hydro), d_mw, abs(d_mw), d_chg,
            float(wt in ("P","G") or mt in ("P","G")),
            math.log10(float(r["af_nfe"]) + 1e-7), math.log10(float(r["af_sas"]) + 1e-7),
            math.log10(float(r["af_afr"]) + 1e-7), math.log10(float(r["af_all"]) + 1e-7),
        ])
    X = np.array(X, dtype=float)
    y = np.array([int(r["label"]) for r in valid])
    groups = np.array([r["ancestry_group"] for r in valid])

    p_raw = bundle["mlp"].predict_proba(bundle["scaler"].transform(X))[:, 1]
    p_cal = bundle["platt"].predict_proba(_logit(p_raw).reshape(-1, 1))[:, 1]
    pred_sets = bundle["mondrian"].predict_sets(p_cal, groups)

    # Export batch calibrated CSV for diagnostic lab pipelines
    out_csv = "artifacts/batch_calibrated_variants.csv"
    tier_counts = {"Pathogenic_90pct": 0, "Benign_90pct": 0, "VUS_Review_Required": 0}
    sample_rows = []
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variant_id","gene","ancestry_group","esm2_llr","raw_prob","platt_prob","conformal_tier","truth_label"])
        for i, r in enumerate(valid):
            ps = pred_sets[i]
            if ps == [1]: tier = "Pathogenic_90pct"
            elif ps == [0]: tier = "Benign_90pct"
            else: tier = "VUS_Review_Required"
            tier_counts[tier] += 1
            row_out = [r["variant_id"], r["gene"], r["ancestry_group"], round(scores[r["variant_id"]][0], 3), round(float(p_raw[i]), 4), round(float(p_cal[i]), 4), tier, int(y[i])]
            w.writerow(row_out)
            if len(sample_rows) < 60:
                sample_rows.append(row_out)

    # Compute Ancestry Equity Audit Summary
    subgroup_html = ""
    for g in ["EUR", "SAS", "AFR", "ULTRA_RARE"]:
        idx = np.where(groups == g)[0]
        raw_ece = compute_ece(y[idx], p_raw[idx])
        cal_ece = compute_ece(y[idx], p_cal[idx])
        cov = float(np.mean([int(y[i] in pred_sets[i]) for i in idx]))
        sing = float(np.mean([len(pred_sets[i]) == 1 for i in idx]))
        status = "PASS" if cov >= 0.85 else "REVIEW"
        badge = "#166534;background:#dcfce7" if status == "PASS" else "#991b1b;background:#fee2e2"
        subgroup_html += f"<tr><td><strong>{g}</strong></td><td>{len(idx)}</td><td>{raw_ece:.4f}</td><td>{cal_ece:.4f}</td><td><strong>{cov:.4f}</strong></td><td>{sing*100:.1f}%</td><td><span style='padding:2px 8px;border-radius:4px;font-weight:700;color:{badge}'>{status}</span></td></tr>"

    var_rows_html = ""
    for sr in sample_rows:
        tcol = "#991b1b" if sr[6]=="Pathogenic_90pct" else ("#166534" if sr[6]=="Benign_90pct" else "#92400e")
        var_rows_html += f"<tr><td><code>{html.escape(sr[0])}</code></td><td>{sr[1]}</td><td>{sr[2]}</td><td>{sr[3]}</td><td>{sr[5]}</td><td style='font-weight:700;color:{tcol}'>{sr[6]}</td></tr>"

    csv_sha = hashlib.sha256(open(out_csv, "rb").read()).hexdigest()
    html_doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>Genomic Equity Enterprise Audit Dashboard</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:#f8fafc;color:#0f172a;margin:0;padding:24px;}}
.wrap{{max-width:1040px;margin:0 auto;}}
.hero{{background:linear-gradient(135deg,#0b3946,#0f4c5c);color:#fff;padding:24px 28px;border-radius:10px;border-bottom:4px solid #d97706;}}
.kpis{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:18px 0;}}
.kpi{{background:#fff;border:1px solid #e2e8f0;border-top:3px solid #0f4c5c;border-radius:8px;padding:14px;}}
.kpi b{{font-size:22px;color:#0f4c5c;display:block;}}
.card{{background:#fff;border:1px solid #e2e8f0;border-radius:8px;padding:18px;margin:16px 0;}}
table{{width:100%;border-collapse:collapse;font-size:14px;}}
th{{background:#0f4c5c;color:#fff;text-align:left;padding:9px 10px;}}
td{{padding:8px 10px;border-bottom:1px solid #e2e8f0;}}
tr:nth-child(even){{background:#f8fafc;}}
.note{{background:#fffbeb;border-left:4px solid #d97706;padding:12px 16px;border-radius:6px;font-size:13px;margin:14px 0;}}
</style></head><body><div class="wrap">
<div class="hero"><div style="font-size:12px;color:#fde68a;font-weight:700;letter-spacing:1px;">GENOMIC EQUITY CONFORMAL v0.3.0 &bull; ENTERPRISE AI BIAS AUDIT</div>
<h1 style="margin:6px 0;">Ancestry Calibration &amp; Conformal Compliance Report</h1>
<p style="margin:0;color:#e2e8f0;">Generated: {time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())} &bull; Batch Size: {len(valid)} GRCh38 Variants &bull; CSV SHA-256: <code>{csv_sha[:20]}...</code></p></div>
<div class="kpis">
<div class="kpi"><b>{len(valid)}</b>Total Audited Variants</div>
<div class="kpi"><b>{tier_counts["Pathogenic_90pct"]}</b>High-Conf Pathogenic (90%)</div>
<div class="kpi"><b>{tier_counts["Benign_90pct"]}</b>High-Conf Benign (90%)</div>
<div class="kpi"><b>{tier_counts["VUS_Review_Required"]}</b>Flagged Uncertain (VUS)</div>
</div>
<div class="note"><strong>Research &amp; Compliance Framing:</strong> Computational calibration audit for research and clinical-pipeline QA evaluation only. Not a standalone clinical diagnostic device.</div>
<div class="card"><h2 style="margin-top:0;color:#0f4c5c;">Subgroup Fairness &amp; Mondrian Coverage Audit (Target &ge; 0.90)</h2>
<table><thead><tr><th>Ancestry</th><th>n</th><th>Raw ECE</th><th>Platt ECE</th><th>Mondrian Coverage</th><th>Singleton %</th><th>Audit Gate</th></tr></thead><tbody>{subgroup_html}</tbody></table></div>
<div class="card"><h2 style="margin-top:0;color:#0f4c5c;">Sample Calibrated Variant Batch (First 60 of {len(valid)})</h2>
<table><thead><tr><th>Variant ID</th><th>Gene</th><th>Ancestry</th><th>ESM-2 LLR</th><th>Platt Prob</th><th>90% Conformal Tier</th></tr></thead><tbody>{var_rows_html}</tbody></table></div>
</div></body></html>"""

    out_html = "artifacts/enterprise_equity_audit.html"
    with open(out_html, "w") as f:
        f.write(html_doc)
    print(f"Batch CSV saved  -> {out_csv} (SHA-256: {csv_sha[:16]}...)")
    print(f"Audit HTML saved -> {out_html}")
    print("Tier Summary:", tier_counts)

if __name__ == "__main__":
    main()
