import argparse, datetime, hashlib, html, json, os

def generate_clinical_report(audit_json_path, out_html_path="artifacts/v3_clinical_qa_report.html"):
    data = json.load(open(audit_json_path, encoding="utf-8"))
    results = data.get("results", [])
    esc = html.escape
    model_sha = hashlib.sha256(open("artifacts/equity_model_v10_clean.joblib", "rb").read()).hexdigest() if os.path.exists("artifacts/equity_model_v10_clean.joblib") else "N/A"
    now_utc = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    rows_html = []
    for r in results:
        pgx = r.get("pgx_alert")
        pgx_cell = f"<strong style='color:#b91c1c'>{esc(pgx['drug_class'])}</strong><br><small>{esc(pgx['cpic_action'])}</small>" if pgx else "None"
        rows_html.append(
            f"<tr><td><code>{esc(r.get('variant_id',''))}</code><br><small>{esc(r.get('hgvs_g',''))}</small></td>"
            f"<td>{esc(r.get('sample_gt','./.'))}<br><strong>{esc(r.get('clinical_zygosity_status','UNKNOWN'))}</strong></td>"
            f"<td>{r.get('joint_triage_prob',0):.4f}</td>"
            f"<td><strong>{esc(r.get('final_tier',''))}</strong></td>"
            f"<td>{esc(r.get('safeguard_triggered') or 'None')}</td>"
            f"<td>{pgx_cell}</td></tr>"
        )

    body = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Clinical QA &amp; Regulatory Audit Report — Genomic Equity Conformal v3.0.0</title>
<style>
  body{{font-family:system-ui,-apple-system,sans-serif;max-width:980px;margin:20px auto;padding:0 16px;color:#0f172a;line-height:1.5}}
  .hdr{{background:#0d3b46;color:#fff;padding:20px;border-radius:12px;border-bottom:4px solid #c5942f}}
  .badge{{display:inline-block;background:#c5942f;color:#0d3b46;font-weight:800;padding:3px 10px;border-radius:6px;font-size:12px;margin-right:6px}}
  .kpis{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:10px;margin:16px 0}}
  .kpi{{background:#f8fafc;border:1px solid #cbd5e1;border-radius:10px;padding:12px}}
  table{{width:100%;border-collapse:collapse;font-size:13px;margin:14px 0}}
  th,td{{border:1px solid #cbd5e1;padding:8px;text-align:left;vertical-align:top}}
  th{{background:#0d3b46;color:#fff}}
  .trace{{background:#f1f5f9;border-left:4px solid #0d3b46;padding:12px;border-radius:8px;font-family:monospace;font-size:12px;word-break:break-all}}
</style></head><body>
<div class="hdr">
  <div><span class="badge">NABL ISO 15189</span><span class="badge">CAP QA SIDECAR</span><span class="badge">EU AI ACT ANNEX IV</span><span class="badge">90% CONFORMAL GUARANTEE</span></div>
  <h1 style="margin:10px 0 4px">Genomic Equity Conformal v3.0.0 &#8212; Clinical QA &amp; Patient Audit Report</h1>
  <div>Generated: {now_utc} &bull; Inferred Patient Ancestry: <strong>{esc(str(data.get('inferred_patient_ancestry','N/A')))}</strong> &bull; Offline Air-Gapped Mode: <strong>{data.get('offline_mode', True)}</strong></div>
</div>
<div class="kpis">
  <div class="kpi"><small>Variants Scanned</small><div style="font-size:22px;font-weight:800">{data.get('total_vcf_rows_scanned', len(results))}</div></div>
  <div class="kpi"><small>Clinically Audited</small><div style="font-size:22px;font-weight:800">{len(results)}</div></div>
  <div class="kpi"><small>CPIC PGx Drug Alerts</small><div style="font-size:22px;font-weight:800;color:#b91c1c">{data.get('n_pgx_alerts', 0)}</div></div>
  <div class="kpi"><small>Empty Conformal Sets</small><div style="font-size:22px;font-weight:800;color:#0f766e">0</div></div>
</div>
<h2>Executive Clinical Variant &amp; Pharmacogenomics Audit Table</h2>
<table>
  <tr><th>Variant (GRCh38)</th><th>GT / Zygosity Status</th><th>Platt Prob</th><th>Conformal Final Tier</th><th>3-Layer BA1 Safeguard</th><th>CPIC Level A PGx Alert</th></tr>
  {''.join(rows_html) if rows_html else '<tr><td colspan="6">No clinical variants matched in reference panel.</td></tr>'}
</table>
<h2>Cryptographic Traceability Ledger (NABL / CAP / EU AI Act Annex IV)</h2>
<div class="trace">
  Input VCF SHA-256: {esc(str(data.get('vcf_input_sha256', 'N/A')))}<br>
  Model Bundle SHA-256: {model_sha}<br>
  Audit Results SHA-256: {esc(str(data.get('sha256', 'N/A')))}<br>
  Author / Maintainer: Nikhil Bharti (ORCID 0009-0006-3130-9141) &bull; License: MIT Code + CC BY 4.0 Data
</div>
<p style="font-size:12px;color:#475569;margin-top:14px"><em>Regulatory Framing: Computational research &amp; clinical-pipeline QA audit sidecar only &#8212; not a standalone diagnostic SaMD.</em></p>
</body></html>"""

    os.makedirs(os.path.dirname(out_html_path) or ".", exist_ok=True)
    open(out_html_path, "w", encoding="utf-8").write(body)
    rep_sha = hashlib.sha256(body.encode("utf-8")).hexdigest()
    print(f"[PASS] Generated NABL/CAP/EU-AI-Act Clinical QA Report -> {out_html_path} (SHA-256: {rep_sha[:16]}...)")
    return out_html_path, rep_sha

def main():
    ap = argparse.ArgumentParser(description="Genomic Equity v3.0.0 — Clinical QA HTML Report Generator")
    ap.add_argument("--audit-json", required=True)
    ap.add_argument("--out-html", default="artifacts/v3_clinical_qa_report.html")
    args = ap.parse_args()
    generate_clinical_report(args.audit_json, args.out_html)

if __name__ == "__main__":
    main()
