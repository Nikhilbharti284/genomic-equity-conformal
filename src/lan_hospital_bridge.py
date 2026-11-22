"""
Genomic Equity Conformal v4.0.0 — Global Air-Gapped Hospital LAN Bridge & Offline UI
Runs on 0.0.0.0:8000 with zero external dependencies (stdlib ThreadingHTTPServer).
Enables hospital LAN computers & EHR/LIMS (HL7 FHIR R4) to talk offline (--network=none).
"""
import argparse
import json
import os
import sqlite3
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_chain import append_audit_entry, verify_audit_chain
from fhir_exporter import build_fhir_r4_bundle
from trio_phasing import analyze_trio_and_compound_het

OFFLINE_DB = os.environ.get("GEC_OFFLINE_DB", "artifacts/v3_offline_clinical_ref.sqlite")
LEDGER_PATH = os.environ.get("GEC_AUDIT_CHAIN", "artifacts/v4_audit_chain.jsonl")
LAN_API_KEY = os.environ.get("GEC_LAN_API_KEY", "").strip()

OFFLINE_UI_HTML = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>GEC v4.0.0 — Offline Hospital LAN Workbench</title>
<style>
body{font-family:system-ui,-apple-system,sans-serif;background:#f1f5f9;color:#0f172a;margin:0;padding:20px}
.wrap{max-width:980px;margin:auto;background:#fff;border:1px solid #cbd5e1;border-radius:14px;padding:22px;box-shadow:0 4px 14px rgba(15,23,42,.06)}
h1{margin:0 0 6px;color:#0f4c5c;font-size:22px} .sub{color:#475569;font-size:13px;margin:0 0 16px}
.row{display:flex;gap:10px;flex-wrap:wrap;margin:10px 0}
button{background:#0f4c5c;color:#fff;border:0;border-radius:8px;padding:10px 14px;font-weight:700;cursor:pointer;font-size:13px}
button.sec{background:#059669} button.alt{background:#334155}
textarea{width:100%;height:140px;font-family:monospace;font-size:12px;padding:10px;border:1px solid #cbd5e1;border-radius:8px;box-sizing:border-box}
pre{background:#0f172a;color:#6ee7b7;padding:14px;border-radius:10px;overflow:auto;font-size:12px;max-height:360px}
.badge{display:inline-block;padding:3px 9px;border-radius:99px;font-size:11px;font-weight:700;background:#e0f2fe;color:#0369a1;margin-right:6px}
</style></head><body><div class="wrap">
<span class="badge">v4.0.0 HOSPITAL LAN</span><span class="badge">100% OFFLINE • ZERO EGRESS</span><span class="badge">HL7 FHIR R4 + TRIO + HASH-CHAIN</span>
<h1>Genomic Equity Conformal — Offline Hospital LAN Workbench</h1>
<p class="sub">Access from any computer on the hospital intranet. No internet, no cloud upload, no login. Not standalone SaMD.</p>
<textarea id="inp">[
  {"variant_id":"ATP7B_H1069Q","gene":"ATP7B","moi":"AR","proband_gt":"0/1","mother_gt":"0/1","father_gt":"0/0","final_tier":"Pathogenic_90pct","platt_calibrated_prob":0.94},
  {"variant_id":"ATP7B_R778L","gene":"ATP7B","moi":"AR","proband_gt":"0/1","mother_gt":"0/0","father_gt":"0/1","final_tier":"Pathogenic_90pct","platt_calibrated_prob":0.91},
  {"variant_id":"NOVELGENE_A1V","gene":"NOVELGENE","moi":"AD","proband_gt":"0/1","mother_gt":"0/0","father_gt":"0/0","final_tier":"VUS_Review_Required"}
]</textarea>
<div class="row">
  <button onclick="runTrio(false)">Run Trio + OOD Audit</button>
  <button class="sec" onclick="runTrio(true)">Export HL7 FHIR R4 Bundle</button>
  <button class="alt" onclick="verifyChain()">Verify SHA-256 Audit Chain</button>
</div>
<pre id="out">Ready. Click any button above to run an offline LAN hospital audit.</pre>
</div><script>
async function runTrio(wantFhir){
  const out=document.getElementById('out'); out.textContent='Running offline audit...';
  try{
    const records=JSON.parse(document.getElementById('inp').value);
    const url=wantFhir?'/fhir/DiagnosticReport':'/trio-audit';
    const body=wantFhir?{results:records}:{records:records};
    const r=await fetch(url,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
    out.textContent=JSON.stringify(await r.json(),null,2);
  }catch(e){out.textContent='Error: '+e.message;}
}
async function verifyChain(){
  const r=await fetch('/audit-chain/verify');
  document.getElementById('out').textContent=JSON.stringify(await r.json(),null,2);
}
</script></body></html>"""


def load_known_genes(db_path: str = OFFLINE_DB) -> tuple:
    if not os.path.exists(db_path):
        return set(), 0
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        tables = [r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()]
        genes, total = set(), 0
        for t in tables:
            cols = [c[1] for c in conn.execute(f"PRAGMA table_info({t})").fetchall()]
            if "gene" in cols:
                genes.update(r[0].upper() for r in conn.execute(f"SELECT DISTINCT gene FROM {t}") if r[0])
                total = max(total, conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0])
        conn.close()
        return genes, total
    except Exception:
        return set(), 0


KNOWN_GENES, OFFLINE_ROWS = load_known_genes()


def annotate_ood(records: list) -> list:
    out = []
    for r in records:
        g = (r.get("gene") or "").upper()
        in_panel = (not KNOWN_GENES) or (g in KNOWN_GENES)
        item = dict(r)
        item["ood_status"] = "IN_OFFLINE_PANEL" if in_panel else "OOD_NOT_IN_OFFLINE_PANEL"
        out.append(item)
    return out


class HospitalLANHandler(BaseHTTPRequestHandler):
    def _send_json(self, code: int, obj: dict):
        body = json.dumps(obj, indent=2).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _authorized(self) -> bool:
        if not LAN_API_KEY:
            return True
        return self.headers.get("X-GEC-API-Key", "") == LAN_API_KEY

    def do_GET(self):
        if self.path in ("/", "/ui"):
            body = OFFLINE_UI_HTML.encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
        elif self.path == "/health":
            self._send_json(200, {
                "status": "ok",
                "version": "4.0.0-hospital-lan",
                "offline_db_rows": OFFLINE_ROWS,
                "known_genes_count": len(KNOWN_GENES),
                "lan_auth_enabled": bool(LAN_API_KEY),
                "audit_chain": verify_audit_chain(LEDGER_PATH),
            })
        elif self.path == "/audit-chain/verify":
            self._send_json(200, verify_audit_chain(LEDGER_PATH))
        else:
            self._send_json(404, {"error": "Endpoints: GET /ui, /health, /audit-chain/verify | POST /trio-audit, /fhir/DiagnosticReport"})

    def do_POST(self):
        if not self._authorized():
            return self._send_json(401, {"error": "Unauthorized LAN request (missing or invalid X-GEC-API-Key)"})
        length = int(self.headers.get("Content-Length", 0))
        data = json.loads(self.rfile.read(length).decode("utf-8")) if length > 0 else {}

        if self.path == "/trio-audit":
            recs = annotate_ood(data.get("records", data if isinstance(data, list) else []))
            trio = analyze_trio_and_compound_het(recs)
            merged = []
            for r in recs:
                t = trio["per_variant"].get(r["variant_id"], {})
                merged.append({**r, "inheritance": t.get("inheritance"), "compound_het_status": t.get("compound_het_status")})
            entry = append_audit_entry("TRIO_LAN_AUDIT", {"n_variants": len(merged), "alerts": len(trio["compound_het_alerts"])}, LEDGER_PATH)
            self._send_json(200, {
                "version": "4.0.0",
                "results": merged,
                "compound_het_alerts": trio["compound_het_alerts"],
                "audit_chain_seq": entry["seq"],
                "audit_chain_hash": entry["entry_hash"],
            })
        elif self.path == "/fhir/DiagnosticReport":
            recs = annotate_ood(data.get("results", data.get("records", [])))
            trio = analyze_trio_and_compound_het(recs)
            for r in recs:
                t = trio["per_variant"].get(r.get("variant_id", ""), {})
                r["compound_het_status"] = t.get("compound_het_status", "None")
            bundle = build_fhir_r4_bundle({"results": recs}, patient_ref=data.get("patient_ref", "Patient/anonymous-local"))
            entry = append_audit_entry("FHIR_R4_EXPORT", {"bundle_id": bundle["id"], "entries": bundle["total"]}, LEDGER_PATH)
            bundle["meta"] = {"tag": [{"code": f"gec-audit-seq-{entry['seq']}", "display": entry["entry_hash"][:16]}]}
            self._send_json(200, bundle)
        else:
            self._send_json(404, {"error": "Use POST /trio-audit or POST /fhir/DiagnosticReport"})


def main():
    ap = argparse.ArgumentParser(description="Genomic Equity Conformal v4.0.0 Offline Hospital LAN Bridge")
    ap.add_argument("--host", default="0.0.0.0")
    ap.add_argument("--port", type=int, default=8000)
    args = ap.parse_args()
    print(f"GEC v4.0.0 Hospital LAN Bridge listening on http://{args.host}:{args.port}/ui (offline rows: {OFFLINE_ROWS})")
    ThreadingHTTPServer((args.host, args.port), HospitalLANHandler).serve_forever()


if __name__ == "__main__":
    main()
