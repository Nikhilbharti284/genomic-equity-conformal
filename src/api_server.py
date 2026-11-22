import json, os, sys, tempfile, uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inference_engine import GenomicEquityEngine
from offline_reference_store import OfflineReferenceStore
from audit_vcf import audit_vcf_offline

OPENAPI_SPEC = {
    "openapi": "3.0.3",
    "info": {
        "title": "Genomic Equity Conformal v3.0.0 — Air-Gapped Clinical QA API",
        "version": "3.0.0",
        "description": "Zero-GPU, air-gapped 7-ancestry clinical variant calibration, 3-Layer BA1 Safeguard, zygosity gating & CPIC PGx sidecar."
    },
    "paths": {
        "/health": {"get": {"summary": "Offline store & model health check", "responses": {"200": {"description": "Healthy"}}}},
        "/openapi.json": {"get": {"summary": "OpenAPI 3.0.3 JSON schema", "responses": {"200": {"description": "Schema"}}}},
        "/docs": {"get": {"summary": "Interactive self-contained API explorer", "responses": {"200": {"description": "HTML UI"}}}},
        "/audit": {"post": {"summary": "Audit single or batch variants (JSON)", "responses": {"200": {"description": "Clinical QA cards"}, "400": {"description": "Validation error"}}}},
        "/audit-vcf": {"post": {"summary": "Audit raw VCF 4.2 text payload offline", "responses": {"200": {"description": "VCF audit JSON + PGx + zygosity"}}}}
    }
}

DOCS_HTML = """<!doctype html><html><head><meta charset="utf-8"><title>Genomic Equity Conformal v3.0.0 API</title>
<style>body{font-family:system-ui,sans-serif;max-width:820px;margin:24px auto;padding:0 16px;color:#0f172a}h1{color:#0d3b46}pre{background:#0b1f26;color:#a7f3d0;padding:12px;border-radius:8px;overflow:auto}</style></head>
<body><h1>Genomic Equity Conformal v3.0.0 — Offline OpenAPI Docs</h1>
<p>Air-gapped hospital REST API (zero external CDN required). Endpoints: <code>GET /health</code>, <code>GET /openapi.json</code>, <code>POST /audit</code>, <code>POST /audit-vcf</code>.</p>
<pre>""" + json.dumps(OPENAPI_SPEC, indent=2) + """</pre></body></html>"""

ENGINE = None
STORE = None

def _get_engine():
    global ENGINE, STORE
    if ENGINE is None:
        ENGINE = GenomicEquityEngine()
    if STORE is None:
        STORE = OfflineReferenceStore()
    return ENGINE, STORE

class ClinicalAPIHandler(BaseHTTPRequestHandler):
    def _send(self, code, body, ctype="application/json"):
        raw = body.encode("utf-8") if isinstance(body, str) else json.dumps(body).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)

    def log_message(self, format, *args):
        pass

    def do_GET(self):
        if self.path == "/health":
            _, st = _get_engine()
            return self._send(200, {"status": "ok", "version": "3.0.0", "offline_store": st.stats()})
        if self.path == "/openapi.json":
            return self._send(200, OPENAPI_SPEC)
        if self.path in ("/docs", "/"):
            return self._send(200, DOCS_HTML, "text/html; charset=utf-8")
        return self._send(404, {"error": "Not found"})

    def do_POST(self):
        req_id = str(uuid.uuid4())
        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length).decode("utf-8", errors="ignore")
        eng, st = _get_engine()

        if self.path == "/audit":
            try:
                data = json.loads(raw)
                items = data if isinstance(data, list) else (data.get("variants") if isinstance(data, dict) and "variants" in data else [data])
                out = []
                for it in items:
                    for req_k in ("gene", "wt", "pos", "mt", "llr"):
                        if req_k not in it:
                            return self._send(400, {"request_id": req_id, "error": f"Missing required field: {req_k}"})
                    out.append(eng.predict(
                        gene=str(it["gene"]), wt=str(it["wt"]), pos=int(it["pos"]), mt=str(it["mt"]),
                        llr=float(it["llr"]), enorm=float(it.get("enorm", 10.0)),
                        ancestry=str(it.get("ancestry", "SAS")),
                        af_nfe=float(it.get("af_nfe", 0.0)), af_sas=float(it.get("af_sas", 0.0)),
                        af_afr=float(it.get("af_afr", 0.0)), af_eas=float(it.get("af_eas", 0.0)),
                        af_amr=float(it.get("af_amr", 0.0)), af_mid=float(it.get("af_mid", 0.0)),
                        af_asj=float(it.get("af_asj", 0.0)), af_fin=float(it.get("af_fin", 0.0)),
                        af_all=float(it.get("af_all", 0.0)), rsid=it.get("rsid")
                    ))
                return self._send(200, {"request_id": req_id, "version": "3.0.0", "n": len(out), "results": out})
            except Exception as e:
                return self._send(400, {"request_id": req_id, "error": str(e)})

        if self.path == "/audit-vcf":
            if "#CHROM" not in raw:
                return self._send(400, {"request_id": req_id, "error": "Invalid VCF payload: missing #CHROM header"})
            with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False, encoding="utf-8") as tf:
                tf.write(raw)
                tmp_path = tf.name
            try:
                res = audit_vcf_offline(tmp_path, out_vcf=None, default_ancestry="auto", eng=eng, store=st)
                res["request_id"] = req_id
                return self._send(200, res)
            finally:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)

        return self._send(404, {"request_id": req_id, "error": "Not found"})

def start_server(host="127.0.0.1", port=8000):
    srv = ThreadingHTTPServer((host, port), ClinicalAPIHandler)
    print(f"Genomic Equity Conformal v3.0.0 API listening on http://{host}:{port}/docs")
    srv.serve_forever()

if __name__ == "__main__":
    start_server()
