import json
from http.server import BaseHTTPRequestHandler, HTTPServer
from inference_engine import GenomicEquityEngine

ENGINE = GenomicEquityEngine()

class EquityAPIHandler(BaseHTTPRequestHandler):
    def _send(self, code, obj):
        body = json.dumps(obj, indent=2).encode("utf-8")
        self.send_response(code); self.send_header("Content-Type", "application/json"); self.send_header("Content-Length", str(len(body))); self.end_headers(); self.wfile.write(body)

    def do_GET(self):
        if self.path == "/health":
            self._send(200, {"status": "ok", "version": "1.0.0-clean", "commercial_safe": True, "shared_engine": True})
        else:
            self._send(404, {"error": "Use GET /health or POST /audit"})

    def do_POST(self):
        if self.path != "/audit":
            return self._send(404, {"error": "Use POST /audit"})
        data = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))).decode("utf-8"))
        items = data if isinstance(data, list) else [data]
        res = [ENGINE.predict(
            gene=it["gene"], wt=it["wt"], pos=int(it["pos"]), mt=it["mt"], llr=float(it["esm2_llr"]),
            enorm=float(it.get("emb_norm", 10.0)), ancestry=it.get("ancestry", "SAS"),
            af_nfe=float(it.get("af_nfe", 0.0)), af_sas=float(it.get("af_sas", 0.0)),
            af_afr=float(it.get("af_afr", 0.0)), af_all=float(it.get("af_all", 0.0)),
        ) for it in items]
        self._send(200, {"n_audited": len(res), "results": res})

if __name__ == "__main__":
    print("Genomic Equity v1.0.0-clean On-Premise API running on http://127.0.0.1:8080")
    HTTPServer(("127.0.0.1", 8080), EquityAPIHandler).serve_forever()
