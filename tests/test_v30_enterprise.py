import json, os, sys, tempfile, threading, time, urllib.error, urllib.request
from http.server import ThreadingHTTPServer

sys.path.insert(0, os.path.abspath("src"))
from offline_reference_store import OfflineReferenceStore
from vcf_patient_profiler import parse_sample_gt, classify_clinical_zygosity, infer_patient_ancestry_from_vcf
from audit_vcf import audit_vcf_offline
from api_server import ClinicalAPIHandler
from clinical_report_generator import generate_clinical_report

def test_offline_sqlite_wal_submillisecond_lookup():
    st = OfflineReferenceStore()
    info = st.stats()
    assert info["n_records"] >= 9000
    t0 = time.perf_counter()
    hit = st.lookup_by_hgvs("chr11:g.5227002T>A")
    dt_ms = (time.perf_counter() - t0) * 1000.0
    assert hit is not None and hit["variant_id"] == "HBB_E7V"
    assert dt_ms < 25.0

def test_sickle_cell_01_carrier_vs_11_affected_zygosity():
    het = parse_sample_gt("GT:DP:GQ", "0/1:80:99", "chr11")
    hom = parse_sample_gt("GT:DP:GQ", "1/1:80:99", "chr11")
    hemi = parse_sample_gt("GT:DP:GQ", "1:50:99", "chrX")
    assert het["zygosity"] == "HETEROZYGOUS_CARRIER"
    assert hom["zygosity"] == "HOMOZYGOUS"
    assert hemi["zygosity"] == "HEMIZYGOUS"
    assert classify_clinical_zygosity("AR", het["zygosity"]) == "AR_Carrier_Only"
    assert classify_clinical_zygosity("AR", hom["zygosity"]) == "AR_Biallelic_Affected"
    assert classify_clinical_zygosity("XL", hemi["zygosity"]) == "XL_Hemizygous_Affected"

def test_aim_maximum_likelihood_ancestry_inference():
    afr_obs = [{"af_afr": 0.047, "af_sas": 0.0006, "af_eas": 0.0, "af_nfe": 0.0001}] * 5
    eas_obs = [{"af_eas": 0.104, "af_sas": 0.01, "af_afr": 0.001, "af_nfe": 0.002}] * 5
    a1, p1 = infer_patient_ancestry_from_vcf(afr_obs)
    a2, p2 = infer_patient_ancestry_from_vcf(eas_obs)
    assert a1 == "AFR" and p1["AFR"] > 0.9
    assert a2 == "EAS" and p2["EAS"] > 0.9

def test_airgapped_rest_api_endpoints():
    srv = ThreadingHTTPServer(("127.0.0.1", 0), ClinicalAPIHandler)
    port = srv.server_address[1]
    th = threading.Thread(target=srv.serve_forever, daemon=True)
    th.start()
    try:
        base = f"http://127.0.0.1:{port}"
        h = json.loads(urllib.request.urlopen(f"{base}/health", timeout=5).read())
        assert h["status"] == "ok" and h["version"] == "3.0.0"

        spec = json.loads(urllib.request.urlopen(f"{base}/openapi.json", timeout=5).read())
        assert spec["openapi"] == "3.0.3" and "/audit-vcf" in spec["paths"]

        docs = urllib.request.urlopen(f"{base}/docs", timeout=5).read().decode()
        assert "Offline OpenAPI Docs" in docs

        # Valid POST /audit
        req = urllib.request.Request(f"{base}/audit", data=json.dumps({"gene":"HBB","wt":"E","pos":7,"mt":"V","llr":-4.5,"af_afr":0.047}).encode(), headers={"Content-Type":"application/json"})
        r1 = json.loads(urllib.request.urlopen(req, timeout=5).read())
        assert r1["n"] == 1 and "request_id" in r1

        # Invalid POST /audit -> HTTP 400
        bad = urllib.request.Request(f"{base}/audit", data=json.dumps({"gene":"HBB"}).encode(), headers={"Content-Type":"application/json"})
        try:
            urllib.request.urlopen(bad, timeout=5)
            assert False, "Expected HTTP 400"
        except urllib.error.HTTPError as e:
            assert e.code == 400

        # Valid POST /audit-vcf
        vcf_txt = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr11\t5227002\trs334\tT\tA\t.\tPASS\t.\tGT:DP\t0/1:60\nchr13\t48045719\trs116855232\tC\tT\t.\tPASS\t.\tGT:DP\t0/1:50\n"
        req2 = urllib.request.Request(f"{base}/audit-vcf", data=vcf_txt.encode(), headers={"Content-Type":"text/plain"})
        r2 = json.loads(urllib.request.urlopen(req2, timeout=5).read())
        assert r2["n_audited"] == 2 and r2["n_pgx_alerts"] >= 1
        assert r2["results"][0]["clinical_zygosity_status"] == "AR_Carrier_Only"
    finally:
        srv.shutdown()

def test_nabl_cap_eu_ai_act_clinical_html_report():
    with tempfile.TemporaryDirectory() as td:
        in_vcf = os.path.join(td, "p.vcf")
        in_json = os.path.join(td, "p.json")
        out_html = os.path.join(td, "report.html")
        open(in_vcf, "w").write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr11\t5227002\trs334\tT\tA\t.\tPASS\t.\tGT\t1/1\nchr13\t48045719\trs116855232\tC\tT\t.\tPASS\t.\tGT\t0/1\n")
        res = audit_vcf_offline(in_vcf, default_ancestry="auto")
        open(in_json, "w").write(json.dumps(res))
        _, sha = generate_clinical_report(in_json, out_html)
        html_txt = open(out_html).read()
        assert len(sha) == 64
        assert "NABL ISO 15189" in html_txt and "EU AI ACT ANNEX IV" in html_txt
        assert "AR_Biallelic_Affected" in html_txt and "Thiopurines" in html_txt
