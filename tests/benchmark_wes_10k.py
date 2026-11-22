import gzip, json, os, resource, sys, tempfile, time
sys.path.insert(0, os.path.abspath("src"))
from audit_vcf import audit_vcf_offline
from inference_engine import GenomicEquityEngine
from offline_reference_store import OfflineReferenceStore
from clinical_report_generator import generate_clinical_report

def run_benchmark():
    eng = GenomicEquityEngine()
    store = OfflineReferenceStore()
    with tempfile.TemporaryDirectory() as td:
        vcf_gz = os.path.join(td, "synthetic_wes_10k.vcf.gz")
        out_vcf = os.path.join(td, "audited_wes_10k.vcf")
        out_json = "artifacts/v3_wes_sample_audit.json"
        with gzip.open(vcf_gz, "wt", encoding="utf-8") as f:
            f.write("##fileformat=VCFv4.2\n##reference=GRCh38\n")
            f.write("\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE1"]) + "\n")
            # 4 high-impact clinical/PGx control variants
            f.write("\t".join(["chr17","7675088","rs28934578","C","T",".","PASS",".","GT:DP:GQ","0/1:85:99"]) + "\n")
            f.write("\t".join(["chr11","5227002","rs334","T","A",".","PASS",".","GT:DP:GQ","0/1:92:99"]) + "\n")
            f.write("\t".join(["chr17","7676154","rs1042522","G","C",".","PASS",".","GT:DP:GQ","1/1:110:99"]) + "\n")
            f.write("\t".join(["chr13","48045719","rs116855232","C","T",".","PASS",".","GT:DP:GQ","0/1:74:99"]) + "\n")
            # 9,996 background exome variants to reach 10,000 rows
            for i in range(9996):
                f.write(f"chr1\t{100000 + i}\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t0/1:40:99\n")

        t0 = time.perf_counter()
        res = audit_vcf_offline(vcf_gz, out_vcf=out_vcf, default_ancestry="auto", eng=eng, store=store)
        elapsed = time.perf_counter() - t0
        peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
        open(out_json, "w", encoding="utf-8").write(json.dumps(res, indent=2) + "\n")
        rep_path, rep_sha = generate_clinical_report(out_json, "artifacts/v3_clinical_qa_report.html")

        print(f"[BENCHMARK 10K WES] Scanned={res['total_vcf_rows_scanned']} | Audited={res['n_audited']} | PGx={res['n_pgx_alerts']} | InferredAnc={res['inferred_patient_ancestry']} | Time={elapsed:.3f}s | PeakRAM={peak_mb:.1f} MB")
        assert res["total_vcf_rows_scanned"] == 10000
        assert res["n_audited"] == 4
        assert elapsed < 5.0, f"Too slow: {elapsed:.2f}s"
        assert peak_mb < 250.0, f"Peak RAM exceeded 250 MB: {peak_mb:.1f} MB"
        return elapsed, peak_mb

if __name__ == "__main__":
    run_benchmark()
