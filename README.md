# Genomic Equity Conformal (v4.0.0)

[[https://github.com/Nikhilbharti284/genomic-equity-conformal/actions/workflows/ci.yml/badge.svg]](https://github.com/Nikhilbharti284/genomic-equity-conformal/actions/workflows/ci.yml)
[[https://img.shields.io/badge/Release-v4.0.0-0f4c5c?style=flat-square]](https://github.com/Nikhilbharti284/genomic-equity-conformal/releases/tag/v4.0.0)
[[https://img.shields.io/badge/ORCID-0009--0006--3130--9141-a6ce39?style=flat-square&logo=orcid&logoColor=white]](https://orcid.org/0009-0006-3130-9141)
[[https://img.shields.io/badge/Code-MIT-0d3b46?style=flat-square]](LICENSE)
[[https://img.shields.io/badge/Data-CC%20BY%204.0-c5942f?style=flat-square]](https://creativecommons.org/licenses/by/4.0/)
[[https://img.shields.io/badge/Pytest-27%2F27%20Passing-16a34a?style=flat-square]](#quickstart--verification)
[[https://img.shields.io/badge/Interop-HL7%20FHIR%20R4-0284c7?style=flat-square]](#whats-new-in-v400-global-hospital-lan--interoperability)
[[https://img.shields.io/badge/Live%20Portal-GitHub%20Pages-0f766e?style=flat-square]](https://nikhilbharti284.github.io/genomic-equity-conformal/)

**Air-Gapped Pan-Global 7-Ancestry Clinical Variant Calibration, Multi-Sample Trio Compound-Het Phasing, CPIC Pharmacogenomics (`PGx`), HL7 FHIR R4 Export & Tamper-Evident Hospital LAN QA Sidecar** — zero-GPU, 100% commercially clean (`MIT + CC BY 4.0`), built & live-verified on `$0`-cost GitHub Codespaces Free Tier.

Created & maintained by **Nikhil Bharti** ([@Nikhilbharti284](https://github.com/Nikhilbharti284) | [ORCID 0009-0006-3130-9141](https://orcid.org/0009-0006-3130-9141)).

## Why This Matters Worldwide

- **Zero Cloud Egress & Sovereign Privacy:** Runs **100% air-gapped (`--network=none`)** on a basic 2-core CPU workstation—compliant with **EU GDPR, UK NHS, US HIPAA, India DPDP Act, Brazil LGPD, Saudi PDPL, and South Africa POPIA** without uploading patient VCFs to external clouds.
- **Eliminates Global Ancestry Bias:** Cuts overall Expected Calibration Error by **−90.4%** (`0.0686 → 0.0066`) with **0 empty sets** across 8 global strata (`SAS, AFR, EAS, AMR, MID, NFE, ASJ+FIN, ULTRA_RARE`) while protecting balancing-selection alleles (`HBB_E7V`, `GJB2_V37I`, `G6PD_V68M`).
- **Plug-and-Play Hospital LAN Interoperability (`v4.0.0`):** Any doctor, geneticist, or LIMS/EHR computer on the hospital intranet can connect offline (`http://<server-ip>:8000/ui`) for **Multi-Sample Trio Compound-Het Phasing**, **OOD Panel Abstention**, **HL7 FHIR R4 Genomics JSON**, and **SHA-256 Hash-Chained Audit Verification**.

## What's New in v4.0.0 (Global Hospital LAN & Interoperability)

1. **Offline Hospital LAN Bridge & Zero-CDN Web Workbench (`src/lan_hospital_bridge.py`):** Binds `0.0.0.0:8000` with optional offline `X-GEC-API-Key` (`GEC_LAN_API_KEY`), serves a built-in zero-login `/ui` dashboard for any computer on the hospital intranet, and flags unseen genes with explicit **`OOD_NOT_IN_OFFLINE_PANEL`** abstention.
2. **Multi-Sample Trio & Compound-Het Phasing Engine (`src/trio_phasing.py`):** Parses Proband + Mother + Father VCF genotypes (or phased `0|1` / `1|0`) to detect **`AR_Compound_Het_Trans_Affected`**, **`AR_Cis_Carrier_Only`**, **`AR_Unphased_Potential_Compound_Het`**, and **`De_Novo_Candidate`** hits.
3. **Global HL7 FHIR R4 Genomics Exporter (`src/fhir_exporter.py`):** Exposes `POST /fhir/DiagnosticReport` emitting standard FHIR R4 `Bundle` (`DiagnosticReport` LOINC `51969-4` + `Observation` LOINC `69548-6`) ready for Epic, Cerner, Bahmni, OpenMRS, UK NHS, EU, and ABDM LIMS.
4. **Tamper-Evident SHA-256 Hash-Chain Audit Ledger (`src/audit_chain.py`):** Appends every LAN/FHIR audit to `artifacts/v4_audit_chain.jsonl` with cryptographic `prev_hash → entry_hash` chaining and live `GET /audit-chain/verify` (verified to catch even 1-byte tampering).
5. **Full Regression Coverage (`tests/test_v4_hospital_lan.py`):** **`27/27` pytest passing in `37.42s`** (preserving 100% of `v3.0.0` tests + 7 new `v4.0.0` tests).

## What's New in v3.0.0 (Enterprise Air-Gapped Clinical Core — Fully Preserved)

1. **Zero-Internet Offline SQLite WAL Store (`src/offline_reference_store.py`):** `9,249` indexed GRCh38 records (`artifacts/v3_offline_clinical_ref.sqlite`) for sub-ms `--network=none` execution (`gec-audit --offline`).
2. **Sample Zygosity (`GT`/`DP`/`GQ`) & AIM 7-Ancestry MLE (`src/vcf_patient_profiler.py`):** Separates `0/1` `AR_Carrier_Only` vs `1/1` `AR_Biallelic_Affected` (e.g., `HBB_E7V`) and auto-infers 7-ancestry background (`--default-ancestry auto`).
3. **Constant-Memory Streaming WES Auditor & Stdlib OpenAPI Server (`src/audit_vcf.py`, `src/api_server.py`):** Benchmarked `10,000`-variant WES streaming in **`0.917s`** at **`181.1 MB` peak RAM**, emitting 7 `GEC_*` VCF 4.2 INFO tags (`GEC_PGX_ALERT` included).
4. **NABL (ISO 15189) / CAP / EU AI Act Annex IV HTML QA Report (`src/clinical_report_generator.py`):** Offline print-ready patient QA report (`gec-report`) with SHA-256 traceability.
5. **Autonomous Discoveries, CPIC PGx & AlphaFold DB v6 3D Structure:** **19 unlisted BS1+ founder candidates** (`ABCA4_G1961E` ASJ 2.39%, `ATP7B_H1069Q` ASJ, `SLC12A3_T180K` EAS, `SLC26A4_F335L`/`V239D` SAS), **6-panel CPIC Level A PGx** (`DPYD`, `NUDT15/TPMT`, `CYP2C19`, `SLCO1B1`, `G6PD`, `HLA-B`), and pure-NumPy **10Å C-Beta RSA Burial** beating raw `ESM2-35M` on **4/5 wet-lab MaveDB genes** (`TP53`, `BRCA1`, `PTEN`, `MSH2`; `LDLR` `−0.0002` honestly disclosed).

## Live-Verified Pan-Global 7-Ancestry Benchmarks (118 Genes, 100 Clusters, n = 9,142)

- **OOF AUROC:** `0.9844` | **OOF PR-AUC:** `0.9898` | **Overall ECE:** `0.0686 → 0.0066` (**−90.4%**) | **Empty Sets:** `0`
- **SHA-256 (`artifacts/v20_global_7ancestry_report.json`):** `08b002cc58beba804c64416f2e4dc38ca2d266a55ee156a3af830f209844e08a`

| Cohort | n | Raw AM ECE | Calibrated ECE | ECE Cut | Mondrian Cov | Singletons | Empty Sets |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **Overall (All 118 Genes)** | 9,142 | 0.0686 | **0.0066** | **−90.4%** | ≥ 0.90 | — | 0 |
| **SAS (South Asian)** | 700 | 0.0530 | **0.0229** | −56.8% | 0.9114 | 97.9% | 0 |
| **AFR (African / African Am.)** | 1,221 | 0.0543 | **0.0141** | −74.0% | 0.9328 | 100.0% | 0 |
| **EAS (East Asian)** | 740 | 0.0685 | **0.0335** | −51.1% | 0.9149 | 95.8% | 0 |
| **AMR (Admixed American)** | 717 | 0.0480 | **0.0288** | −40.0% | 0.9121 | 94.3% | 0 |
| **MID (Middle Eastern)\*** | 122 | 0.1127 | **0.0566** | −49.8% | 0.8934\* | 100.0% | 0 |
| **NFE (Non-Finnish European)** | 1,192 | 0.0941 | **0.0180** | −80.9% | 0.9077 | 90.5% | 0 |
| **BOTTLENECK_FOUNDER (ASJ+FIN)** | 622 | 0.0669 | **0.0179** | −73.3% | 0.9277 | 100.0% | 0 |
| **ULTRA_RARE (AF = 0)** | 3,828 | 0.1333 | **0.0054** | −96.0% | 0.9875 | 100.0% | 0 |

*\*Honest finite-sample disclosure: `MID` has `n = 122` (`0.8934` is within 1 variant of `0.90`). All other 7 strata exceed `0.90`.*

## Quickstart & Verification

```bash
# 1. Install editable package (zero-GPU CPU)
pip install -e . --no-deps --no-build-isolation

# 2. Run full v3 + v4 test suite (Verified: 27 passed in 37.42s)
python3 -m pytest -q tests/

# 3. Launch Offline Hospital LAN Bridge (open http://<server-ip>:8000/ui from any LAN computer)
export GEC_LAN_API_KEY="hospital-local-secret"   # Optional LAN auth
python3 src/lan_hospital_bridge.py --host 0.0.0.0 --port 8000

# 4. Run air-gapped CLI WES audit & NABL/CAP/EU AI Act HTML report
gec-audit --offline --default-ancestry auto --vcf sample.vcf.gz --out-vcf audited.vcf
gec-report --audit-json artifacts/raw_vcf_audit_output.json --out-html report.html

# 5. Zero-Egress Hospital Docker
docker build -t genomic-equity:4.0.0 .
docker run --network=none -p 8000:8000 genomic-equity:4.0.0
