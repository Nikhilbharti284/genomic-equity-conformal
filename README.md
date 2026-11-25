# Genomic Equity Conformal (v5.0.0)

<p>
  <a href="https://github.com/Nikhilbharti284/genomic-equity-conformal/releases/tag/v5.0.0"><img src="https://img.shields.io/badge/Release-v5.0.0-111111?style=flat-square" alt="Release v5.0.0"></a>
  <a href="https://doi.org/10.5281/zenodo.24465257"><img src="https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.24465257-111111?style=flat-square" alt="Zenodo DOI"></a>
  <a href="https://orcid.org/0009-0006-3130-9141"><img src="https://img.shields.io/badge/ORCID-0009--0006--3130--9141-111111?style=flat-square" alt="ORCID"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/Code-MIT-111111?style=flat-square" alt="MIT License"></a>
  <a href="https://creativecommons.org/licenses/by/4.0/"><img src="https://img.shields.io/badge/Data-CC%20BY%204.0-555555?style=flat-square" alt="CC BY 4.0"></a>
  <a href="#quickstart--verification"><img src="https://img.shields.io/badge/Pytest-33%2F33%20Passing-111111?style=flat-square" alt="Pytest 33/33"></a>
  <a href="https://nikhilbharti284.github.io/genomic-equity-conformal/"><img src="https://img.shields.io/badge/Live%20Portal-GitHub%20Pages-111111?style=flat-square" alt="Live Portal"></a>
</p>

**100% Pan-Human Exome (`19,276` HGNC genes, `67,402` GRCh38 ClinVar >=2-star missense rows), Air-Gapped (`--network=none`) 7-Ancestry Mondrian Conformal Calibration, Multi-Sample Trio Compound-Het Phasing, CPIC PGx, HL7 FHIR R4 Export & Tamper-Evident Hospital LAN QA Sidecar** — zero-GPU, commercially clean (`MIT + CC BY 4.0`).

Created & maintained by **Nikhil Bharti** ([@Nikhilbharti284](https://github.com/Nikhilbharti284) | [ORCID 0009-0006-3130-9141](https://orcid.org/0009-0006-3130-9141) | [LinkedIn](https://www.linkedin.com/in/nikhil-bharti-063b9a225)).

> **Non-SaMD Notice:** Open-source computational research and clinical-pipeline QA audit sidecar — **not** a standalone diagnostic Software as a Medical Device (SaMD).

## What's New in v5.0.0 (100% Pan-Human Exome Edition)

1. **Pan-Human Exome SQLite WAL (`src/build_v5_exome_store.py`):** Indexes all **19,276 HGNC protein-coding genes**, preserves all `9,249` calibrated `v4` rows, and streams NCBI ClinVar >=2-star GRCh38 missense SNVs into `artifacts/v5_panhuman_exome_ref.sqlite` (**`67,402` rows across `7,336` unique genes**, **`15.75 MB`**, SHA-256 `d30f06f90be50f51aacd6e6ccc752031bf3951758524cd362d4fe8f292df9efb`).
2. **Two-Tier Hybrid OOD Fallback (`src/v5_hybrid_fallback.py`):** `TIER1_SQLITE_INDEXED`, `TIER2_CALIBRATED_NOVEL_MISSENSE` (zero-drop across all 19,276 HGNC genes), and `OOD_NON_MISSENSE_OR_UNMAPPED` abstention.
3. **Multi-Continental Benchmark (`src/v5_global_cohort_audit.py`):** Verified across `SAS, AFR, EAS, AMR, EUR` at **`184.79 MB` peak RSS RAM** (`<250 MB` edge target) with Trio trans compound-het verification.
4. **Full Verification (`tests/test_v5_exome_and_global.py`):** **`33/33` pytest passing (`4.58s`)** and W3C RO-Crate 1.1 (`ro-crate-metadata.json`, `71` SHA-256 verified files).

## Preserved Core Capabilities (v0.1.0 -> v4.0.0)

- **v4.0.0:** Hospital LAN bridge (`0.0.0.0:8000/ui`, optional `X-GEC-API-Key`), Multi-Sample Trio Compound-Het phasing, HL7 FHIR R4 Genomics Bundle (`LOINC 51969-4` / `69548-6`), SHA-256 tamper-evident hash chain, Zenodo DOI `10.5281/zenodo.24465257`.
- **v3.0.0:** 10,000 WES streaming (`0.917s`, `181.1 MB` RAM), GT zygosity, AIM 7-ancestry MLE, stdlib OpenAPI REST API, NABL/CAP/EU AI Act HTML QA reports.
- **v2.0.0:** 7-ancestry Non-Empty Mondrian Conformal (`0.9844` OOF AUROC, `0.9898` PR-AUC, `n=9,142`, ECE `-90.4%`, `0` empty sets), 19 unlisted BS1+ founders, 6-panel CPIC PGx, AlphaFold DB v6 C-Beta RSA (`4/5` MaveDB wins).
- **v0.1.0 -> v1.1.0:** 118-gene / 100-cluster homology core, `GenomicEquityEngine`, VCF 4.2 INFO tags, Pejaver PP3/BP4 tiers, 11 ClinGen BA1 guards (`HBB_E7V`).

## Quickstart & Verification

```bash
pip install -e . --no-deps --no-build-isolation
python3 src/build_v5_exome_store.py
python3 src/v5_global_cohort_audit.py
python3 -m pytest -v
docker build -t genomic-equity:5.0.0 .
docker run --network=none -p 8000:8000 genomic-equity:5.0.0


```
