# Genomic Equity Conformal (v2.0.0)

[[https://github.com/Nikhilbharti284/genomic-equity-conformal/actions/workflows/ci.yml/badge.svg]](https://github.com/Nikhilbharti284/genomic-equity-conformal/actions/workflows/ci.yml)
[[https://img.shields.io/badge/ORCID-0009--0006--3130--9141-a6ce39?logo=orcid&logoColor=white]](https://orcid.org/0009-0006-3130-9141)
[[https://img.shields.io/badge/Code-MIT-0d3b46.svg]](LICENSE)
[[https://img.shields.io/badge/Data-CC%20BY%204.0-c5942f.svg]](https://creativecommons.org/licenses/by/4.0/)
[[https://img.shields.io/badge/Live%20Portal-GitHub%20Pages-0f766e.svg]](https://nikhilbharti284.github.io/genomic-equity-conformal/)

**Pan-Global 7-Ancestry Clinical Variant Calibration, CPIC Pharmacogenomics (`PGx`) & Autonomous Founder-Variant Discovery Engine** — air-gapped, zero-GPU, built & live-verified on $0-cost GitHub Codespaces Free Tier.

Created & maintained by **Nikhil Bharti** ([@Nikhilbharti284](https://github.com/Nikhilbharti284) | [ORCID 0009-0006-3130-9141](https://orcid.org/0009-0006-3130-9141)).

## Live-Verified v2.0.0 Pan-Global 7-Ancestry Benchmarks (118 Genes, n = 9,142)

- **OOF AUROC:** `0.9844` | **OOF PR-AUC:** `0.9898` | **Overall ECE:** `0.0686 -> 0.0066` (**-90.4%**) | **Empty Sets:** `0`
- **SHA-256 (`v20_global_7ancestry_report.json`):** `08b002cc58beba804c64416f2e4dc38ca2d266a55ee156a3af830f209844e08a`

| Cohort | n | Raw AM ECE | v2.0.0 ECE | ECE Cut | Mondrian Cov | Singletons | Empty Sets |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **SAS** | 700 | 0.0530 | **0.0229** | -56.8% | 0.9114 | 97.9% | 0 |
| **AFR** | 1,221 | 0.0543 | **0.0141** | -74.0% | 0.9328 | 100.0% | 0 |
| **EAS** | 740 | 0.0685 | **0.0335** | -51.1% | 0.9149 | 95.8% | 0 |
| **AMR** | 717 | 0.0480 | **0.0288** | -40.0% | 0.9121 | 94.3% | 0 |
| **MID** | 122 | 0.1127 | **0.0566** | -49.8% | 0.8934 | 100.0% | 0 |
| **NFE** | 1,192 | 0.0941 | **0.0180** | -80.9% | 0.9077 | 90.5% | 0 |
| **BOTTLENECK_FOUNDER (ASJ+FIN)** | 622 | 0.0669 | **0.0179** | -73.3% | 0.9277 | 100.0% | 0 |
| **ULTRA_RARE** | 3,828 | 0.1333 | **0.0054** | -96.0% | 0.9875 | 100.0% | 0 |

## What's New in v2.0.0
1. **Autonomous Unlisted Founder Discovery (`src/discover_unlisted_ba1_candidates.py`):** Discovered **19 unlisted BS1+ pathogenic founder/balancing-selection variants** (`ABCA4_G1961E` ASJ 2.39%, `ATP7B_H1069Q` ASJ, `SLC12A3_T180K` EAS, `SLC26A4_F335L` & `SLC26A4_V239D` SAS, etc.) saved in `artifacts/v20_novel_founder_discoveries.json` & `.csv`.
2. **Global CPIC Pharmacogenomics Sidecar (`src/pgx_drug_safety.py`):** Drug-toxicity alerts for `DPYD`, `NUDT15`, `TPMT`, `CYP2C19`, `SLCO1B1`, `G6PD`, and `HLA-B*57:01`/`B*15:02`, written directly into annotated VCF 4.2 (`GEC_PGX_ALERT`).
3. **AlphaFold DB v6 C-Beta RSA Burial Upgrade (`src/v2_mavedb_rsa_upgrade.py`):** Flipped `TP53` MaveDB DMS from miss to **WIN** (`0.2485 -> 0.2488`), alongside wins on `BRCA1`, `PTEN`, and `MSH2` (`LDLR` `-0.0002` honestly reported).

## Quickstart

    pip install -e . --no-deps --no-build-isolation
    gec-audit --vcf path/to/hospital.vcf --out-vcf artifacts/audited_output.vcf
    pytest -v tests/

## Disclaimer
Computational research and clinical-pipeline QA audit sidecar only — not a standalone diagnostic SaMD.
