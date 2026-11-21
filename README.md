# Genomic Equity Conformal (v1.1.0)

[[https://github.com/Nikhilbharti284/genomic-equity-conformal/actions/workflows/ci.yml/badge.svg]](https://github.com/Nikhilbharti284/genomic-equity-conformal/actions/workflows/ci.yml)
[[https://img.shields.io/badge/ORCID-0009--0006--3130--9141-a6ce39?logo=orcid&logoColor=white]](https://orcid.org/0009-0006-3130-9141)
[[https://img.shields.io/badge/Code-MIT-0d3b46.svg]](LICENSE)
[[https://img.shields.io/badge/Data-CC%20BY%204.0-c5942f.svg]](https://creativecommons.org/licenses/by/4.0/)
[[https://img.shields.io/badge/Live%20Portal-GitHub%20Pages-0f766e.svg]](https://nikhilbharti284.github.io/genomic-equity-conformal/)

**Air-Gapped, Zero-GPU Clinical Variant Calibration, Ancestry Equity & 3-Layer BA1 Regulatory QA Sidecar** for South Asian (`SAS`) and African (`AFR`) genomic cohorts.

Created & maintained by **Nikhil Bharti** ([@Nikhilbharti284](https://github.com/Nikhilbharti284) | [ORCID 0009-0006-3130-9141](https://orcid.org/0009-0006-3130-9141)).

## Live-Verified Benchmarks (118 Genes, 100 Homology Clusters, n = 9,142)

| Cohort | n | Raw AM ECE | v1.0.0 ECE | ECE Cut | Mondrian Cov | Singletons | Empty Sets |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **Overall** | 9,142 | 0.0686 | **0.0122** | **-82.2%** | >= 0.90 | - | **0** |
| **EUR** | 2,209 | 0.0622 | 0.0343 | -44.9% | 0.9086 | 92.6% | 0 |
| **SAS** | 1,098 | 0.0527 | 0.0298 | -43.5% | 0.9144 | 97.9% | 0 |
| **AFR** | 1,412 | 0.0637 | 0.0233 | -63.4% | 0.9334 | 100.0% | 0 |
| **ULTRA_RARE** | 4,423 | 0.1192 | 0.0072 | -94.0% | 0.9767 | 100.0% | 0 |

- **118-Gene OOF AUROC:** `0.9839` (`95% CI [0.9818, 0.9858]`, `p < 0.001` vs AlphaMissense `0.9583`) | **Post-2023 Holdout (`n=7,420`):** `0.9846`
- **Wet-Lab MaveDB DMS Suite:** 4/5 wins (`PTEN`, `MSH2`, `LDLR`, `BRCA1`); `TP53` explicitly reported as **HONEST MISS**.
- **3-Layer ClinGen BA1 Safeguard (`v1.1.0`):** Protects **11 registered ClinGen BA1 exceptions** (`HBB_E7V`, `G6PD_V68M`, `HFE_C282Y`, `BTD_D444H`, `MEFV_E148Q`, `CFTR_R117H`, etc.) via Layer-1, enforces >=5% stand-alone BA1 Benign via Layer-2, and catches BS1-BA1 conflicts via Layer-3.

## Quickstart

    pip install -e . --no-deps --no-build-isolation
    gec-audit --vcf path/to/hospital.vcf --out-vcf artifacts/audited_output.vcf
    pytest -v tests/

## Disclaimer
Computational research and clinical-pipeline QA audit sidecar only — not a standalone diagnostic SaMD.
