const VERSIONS = {
  genesis: {
    badge: "GENESIS • NOV 2026 AUDIT",
    title: "8-Track Global-Health Blueprint & Fungal Pma1 Audit",
    why: "Started by auditing an early Fungal Pma1 molecular-glue proposal and catching 8 critical/high defects (8-core hardware mismatch, pooled ROC-AUC confounding across 18 PPI systems, fusicoccin 2O98 opposite activation direction, vacuolar sink omission, unverified patents). Scored 8 CPU-tractable tracks on N × F × Fu × I (max 625).",
    built: [
      "Track 7 (Genomic Equity) tied #1 overall at 500/625 with the highest Novelty (5/5), 5/5 Feasibility, 5/5 Fundability, and ~10 core-h/mo compute.",
      "Selected Track 7 to build first because it achieves $0 digital validation (ClinVar, gnomAD, MaveDB) without waiting for $1.5k–$2k wet-lab compound procurement.",
      "Adopted AlphaFold DB v6 and W3C RO-Crate SHA-256 provenance concepts from the OpenRepurpose-CPU hybrid blueprint."
    ],
    metrics: "Rubric Score: 500/625 • Monthly Budget: ~10 core-h • Wet-Lab Cost: $0"
  },
  v010: {
    badge: "v0.1.0 • 20-GENE PILOT",
    title: "Frozen ESM-2 (8M) + Platt + Mondrian Proof-of-Concept",
    why: "Proved that ancestry-conditional calibration works on free GitHub Codespaces (2-core CPU / 8 GB RAM) with strict zero-gene-leakage GroupKFold splits.",
    built: [
      "Stream-curated 1,206 verified ClinVar GRCh38 missense variants across 20 actionable genes (34 MB RAM, 90.06s).",
      "Cached Frozen facebook/esm2_t6_8M_UR50D zero-shot LLRs in SQLite WAL with SIGTERM timeout guards (507.06 MB peak RAM).",
      "Diagnosed & rejected Isotonic Regression overfitting (overall ECE worsened 0.0817 → 0.0982; AFR 0.1125 → 0.2101) and replaced it with 2-parameter Platt Scaling + L2 MLP (32,16)."
    ],
    metrics: "OOF AUROC: 0.8038 (+13.65% vs Raw ESM-2 0.6673) • PR-AUC: 0.8342 • SAS Cov: 0.9053 • AFR Cov: 0.8933"
  },
  v100: {
    badge: "v1.0.0 • COMMERCIAL-SAFE SCALE-UP",
    title: "118-Gene / 100-Cluster Track B (MIT + CC BY 4.0) & MaveDB Suite",
    why: "Upgraded because 20 genes had small calibration folds, REVEL carried non-commercial B2B licensing liability, and naive Mondrian could emit empty sets when q_g < 0.5.",
    built: [
      "Removed REVEL, CADD v1.7, and SpliceAI; streamed 71M-row CC BY 4.0 AlphaMissense in 55.25 MB RAM (matched 9,142/9,247 variants).",
      "Built 5-mer containment + Union-Find <30% homology clustering (119 genes → 101 clusters, merging 45 real paralogs like ARID1A/B, COL1A1–COL7A1).",
      "Introduced 3-Model Soft-Voting Stack (MLP + HistGBM + LogReg), Nested OOF Platt, Non-Empty Mondrian floor (0 empty sets), 5-gene wet-lab MaveDB DMS (4/5 wins; TP53 miss), and W3C RO-Crate 1.1 (42 SHA-256 hashes)."
    ],
    metrics: "OOF AUROC: 0.9839 (p < 0.001 vs AM 0.9583) • Post-2023 Holdout: 0.9846 (n=7,420) • ECE Cut: −82.2% (0.0122)"
  },
  v110: {
    badge: "v1.1.0 • CLINICAL PIPELINE HARDENING",
    title: "Unified GenomicEquityEngine, VCF 4.2 & 11 ClinGen BA1 Exceptions",
    why: "Upgraded to eliminate fragmented safety guards between CLI, VCF auditor, and REST API, and to fix Bug #11 (TP53_P72R stand-alone BA1 >=5% rule).",
    built: [
      "Centralized all inference, ClinGen Pejaver PP3/BP4 Bayesian tiers, and 3-Layer BA1 Safeguard inside src/inference_engine.py (zero REVEL).",
      "Expanded data/ba1_exception_registry.json to 11 ClinGen SVI BA1 exceptions + gene-specific MOI thresholds.",
      "Packaged pip-installable gec-audit & gec-predict CLI with hospital-ready annotated VCF 4.2 output (--out-vcf) and 12/12 passing pytest."
    ],
    metrics: "11 ClinGen BA1 Exceptions • Layer-1/2/3 Hierarchy • 12/12 Pytest Passing"
  },
  v200: {
    badge: "v2.0.0 • PAN-GLOBAL 7-ANCESTRY & PGX",
    title: "7-Ancestry Conformal, 19 Autonomous Discoveries, CPIC PGx & C-Beta RSA",
    why: "Upgraded from 4 cohorts to full Pan-Global 7-Ancestry coverage, added drug-safety pharmacogenomics, and fixed v1.0.0's TP53 wet-lab miss by replacing aggressive pLDDT scaling with 10Å C-Beta RSA burial.",
    built: [
      "Expanded to 8 strata (SAS, AFR, EAS, AMR, MID, NFE, BOTTLENECK_FOUNDER [ASJ+FIN], ULTRA_RARE) with inner GroupKFold(3) Nested OOF cutting overall ECE by −90.4% (0.0686 → 0.0066).",
      "Autonomously discovered 19 unlisted BS1+ pathogenic founder/balancing-selection variants absent from ClinGen's static BA1 list (top: ABCA4_G1961E, ATP7B_H1069Q, SLC12A3_T180K).",
      "Added 6-panel Global CPIC Level A PGx Sidecar (GEC_PGX_ALERT) and pure-NumPy 10Å C-Beta RSA Burial flipping TP53 from miss to WIN (+0.0003; 4/5 wins, LDLR honest miss)."
    ],
    metrics: "OOF AUROC: 0.9844 (PR-AUC 0.9898) • ECE: −90.4% (0.0066) • 19 Founders • TP53 MaveDB WIN"
  },
  v300: {
    badge: "v3.0.0 • ENTERPRISE AIR-GAPPED EDITION",
    title: "Offline SQLite WAL, GT Zygosity Gating, 10K WES Streaming & NABL/CAP QA",
    why: "Upgraded for zero-egress hospital production under DPDP Act / NABL ISO 15189 / CAP / EU AI Act: eliminated live API dependency, added sample-level genotype zygosity gating, and benchmarked full-exome streaming.",
    built: [
      "Indexed 9,249 GRCh38 records in artifacts/v3_offline_clinical_ref.sqlite (WAL mode) for sub-ms --network=none execution.",
      "Built src/vcf_patient_profiler.py parsing VCF GT/DP/GQ (separating 0/1 AR_Carrier_Only from 1/1 AR_Biallelic_Affected) + AIM 7-Ancestry MLE auto-detection.",
      "Verified 10,000-variant WES streaming in 0.917s at 181.1 MB peak RAM, stdlib OpenAPI REST API (gec-api), NABL/CAP/EU AI Act HTML report (gec-report), and 20/20 pytest."
    ],
    metrics: "9,249 Offline DB Rows • 10K WES: 0.917s / 181.1 MB RAM • 20/20 Pytest Passing"
  },
  v400: {
    badge: "v4.0.0 • GLOBAL HOSPITAL LAN EDITION",
    title: "Offline LAN Bridge (/ui), Trio Compound-Het, HL7 FHIR R4 & SHA-256 Hash-Chain",
    why: "Upgraded on top of v3.0.0 so hospitals worldwide (EU GDPR, UK NHS, US HIPAA, India DPDP, Brazil LGPD, Saudi PDPL) can connect multiple LAN computers and EHR/LIMS offline with zero cloud egress.",
    built: [
      "Built src/lan_hospital_bridge.py (0.0.0.0:8000, optional X-GEC-API-Key, zero-CDN /ui dashboard) with explicit OOD_NOT_IN_OFFLINE_PANEL abstention.",
      "Built src/trio_phasing.py detecting Trans Compound Heterozygotes (AR_Compound_Het_Trans_Affected), Cis Carriers, and De Novo candidates.",
      "Built src/fhir_exporter.py (HL7 FHIR R4 Bundle, LOINC 51969-4 / 69548-6) & src/audit_chain.py (tamper-evident SHA-256 ledger catching 1-byte edits); 27/27 pytest passing."
    ],
    metrics: "HL7 FHIR R4 • Trio Compound-Het • SHA-256 Hash-Chain • 27/27 Pytest (37.42s)"
  }
};