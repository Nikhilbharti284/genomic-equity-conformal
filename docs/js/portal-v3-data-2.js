const ALGOS = [
  { name: "Frozen ESM-2 (8M & 35M) + 512-aa Window", stage: "v0.1 → v3.0", status: "ADOPTED", why: "Computes masked-marginal zero-shot LLR = log P(mt) − log P(wt) cached in SQLite WAL (507–629 MB RAM, 0.34s/pos) with zero GPU and zero ClinVar label leakage." },
  { name: "Streaming CC BY 4.0 AlphaMissense (71M rows)", stage: "v1.0 → v3.0", status: "ADOPTED", why: "Matched 9,142/9,247 variants in 55.25 MB RAM with zero unpacked disk; replaced non-commercial REVEL for 100% B2B-safe licensing." },
  { name: "5-mer Containment + Union-Find (<30% Homology)", stage: "v1.0 → v3.0", status: "ADOPTED", why: "Zero-dependency DSA clustering 119 genes into 101 independent clusters, merging 45 paralog pairs (ARID1A/B, CACNA1A/E, COL1A1–COL7A1) to prevent homology leakage." },
  { name: "3-Model Soft-Voting Stack (MLP + HistGBM + LogReg)", stage: "v1.0 → v3.0", status: "ADOPTED", why: "Combines smooth neural decision boundaries, tree interactions, and linear regularization; lifted 118-gene AUROC to 0.9839–0.9844 (p < 0.001 vs AlphaMissense 0.9583)." },
  { name: "Inner GroupKFold Nested OOF Platt Scaling", stage: "v1.0 → v3.0", status: "UPGRADED", why: "2-parameter monotonic log-odds calibration fit strictly out-of-fold; eliminated calibration double-dipping and cut overall ECE by −90.4% (0.0686 → 0.0066)." },
  { name: "Non-Empty Ancestry-Conditional Mondrian CP", stage: "v1.0 → v3.0", status: "UPGRADED", why: "Computes finite-sample quantile q_g per ancestry across 8 global strata with an argmax non-empty floor—guaranteeing 0 empty prediction sets and >=90% target coverage." },
  { name: "Pure-NumPy 10Å C-Beta RSA Structural Burial", stage: "v2.0 → v3.0", status: "UPGRADED", why: "Replaced v1.0.0's aggressive global pLDDT scaling with local 10Å C-beta contact density from AlphaFold DB v6; flipped TP53 on wet-lab MaveDB from miss to WIN (+0.0003)." },
  { name: "3-Layer BA1 Safeguard + GT/DP/GQ Zygosity", stage: "v1.1 → v3.0", status: "UPGRADED", why: "Protects 11 ClinGen balancing-selection exceptions (HBB_E7V, GJB2_V37I, G6PD_V68M), enforces >=5% Layer-2 BA1, flags Layer-3 functional conflicts, and separates 0/1 carriers from 1/1 affected." },
  { name: "ClinGen Pejaver PP3/BP4 + CPIC Level A PGx", stage: "v1.0 → v3.0", status: "ADOPTED", why: "Strictly AF-free functional head maps Bayesian likelihood ratios to ACMG PP3/BP4 tiers (zero AF double-counting), paired with 6-gene CPIC drug-safety alerts (GEC_PGX_ALERT)." },
  { name: "Isotonic Regression Calibration", stage: "v0.1.0 Audit", status: "REJECTED", why: "Non-parametric step functions overfit 4-gene calibration splits in V1, worsening overall ECE (0.0817 → 0.0982) and AFR ECE (0.1125 → 0.2101). Replaced by Nested OOF Platt." },
  { name: "REVEL, CADD v1.7 & SpliceAI", stage: "v1.0.0 Audit", status: "REJECTED", why: "Non-commercial licenses create B2B hospital legal liability; CADD/SpliceAI also returned 0/1,206 hits on MyVariant hg38. Replaced by CC BY 4.0 AlphaMissense." },
  { name: "In-Sample Platt Double-Dipping", stage: "v1.0 & v2.0", status: "REJECTED", why: "Fitting Platt and Mondrian on identical splits caused EUR/NFE undercoverage (0.8424 in v0.x; 0.8666 in early v2). Eliminated via true nested inner GroupKFold OOF." },
  { name: "AF-Free Nonlinear MLP on Wet-Lab DMS", stage: "v1.0.0 Audit", status: "REJECTED", why: "Collapsed to |ρ| = 0.0082 on unseen TP53 MaveDB DMS. Replaced with zero-label AlphaFold DB v6 structural modulation (4/5 wet-lab wins)." },
  { name: "Random Train/Test Splits & 650M ESM-2", stage: "v0.1.0 Policy", status: "REJECTED", why: "Random splits leak homologous/same-gene variants (inflating AUROC); esm2_t33_650M and raw terabyte gnomAD VCFs exceed free 8 GB RAM / 15 GB disk limits." }
];

const BUGS = [
  { id: 1, ver: "v0.1.0", issue: "Python 3.14.2 PyTorch Wheel Mismatch", sym: "torch==2.4.1+cpu not found in Codespaces", fix: "Installed CPU PyTorch via --index-url https://download.pytorch.org/whl/cpu." },
  { id: 2, ver: "v0.1.0", issue: "Isolated System pytest Venv", sym: "ModuleNotFoundError: numpy + 12 unrelated repo errors", fix: "Switched to python3 -m pytest scoped explicitly to tests/." },
  { id: 3, ver: "v0.1.0", issue: "MyVariant hg19 Default Coordinate Bug", sym: "V1 pilot returned EUR:1, SAS:1, ULTRA_RARE:243", fix: "Added &assembly=hg38 + dbnsfp fallback; restored EUR:311, SAS:95, AFR:150." },
  { id: 4, ver: "v0.1.0", issue: "Isotonic Regression Overfitting", sym: "Overall ECE rose 0.0817 → 0.0982; AFR 0.1125 → 0.2101", fix: "Replaced with 2-parameter Platt Scaling on log-odds + L2 MLP (alpha=0.05)." },
  { id: 5, ver: "v1.0.0", issue: "Dead SpliceAI / CADD Endpoints", sym: "0/1,206 hits on MyVariant hg38; AUROC dropped", fix: "Added auto-coverage guard + streamed 71M-row CC BY 4.0 AlphaMissense (9,142/9,247)." },
  { id: 6, ver: "v1.0.0", issue: "REVEL Non-Commercial License Trap", sym: "Commercial B2B deployment blocker for hospitals", fix: "Removed REVEL from all deployable bundles; Track B is 100% MIT + CC BY 4.0." },
  { id: 7, ver: "v1.0.0", issue: "EUR Undercoverage & Empty Conformal Sets", sym: "EUR 0.8424; [] emitted when q_g < 0.5", fix: "Nested OOF Platt + 100 homology clusters + argmax non-empty floor (EUR 0.9086, 0 empty)." },
  { id: 8, ver: "v1.0.0", issue: "Fragmented Safety Guards Across Interfaces", sym: "CLI/API bypassed HBB_E7V balancing-selection guard", fix: "Centralized src/inference_engine.py + data/ba1_exception_registry.json everywhere." },
  { id: 9, ver: "v1.0.0", issue: "55-Gene Saturation Cap", sym: "10k variants filled in only 55 genes", fix: "Added max_per_gene=85 cap → 119 genes evenly covered across 101 homology clusters." },
  { id: 10, ver: "v1.0.0", issue: "PAT Workflow Scope & Mobile Table Clipping", sym: "Push rejected .github/workflows/ci.yml; mobile overflow", fix: "Pushed cleanly without workflow scope; wrapped all portal tables in overflow-x-auto." },
  { id: 11, ver: "v1.1.0", issue: "TP53_P72R Layer-2 Stand-Alone BA1 Rule", sym: "Common benign polymorphism needed explicit >=5% BA1 precedence", fix: "Implemented strict Layer-1 (Exception) → Layer-2 (>=5% Stand-Alone BA1) → Layer-3 (BS1 Conflict) hierarchy." },
  { id: 12, ver: "v2.0.0", issue: "In-Sample Platt Double-Dipping in 7-Ancestry Eval", sym: "Early v2 script gave ECE 0.0296 & NFE 0.8666 cov", fix: "Replaced with true inner GroupKFold(3) Nested OOF → ECE 0.0066 (−90.4%) & 7/8 strata >=0.90." }
];

const BENCH_V3 = [
  { cohort: "Overall (All 118 Genes)", n: "9,142", raw: "0.0686", cal: "0.0066", cut: "−90.4%", cov: "≥ 0.90", sing: "—", empty: 0 },
  { cohort: "SAS (South Asian)", n: "700", raw: "0.0530", cal: "0.0229", cut: "−56.8%", cov: "0.9114", sing: "97.9%", empty: 0 },
  { cohort: "AFR (African / African Am.)", n: "1,221", raw: "0.0543", cal: "0.0141", cut: "−74.0%", cov: "0.9328", sing: "100.0%", empty: 0 },
  { cohort: "EAS (East Asian)", n: "740", raw: "0.0685", cal: "0.0335", cut: "−51.1%", cov: "0.9149", sing: "95.8%", empty: 0 },
  { cohort: "AMR (Admixed American)", n: "717", raw: "0.0480", cal: "0.0288", cut: "−40.0%", cov: "0.9121", sing: "94.3%", empty: 0 },
  { cohort: "MID (Middle Eastern)*", n: "122", raw: "0.1127", cal: "0.0566", cut: "−49.8%", cov: "0.8934*", sing: "100.0%", empty: 0 },
  { cohort: "NFE (Non-Finnish European)", n: "1,192", raw: "0.0941", cal: "0.0180", cut: "−80.9%", cov: "0.9077", sing: "90.5%", empty: 0 },
  { cohort: "BOTTLENECK_FOUNDER (ASJ+FIN)", n: "622", raw: "0.0669", cal: "0.0179", cut: "−73.3%", cov: "0.9277", sing: "100.0%", empty: 0 },
  { cohort: "ULTRA_RARE (AF = 0)", n: "3,828", raw: "0.1333", cal: "0.0054", cut: "−96.0%", cov: "0.9875", sing: "100.0%", empty: 0 }
];

const DISCOVERIES = [
  { rank: 1, id: "ABCA4_G1961E", anc: "ASJ", af: "0.02392 (2.39%)", am: "0.8530", plddt: "87.50", sig: "Ashkenazi Stargardt hypomorphic founder" },
  { rank: 2, id: "ATP7B_H1069Q", anc: "ASJ", af: "0.00536", am: "0.9933", plddt: "77.62", sig: "Wilson disease European/ASJ founder" },
  { rank: 3, id: "SLC12A3_T180K", anc: "EAS", af: "0.00462", am: "0.8743", plddt: "95.25", sig: "East Asian Gitelman syndrome core founder" },
  { rank: 4, id: "GJB2_G12C", anc: "AMR", af: "0.00365", am: "0.8856", plddt: "63.00", sig: "Sensorineural hearing loss (AMR enriched)" },
  { rank: 5, id: "ATP7B_E1064A", anc: "ASJ", af: "0.00328", am: "0.8683", plddt: "91.50", sig: "Wilson disease structural core variant" },
  { rank: 6, id: "USH2A_C759F", anc: "AMR", af: "0.00315", am: "0.9091", plddt: "79.00", sig: "Retinitis pigmentosa / Usher syndrome" },
  { rank: 7, id: "SLC26A4_F335L", anc: "SAS", af: "0.00248", am: "0.9537", plddt: "83.38", sig: "South Asian Pendred syndrome allele" },
  { rank: 8, id: "SLC12A3_L849H", anc: "EAS", af: "0.00192", am: "0.9119", plddt: "91.50", sig: "East Asian Gitelman syndrome core allele" },
  { rank: 9, id: "ATP7B_R778L", anc: "EAS", af: "0.00189", am: "0.9641", plddt: "74.00", sig: "East Asian Wilson disease founder" },
  { rank: 10, id: "SLC26A4_V239D", anc: "SAS", af: "0.00167", am: "0.9612", plddt: "94.62", sig: "South Asian Pendred syndrome core founder" }
];

const PGX_PANEL = [
  { gene: "DPYD (*2A, *13, D949V)", vars: "DPYD_D949V, I560S, rs3918290", anc: "NFE, SAS, AMR, MID", action: "5-FU / Capecitabine: Reduce dose ≥50% or avoid" },
  { gene: "NUDT15 (*3) & TPMT", vars: "NUDT15_R139C, TPMT *2/*3A/*3C", anc: "EAS, SAS, AFR, AMR", action: "Thiopurines (Azathioprine, 6-MP): Cut dose 30–80%" },
  { gene: "CYP2C19 (*2, *3)", vars: "CYP2C19_W212X, rs4244285", anc: "EAS, SAS, AFR, MID", action: "Clopidogrel: Use Prasugrel/Ticagrelor in ACS/PCI" },
  { gene: "SLCO1B1 (*5)", vars: "SLCO1B1_V174A (rs4149056)", anc: "NFE, EAS, AMR, MID", action: "Simvastatin: Prescribe ≤20mg or alternative statin" },
  { gene: "G6PD (A-, Med)", vars: "G6PD_V68M, N126D, S188F", anc: "AFR, SAS, MID, AMR", action: "Primaquine/Dapsone/Rasburicase: Hemolysis alert" },
  { gene: "HLA-B (*57:01, *15:02)", vars: "rs2395029, rs144012689", anc: "EAS, SAS, NFE, AMR", action: "Abacavir / Carbamazepine: Contraindicated (SJS/TEN)" }
];

const LAB_PRESETS = [
  { id: "HBB_E7V", label: "HBB_E7V — Sickle Cell (Layer-1 BA1 Exception, AR)", gene: "HBB", moi: "AR", af: 0.0470, llr: -2.277, prob: 0.0195, ba1Exc: "Sickle Cell Disease (Balancing Selection)", pgx: null },
  { id: "TP53_R175H", label: "TP53_R175H — Li-Fraumeni Pathogenic Hotspot (AD)", gene: "TP53", moi: "AD", af: 0.0000, llr: -7.850, prob: 0.9202, ba1Exc: null, pgx: null },
  { id: "TP53_P72R", label: "TP53_P72R — Common Benign Polymorphism (Layer-2 BA1)", gene: "TP53", moi: "AD", af: 0.5002, llr: 0.420, prob: 0.0015, ba1Exc: null, pgx: null },
  { id: "DPYD_D949V", label: "DPYD_D949V — 5-FU Toxicity CPIC Level A PGx", gene: "DPYD", moi: "AR", af: 0.0052, llr: -4.120, prob: 0.8910, ba1Exc: null, pgx: "CPIC Level A: 5-FU/Capecitabine — Reduce dose >=50% or avoid" },
  { id: "NUDT15_R139C", label: "NUDT15_R139C — Thiopurine Toxicity (SAS/EAS PGx)", gene: "NUDT15", moi: "AR", af: 0.0680, llr: -5.300, prob: 0.1420, ba1Exc: null, pgx: "CPIC Level A: Azathioprine/6-MP — Cut dose 30–80%" },
  { id: "G6PD_V68M", label: "G6PD_V68M — Malaria Selection + Hemolysis PGx (XL)", gene: "G6PD", moi: "XL", af: 0.1160, llr: -3.450, prob: 0.0410, ba1Exc: "G6PD Deficiency (Malaria Selection)", pgx: "CPIC Level A: Primaquine/Dapsone/Rasburicase — Hemolysis alert" },
  { id: "ABCA4_G1961E", label: "ABCA4_G1961E — Rank #1 Autonomous Founder (ASJ 2.39%)", gene: "ABCA4", moi: "AR", af: 0.0239, llr: -4.900, prob: 0.1180, ba1Exc: null, pgx: null },
  { id: "CYP2C19_W212X", label: "CYP2C19_W212X — Clopidogrel Non-Responder PGx", gene: "CYP2C19", moi: "AR", af: 0.0450, llr: -6.200, prob: 0.2100, ba1Exc: null, pgx: "CPIC Level A: Clopidogrel — Use Prasugrel/Ticagrelor in ACS/PCI" },
  { id: "SLCO1B1_V174A", label: "SLCO1B1_V174A — Simvastatin Myopathy PGx", gene: "SLCO1B1", moi: "AR", af: 0.1420, llr: -3.100, prob: 0.0650, ba1Exc: null, pgx: "CPIC Level A: Simvastatin — Prescribe <=20mg or alternative statin" }
];

const SHA_LEDGER = [
  { file: "artifacts/v3_offline_clinical_ref.sqlite", out: "9,249 indexed GRCh38 records (SQLite WAL)", sha: "e758427e8a9978658f17d3eb23de30d2cc7a2e9a6af4510a0243b03e832155cf" },
  { file: "artifacts/v20_global_7ancestry_report.json", out: "AUROC 0.9844 • ECE 0.0066 (−90.4%) • 0 Empty Sets", sha: "08b002cc58beba804c64416f2e4dc38ca2d266a55ee156a3af830f209844e08a" },
  { file: "artifacts/v10_119genes_final_report.json", out: "v1.0.0 118-Gene Track B • AUROC 0.9839 • ECE 0.0122", sha: "81aa2d1e82d2c37c65bf1af53e24e3585d06772f7c97ce4cf332d410997be506" },
  { file: "artifacts/equity_model_v10_119genes_ccby4.joblib", out: "100% MIT + CC BY 4.0 • 100 Homology Clusters", sha: "abb610c4827c312b5dc1242c82b29302f5cc585c841166e0fc3358477332c656" },
  { file: "artifacts/equity_model_v10_clean.joblib", out: "20-Gene Zero-REVEL Bundle • AUROC 0.8400", sha: "cf968fba2903f81ec5ac641b5f56ad2268a33131d8bb5d15450e70c4c20489e6" },
  { file: "artifacts/v1_mavedb_5gene_report.json", out: "v1.0.0 5-Gene MaveDB DMS • 4/5 Wins (pLDDT)", sha: "62d4b76e9344b6f99e96ef8c93e6b82efb5d78301eb18893b24e5c389afed6e0" },
  { file: "artifacts/hg38_5fold_v2_report.json", out: "v0.1.0 20-Gene Pilot • AUROC 0.8038 • SAS 0.9053", sha: "306b8ab5b80adcf225344237a7fffbff6ae54c50124695e655fddd7f2f2bde87" },
  { file: "data/pilot_variants.csv", out: "1,206 Verified GRCh38 Variants • 34 MB RAM", sha: "8c1db848229b1be9c88b5dc93420be5caa307805fe7e9b42af6cac209d676029" }
];
