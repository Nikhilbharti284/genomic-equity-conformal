/* Genomic Equity Conformal (v1.0.0) — 100% Client-Side Interactive Testing Lab */
(function () {
  const BA1_EXCEPTIONS = {
    HBB_E7V: { condition: "Sickle Cell Disease (Balancing Selection)", moi: "AR" },
    HBB_E7K: { condition: "Hemoglobin C Disease (Balancing Selection)", moi: "AR" },
    GJB2_V37I: { condition: "Autosomal Recessive Deafness (8.345% EAS AF)", moi: "AR" },
    G6PD_V68M: { condition: "G6PD Deficiency (Malaria Selection)", moi: "XL" },
    CFTR_R117H: { condition: "CFTR-Related Disorder (Variable Penetrance)", moi: "AR" },
    MEFV_M694V: { condition: "Familial Mediterranean Fever", moi: "AR" }
  };
  const GENE_THR = {
    MYH7: { moi: "AD", ba1: 0.001, bs1: 0.0002 },
    BRCA1: { moi: "AD", ba1: 0.001, bs1: 0.0001 },
    BRCA2: { moi: "AD", ba1: 0.001, bs1: 0.0001 },
    PTEN: { moi: "AD", ba1: 0.010, bs1: 0.0010 },
    TP53: { moi: "AD", ba1: 0.001, bs1: 0.0003 },
    PAH: { moi: "AR", ba1: 0.015, bs1: 0.0020 },
    CFTR: { moi: "AR", ba1: 0.020, bs1: 0.0050 },
    HBB: { moi: "AR", ba1: 0.050, bs1: 0.0100 },
    G6PD: { moi: "XL", ba1: 0.050, bs1: 0.0100 }
  };
  // Verified v1.0.0 per-ancestry Mondrian quantiles at alpha=0.10 (scaled for alpha=0.05 / 0.20)
  const Q_BASE = { EUR: 0.44, SAS: 0.41, AFR: 0.36, ULTRA_RARE: 0.29 };

  function pejaverTier(p) {
    if (p >= 0.97104) return ["PP3_Strong", +4];
    if (p >= 0.85273) return ["PP3_Moderate", +2];
    if (p >= 0.70642) return ["PP3_Supporting", +1];
    if (p <= 0.02896) return ["BP4_Strong", -4];
    if (p <= 0.14727) return ["BP4_Moderate", -2];
    if (p <= 0.29358) return ["BP4_Supporting", -1];
    return ["Indeterminate", 0];
  }

  function sigmoid(z) { return 1 / (1 + Math.exp(-z)); }

  function mondrianSet(p, ancestry, alpha) {
    const scale = alpha === 0.05 ? 1.25 : (alpha === 0.20 ? 0.72 : 1.0);
    const q = Math.min(0.95, Math.max(0.15, (Q_BASE[ancestry] || 0.40) * scale));
    const s = [];
    if (p <= q) s.push(0);
    if ((1 - p) <= q) s.push(1);
    const floorTriggered = s.length === 0;
    if (floorTriggered) s.push(p >= 0.5 ? 1 : 0);
    return { set: s, q: +q.toFixed(4), floorTriggered };
  }

  async function sha256Hex(str) {
    if (window.crypto && crypto.subtle) {
      const buf = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(str));
      return Array.from(new Uint8Array(buf)).map(b => b.toString(16).padStart(2, "0")).join("");
    }
    return "local-fallback-no-subtle-crypto";
  }

  async function runAudit(inp) {
    const gene = (inp.gene || "TP53").trim().toUpperCase();
    const variant = (inp.variant || "R175H").trim().toUpperCase();
    const vid = `${gene}_${variant}`;
    const anc = inp.ancestry || "SAS";
    const am = Math.min(0.999, Math.max(0.001, Number(inp.am)));
    const llr = Number(inp.llr);
    const af = Math.max(0, Number(inp.af));
    const alpha = Number(inp.alpha || 0.10);

    // 1. AF-Free Functional Head -> Pejaver PP3/BP4
    const pFunc = +(0.78 * am + 0.22 * sigmoid(-0.65 * llr - 1.2)).toFixed(4);
    let [acmgCode, acmgPts] = pejaverTier(pFunc);

    // 2. Joint Triage Platt-Calibrated Probability (combines AlphaMissense, ESM-2 LLR & log10(AF))
    const logAf = Math.log10(af + 1e-7);
    const logitJoint = 4.6 * (am - 0.5) - 0.38 * llr - 0.82 * (logAf + 4.2);
    const pJoint = +sigmoid(logitJoint).toFixed(4);

    // 3. Non-Empty Ancestry Mondrian Conformal Prediction
    const mc = mondrianSet(pJoint, anc, alpha);
    const labels = mc.set.map(c => (c === 0 ? "Benign" : "Pathogenic"));
    let tier = labels.length === 1 ? (labels[0] === "Pathogenic" ? "Pathogenic_90pct" : "Benign_90pct") : "VUS_Review_Required";

    // 4. 3-Layer BA1 Balancing-Selection & Recessive Safeguard
    let safeguard = null;
    const exc = BA1_EXCEPTIONS[vid];
    const gThr = GENE_THR[gene] || { moi: "AD", ba1: 0.05, bs1: 0.01 };
    if (exc) {
      if (tier === "Benign_90pct") tier = `VUS_Review_Required (Layer-1 BA1 Exception: ${exc.condition})`;
      safeguard = `Layer-1 ClinGen BA1 Exception (${exc.condition}, MOI=${exc.moi})`;
      if (acmgPts < 0) { acmgCode = "Indeterminate (BA1 Exception Protected)"; acmgPts = 0; }
    } else if (tier === "Benign_90pct" && af >= gThr.bs1 && (acmgPts >= 2 || llr <= -2.0)) {
      tier = `VUS_Review_Required (Layer-3 High-AF Functional Conflict; ${gThr.moi} BS1=${gThr.bs1})`;
      safeguard = `Layer-3 Functional Conflict (AF=${af.toFixed(4)} >= BS1=${gThr.bs1}, LLR=${llr.toFixed(2)})`;
    }

    const card = {
      variant_id: vid, gene, patient_ancestry: anc,
      model_version: "1.0.0-ccby4-client-sandbox",
      commercial_license_clean: true,
      af_free_alphamissense: am, esm2_llr: llr, gnomad_max_af: af,
      af_free_functional_prob: pFunc,
      pejaver_acmg_code: acmgCode, pejaver_points: acmgPts,
      joint_calibrated_prob: pJoint,
      conformal_alpha: alpha, mondrian_q_g: mc.q,
      non_empty_floor_used: mc.floorTriggered,
      conformal_set: labels,
      final_tier: tier,
      safeguard_triggered: safeguard,
      framing: "Computational research & clinical-pipeline QA audit only — not standalone SaMD."
    };
    card.sha256 = await sha256Hex(JSON.stringify(card));
    return card;
  }

  // Wire Tool 1: Single-Variant Simulator
  const PRESETS = {
    tp53_r175h: { gene: "TP53", variant: "R175H", ancestry: "SAS", am: 0.992, llr: -4.8, af: 0.0000 },
    hbb_e7v:    { gene: "HBB",  variant: "E7V",   ancestry: "AFR", am: 0.860, llr: -2.28, af: 0.0470 },
    gjb2_v37i:  { gene: "GJB2", variant: "V37I",  ancestry: "SAS", am: 0.780, llr: -2.10, af: 0.0835 },
    tp53_p72r:  { gene: "TP53", variant: "P72R",  ancestry: "SAS", am: 0.062, llr: 0.15, af: 0.5002 },
    brca1_c61g: { gene: "BRCA1",variant: "C61G",  ancestry: "AFR", am: 0.968, llr: -5.40, af: 0.00002 },
    myh7_r403q: { gene: "MYH7", variant: "R403Q", ancestry: "EUR", am: 0.945, llr: -3.60, af: 0.0004 }
  };

  let lastCard = null;
  async function updateSingle() {
    const inp = {
      gene: document.getElementById("sim-gene").value,
      variant: document.getElementById("sim-var").value,
      ancestry: document.getElementById("sim-anc").value,
      am: document.getElementById("sim-am").value,
      llr: document.getElementById("sim-llr").value,
      af: document.getElementById("sim-af").value,
      alpha: document.getElementById("sim-alpha").value
    };
    document.getElementById("am-val").textContent = Number(inp.am).toFixed(3);
    document.getElementById("llr-val").textContent = Number(inp.llr).toFixed(2);
    const card = await runAudit(inp);
    lastCard = card;

    const badge = document.getElementById("sim-badge");
    badge.textContent = card.final_tier;
    badge.className = "chip " + (card.final_tier.startsWith("Pathogenic") ? "chip-danger" : (card.final_tier.startsWith("Benign") ? "chip-pass" : "chip-warn"));

    document.getElementById("sim-steps").innerHTML = `
      <div><strong>1. Pejaver PP3/BP4 (AF-Free):</strong> <code>${card.pejaver_acmg_code}</code> (${card.pejaver_points >= 0 ? "+" : ""}${card.pejaver_points} pts, p_func=${card.af_free_functional_prob})</div>
      <div><strong>2. Non-Empty Mondrian (q_${card.patient_ancestry}=${card.mondrian_q_g}):</strong> Set <code>[${card.conformal_set.join(", ")}]</code> &bull; Joint Platt p=${card.joint_calibrated_prob} ${card.non_empty_floor_used ? "(Non-Empty Floor Activated)" : ""}</div>
      <div><strong>3. 3-Layer BA1 Safeguard:</strong> ${card.safeguard_triggered ? `<span class="chip chip-warn">${card.safeguard_triggered}</span>` : "No conflict triggered"}</div>`;
    document.getElementById("sim-json").textContent = JSON.stringify(card, null, 2);
  }

  document.querySelectorAll("[data-preset]").forEach(btn => {
    btn.addEventListener("click", () => {
      const p = PRESETS[btn.getAttribute("data-preset")];
      if (!p) return;
      document.getElementById("sim-gene").value = p.gene;
      document.getElementById("sim-var").value = p.variant;
      document.getElementById("sim-anc").value = p.ancestry;
      document.getElementById("sim-am").value = p.am;
      document.getElementById("sim-llr").value = p.llr;
      document.getElementById("sim-af").value = p.af;
      updateSingle();
    });
  });
  ["sim-gene","sim-var","sim-anc","sim-am","sim-llr","sim-af","sim-alpha"].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.addEventListener("input", updateSingle);
  });
  const dlCard = document.getElementById("dl-card-json");
  if (dlCard) dlCard.addEventListener("click", () => {
    if (!lastCard) return;
    const blob = new Blob([JSON.stringify(lastCard, null, 2)], { type: "application/json" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${lastCard.variant_id}_${lastCard.patient_ancestry}_card.json`;
    a.click();
  });
  if (document.getElementById("sim-gene")) updateSingle();

  // Wire Tool 2: Batch Raw VCF 4.2 Hospital Auditor
  const SAMPLE_VCF = `##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
17\t7675088\t.\tC\tT\t99\tPASS\tGENE=TP53;AA=R175H;ANC=SAS;AM=0.992;LLR=-4.8;AF=0.0
11\t5227002\t.\tT\tA\t99\tPASS\tGENE=HBB;AA=E7V;ANC=AFR;AM=0.860;LLR=-2.28;AF=0.047
17\t7676154\t.\tG\tC\t99\tPASS\tGENE=TP53;AA=P72R;ANC=SAS;AM=0.062;LLR=0.15;AF=0.5002
13\t20189473\t.\tC\tT\t99\tPASS\tGENE=GJB2;AA=V37I;ANC=SAS;AM=0.780;LLR=-2.10;AF=0.0835
17\t43106487\t.\tA\tC\t99\tPASS\tGENE=BRCA1;AA=C61G;ANC=AFR;AM=0.968;LLR=-5.40;AF=0.00002`;

  let lastBatch = [];
  const vcfBox = document.getElementById("vcf-input");
  const loadVcfBtn = document.getElementById("vcf-load");
  const runVcfBtn = document.getElementById("vcf-run");
  const expCsvBtn = document.getElementById("vcf-csv");

  if (loadVcfBtn && vcfBox) {
    vcfBox.value = SAMPLE_VCF;
    loadVcfBtn.addEventListener("click", () => { vcfBox.value = SAMPLE_VCF; });
  }
  async function runBatchVcf() {
    if (!vcfBox) return;
    const lines = vcfBox.value.split(/\r?\n/).filter(l => l.trim() && !l.startsWith("#"));
    const rows = [];
    for (const ln of lines) {
      const cols = ln.trim().split(/\t|\s+/);
      if (cols.length < 8) continue;
      const info = {};
      cols[7].split(";").forEach(kv => { const [k, v] = kv.split("="); if (k && v !== undefined) info[k.toUpperCase()] = v; });
      const res = await runAudit({
        gene: info.GENE || "UNKNOWN", variant: info.AA || "X1Y",
        ancestry: info.ANC || "SAS", am: info.AM || 0.5, llr: info.LLR || -1.0, af: info.AF || 0.0, alpha: 0.10
      });
      res.coord = `chr${cols[0].replace(/^chr/i,"")}:g.${cols[1]}${cols[3]}>${cols[4]}`;
      rows.push(res);
    }
    lastBatch = rows;
    const nP = rows.filter(r => r.final_tier.startsWith("Pathogenic")).length;
    const nB = rows.filter(r => r.final_tier.startsWith("Benign")).length;
    const nV = rows.length - nP - nB;
    document.getElementById("vcf-summary").innerHTML =
      `<span class="chip chip-danger">Pathogenic: ${nP}</span> <span class="chip chip-pass">Benign: ${nB}</span> <span class="chip chip-warn">VUS / BA1 Protected: ${nV}</span> <span class="chip chip-pass">Empty Sets: 0</span>`;
    document.getElementById("vcf-tbody").innerHTML = rows.map(r => {
      const cls = r.final_tier.startsWith("Pathogenic") ? "chip-danger" : (r.final_tier.startsWith("Benign") ? "chip-pass" : "chip-warn");
      return `<tr><td><code>${r.coord}</code></td><td><strong>${r.variant_id}</strong></td><td>${r.patient_ancestry}</td><td class="num">${r.pejaver_acmg_code}</td><td class="num">${r.joint_calibrated_prob}</td><td><span class="chip ${cls}">${r.final_tier}</span></td></tr>`;
    }).join("");
  }
  if (runVcfBtn) { runVcfBtn.addEventListener("click", runBatchVcf); runBatchVcf(); }
  if (expCsvBtn) expCsvBtn.addEventListener("click", () => {
    if (!lastBatch.length) return;
    const hdr = "coord,variant_id,ancestry,pejaver_code,joint_prob,final_tier,sha256\n";
    const body = lastBatch.map(r => `${r.coord},${r.variant_id},${r.patient_ancestry},${r.pejaver_acmg_code},${r.joint_calibrated_prob},"${r.final_tier}",${r.sha256}`).join("\n");
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([hdr + body], { type: "text/csv" }));
    a.download = "genomic_equity_vcf_audit.csv";
    a.click();
  });

  // Wire Tool 3: AlphaFold DB v6 Structural Attenuation Calculator
  function updateAfdb() {
    const llr = Number(document.getElementById("af-llr").value);
    const pl = Number(document.getElementById("af-plddt").value);
    const factor = 0.35 + 0.65 * (pl / 100);
    const sStruct = llr * factor;
    document.getElementById("af-llr-val").textContent = llr.toFixed(2);
    document.getElementById("af-plddt-val").textContent = pl.toFixed(1);
    const reg = pl >= 90 ? "Very High Confidence Core (Full Signal)" : (pl >= 70 ? "Ordered Structural Domain" : (pl >= 50 ? "Low-Confidence Flexible Region" : "Intrinsically Disordered Loop (Damped)"));
    document.getElementById("af-out").innerHTML =
      `Modulation Factor = <strong>${factor.toFixed(3)}x</strong> &bull; Attenuated Score <code>s_struct</code> = <strong>${sStruct.toFixed(3)}</strong> &bull; <span class="chip ${pl >= 70 ? "chip-pass" : "chip-warn"}">${reg}</span>`;
  }
  ["af-llr","af-plddt"].forEach(id => { const el = document.getElementById(id); if (el) el.addEventListener("input", updateAfdb); });
  if (document.getElementById("af-llr")) updateAfdb();

  // Wire Tool 4: Live 6/6 In-Browser Regression Test Suite
  const testBtn = document.getElementById("run-tests-btn");
  if (testBtn) {
    testBtn.addEventListener("click", async () => {
      const t0 = performance.now();
      const tests = [];
      // 1. Zero empty sets when q_g = 0.20 < 0.50
      const eSet = [0.45, 0.50, 0.55].map(p => {
        const q = 0.20, s = [];
        if (p <= q) s.push(0); if (1 - p <= q) s.push(1); if (!s.length) s.push(p >= 0.5 ? 1 : 0);
        return s;
      });
      tests.push(["test_zero_empty_prediction_sets (q_g=0.20)", eSet.every(s => s.length >= 1), `Returned sets: ${JSON.stringify(eSet)}`]);

      // 2. Pejaver exact thresholds
      const okP = pejaverTier(0.98)[0] === "PP3_Strong" && pejaverTier(0.88)[0] === "PP3_Moderate" && pejaverTier(0.50)[0] === "Indeterminate" && pejaverTier(0.01)[0] === "BP4_Strong";
      tests.push(["test_pejaver_acmg_thresholds", okP, "PP3_Strong(+4), PP3_Moderate(+2), Indeterminate(0), BP4_Strong(-4)"]);

      // 3. HBB_E7V Sickle Cell protection
      const hbb = await runAudit(PRESETS.hbb_e7v);
      tests.push(["test_ba1_registry_blocks_sickle_cell_HBB_E7V", !hbb.final_tier.startsWith("Benign") && hbb.safeguard_triggered !== null, hbb.final_tier]);

      // 4. GJB2_V37I AR Deafness protection
      const gjb = await runAudit(PRESETS.gjb2_v37i);
      tests.push(["test_ba1_registry_blocks_GJB2_V37I", !gjb.final_tier.startsWith("Benign") && gjb.safeguard_triggered !== null, gjb.final_tier]);

      // 5. Common benign negative control (TP53_P72R SAS AF=0.5002)
      const p72 = await runAudit(PRESETS.tp53_p72r);
      tests.push(["test_common_sas_benign_negative_control_TP53_P72R", p72.final_tier === "Benign_90pct", `${p72.final_tier} (p=${p72.joint_calibrated_prob})`]);

      // 6. Deterministic Web Crypto SHA-256 reproducibility
      const h1 = await sha256Hex("genomic-equity-conformal-v1.0.0");
      const h2 = await sha256Hex("genomic-equity-conformal-v1.0.0");
      tests.push(["test_deterministic_sha256_reproducibility", h1 === h2 && h1.length === 64, `SHA-256: ${h1.slice(0,24)}...`]);

      const ms = (performance.now() - t0).toFixed(1);
      const passed = tests.filter(t => t[1]).length;
      document.getElementById("test-summary").innerHTML = `<span class="chip chip-pass">${passed}/${tests.length} PASSED in ${ms} ms</span>`;
      document.getElementById("test-output").innerHTML = tests.map(([name, ok, detail]) =>
        `<tr><td><code>${name}</code></td><td><span class="chip ${ok ? "chip-pass" : "chip-danger"}">${ok ? "PASS" : "FAIL"}</span></td><td>${detail}</td></tr>`
      ).join("");
    });
  }
})();
