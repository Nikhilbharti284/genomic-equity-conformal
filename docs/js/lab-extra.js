/* Genomic Equity Conformal (v1.0.0) — Extra Client-Side Tools (Tools 5–8) */
(function () {
  // Known v1.0.0 SHA-256 ledger for instant verification
  const KNOWN_HASHES = {
    "81aa2d1e82d2c37c65bf1af53e24e3585d06772f7c97ce4cf332d410997be506": "v10_119genes_final_report.json (0.9839 AUROC, 118 Genes)",
    "abb610c4827c312b5dc1242c82b29302f5cc585c841166e0fc3358477332c656": "equity_model_v10_119genes_ccby4.joblib (Commercial-Safe Bundle)",
    "62d4b76e9344b6f99e96ef8c93e6b82efb5d78301eb18893b24e5c389afed6e0": "v1_mavedb_5gene_report.json (4/5 MaveDB Wins)",
    "ccd6d4a040dbf7255ab40eb0d7f307e10b0e2fb015812070159833acd41d2912": "v1_clinvar_120genes.csv (9,247 >=2-Star Variants)",
    "cf968fba2903f81ec5ac641b5f56ad2268a33131d8bb5d15450e70c4c20489e6": "equity_model_v10_clean.joblib (20-Gene Zero-REVEL Bundle)",
    "e32c23e2a93c1558434fc77cafc4f9d76b322ab6a140cd84f2c590cdc8a74021": "raw_vcf_audit_output.json (3/3 Raw VCF Verified)",
    "306b8ab5b80adcf225344237a7fffbff6ae54c50124695e655fddd7f2f2bde87": "hg38_5fold_v2_report.json (V2 Pure Zero-Circularity)",
    "8c1db848229b1be9c88b5dc93420be5caa307805fe7e9b42af6cac209d676029": "data/pilot_variants.csv (1,206 Pilot Variants)"
  };

  // Tool 5: 4-Ancestry Side-by-Side Equity Comparator
  const cmpBtn = document.getElementById("run-anc-cmp");
  const cmpBody = document.getElementById("anc-cmp-body");
  const COHORTS = [
    { id: "EUR", name: "EUR (European)", q: 0.44, ece: "-44.9%", cov: "0.9086" },
    { id: "SAS", name: "SAS (South Asian)", q: 0.41, ece: "-43.5%", cov: "0.9144" },
    { id: "AFR", name: "AFR (African)", q: 0.36, ece: "-63.4%", cov: "0.9334" },
    { id: "ULTRA_RARE", name: "ULTRA_RARE", q: 0.29, ece: "-94.0%", cov: "0.9767" }
  ];
  function runAncestryCompare() {
    if (!cmpBody) return;
    const am = Number(document.getElementById("sim-am")?.value || 0.86);
    const llr = Number(document.getElementById("sim-llr")?.value || -2.28);
    const af = Number(document.getElementById("sim-af")?.value || 0.047);
    const gene = (document.getElementById("sim-gene")?.value || "HBB").toUpperCase();
    const vari = (document.getElementById("sim-var")?.value || "E7V").toUpperCase();
    const vid = `${gene}_${vari}`;
    const logAf = Math.log10(Math.max(0, af) + 1e-7);
    const pJoint = 1 / (1 + Math.exp(-(4.6 * (am - 0.5) - 0.38 * llr - 0.82 * (logAf + 4.2))));

    cmpBody.innerHTML = COHORTS.map(c => {
      const s = [];
      if (pJoint <= c.q) s.push("Benign");
      if ((1 - pJoint) <= c.q) s.push("Pathogenic");
      if (!s.length) s.push(pJoint >= 0.5 ? "Pathogenic" : "Benign");
      let tier = s.length === 1 ? (s[0] === "Pathogenic" ? "Pathogenic_90pct" : "Benign_90pct") : "VUS_Review_Required";
      if (["HBB_E7V","HBB_E7K","GJB2_V37I","G6PD_V68M","CFTR_R117H","MEFV_M694V"].includes(vid) && tier === "Benign_90pct") {
        tier = "VUS (Layer-1 BA1 Protected)";
      }
      const cls = tier.startsWith("Pathogenic") ? "chip-danger" : (tier.startsWith("Benign") ? "chip-pass" : "chip-warn");
      return `<tr><td><strong>${c.name}</strong></td><td class="num">${c.q.toFixed(2)}</td><td class="num">${c.ece}</td><td class="num">${c.cov}</td><td><code>[${s.join(", ")}]</code></td><td><span class="chip ${cls}">${tier}</span></td></tr>`;
    }).join("");
  }
  if (cmpBtn) { cmpBtn.addEventListener("click", runAncestryCompare); runAncestryCompare(); }

  // Tool 6: Local .vcf File Upload (Zero-Egress) + Print QA Report
  const vcfFile = document.getElementById("vcf-file-upload");
  if (vcfFile) {
    vcfFile.addEventListener("change", () => {
      const f = vcfFile.files && vcfFile.files[0];
      if (!f) return;
      const reader = new FileReader();
      reader.onload = () => {
        const box = document.getElementById("vcf-input");
        if (box) box.value = String(reader.result || "");
        document.getElementById("vcf-run")?.click();
      };
      reader.readAsText(f);
    });
  }
  const printBtn = document.getElementById("print-qa-btn");
  if (printBtn) printBtn.addEventListener("click", () => window.print());

  // Tool 7: Live Web-Crypto SHA-256 Provenance Verifier
  const hashFile = document.getElementById("sha-file-input");
  const hashOut = document.getElementById("sha-verify-out");
  if (hashFile && hashOut) {
    hashFile.addEventListener("change", async () => {
      const f = hashFile.files && hashFile.files[0];
      if (!f) return;
      hashOut.innerHTML = "Computing local Web-Crypto SHA-256 (zero upload)...";
      const buf = await f.arrayBuffer();
      const digest = await crypto.subtle.digest("SHA-256", buf);
      const hex = Array.from(new Uint8Array(digest)).map(b => b.toString(16).padStart(2, "0")).join("");
      const match = KNOWN_HASHES[hex];
      hashOut.innerHTML = match
        ? `<span class="chip chip-pass">VERIFIED RO-CRATE MATCH</span> <strong>${match}</strong><br><code class="hash-cell">${hex}</code>`
        : `<span class="chip chip-warn">COMPUTED LOCALLY (Not in v1.0.0 core 8-file ledger)</span><br>File: <strong>${f.name}</strong> (${(f.size/1024).toFixed(1)} KB)<br><code class="hash-cell">${hex}</code>`;
    });
  }

  // Tool 8: Live 5-mer Homology Paralog Leakage Checker (v1_cluster_genes.py)
  function kmerSet(seq, k = 5) {
    const s = seq.toUpperCase().replace(/[^A-Z]/g, "");
    const out = new Set();
    for (let i = 0; i <= s.length - k; i++) out.add(s.slice(i, i + k));
    return out;
  }
  function runHomology() {
    const s1 = document.getElementById("hom-s1")?.value || "";
    const s2 = document.getElementById("hom-s2")?.value || "";
    const out = document.getElementById("hom-out");
    if (!out) return;
    const k1 = kmerSet(s1, 5), k2 = kmerSet(s2, 5);
    if (!k1.size || !k2.size) { out.innerHTML = "Enter two protein sequences (&ge;5 amino acids each)."; return; }
    const shared = [...k1].filter(x => k2.has(x));
    const overlap = shared.length / Math.min(k1.size, k2.size);
    const merged = overlap >= 0.08;
    out.innerHTML = `5-mer Containment = <strong>${(overlap * 100).toFixed(1)}%</strong> (${shared.length} shared 5-mers, threshold &ge; 8.0%) &rarr; ` +
      (merged
        ? `<span class="chip chip-danger">PARALOG DETECTED — MERGED INTO SAME CLUSTER</span> (Prevents Train/Test Leakage!)`
        : `<span class="chip chip-pass">INDEPENDENT GENES — SEPARATE CLUSTERS</span> (Safe to split across folds)`);
  }
  const homBtn = document.getElementById("run-hom-btn");
  if (homBtn) { homBtn.addEventListener("click", runHomology); runHomology(); }
  document.getElementById("hom-paralog")?.addEventListener("click", () => {
    document.getElementById("hom-s1").value = "MAAQVAPAAASSLGNPPPPPPSELKKAEQQQREEAGGEAAAAAAAERGEMKAAAGQESEGPAVGPPQPLGKELQDGAESNGGGGGGGAGSGGGPGAEPDLKNSNGNAGPRPALNNNLTEPPGGGGGGSSDGVGAPPHSAAAALPPPAYGFGQPYGRS";
    document.getElementById("hom-s2").value = "MAAQVAPAAASSLGNPPPPPPSELKKAEQQQREEAGGEAAAVAAAERGEMKAAAGQESEGPAVGPPQPLGKELQDGAESNGGGGGGAGSGGGPGAEPDLKNSNGNAGPRPALNNNLTEPPGGGGGGSSDGVGAPPHSAAAALPPPAYGFGQPYGRS";
    runHomology();
  });
  document.getElementById("hom-indep")?.addEventListener("click", () => {
    document.getElementById("hom-s1").value = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRV";
    document.getElementById("hom-s2").value = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH";
    runHomology();
  });
})();
