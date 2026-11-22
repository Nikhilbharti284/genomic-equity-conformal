function showToast(msg) {
  const t = document.getElementById('toast');
  t.textContent = msg;
  t.classList.remove('opacity-0', 'translate-y-2');
  t.classList.add('opacity-100', 'translate-y-0');
  clearTimeout(window._toastTimer);
  window._toastTimer = setTimeout(() => { t.classList.add('opacity-0', 'translate-y-2'); t.classList.remove('opacity-100', 'translate-y-0'); }, 2000);
}

function copyText(text, label) {
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(() => showToast(label + ' copied!')).catch(() => fallbackCopy(text, label));
  } else fallbackCopy(text, label);
}

function fallbackCopy(text, label) {
  const ta = document.createElement('textarea');
  ta.value = text; ta.style.position = 'fixed'; ta.style.opacity = '0';
  document.body.appendChild(ta); ta.select();
  try { document.execCommand('copy'); showToast(label + ' copied!'); } catch (e) { showToast('Copy unavailable'); }
  document.body.removeChild(ta);
}

function simpleHash(str) {
  let h1 = 0xdeadbeef ^ str.length, h2 = 0x41c6ce57 ^ str.length;
  for (let i = 0; i < str.length; i++) {
    const ch = str.charCodeAt(i);
    h1 = Math.imul(h1 ^ ch, 2654435761);
    h2 = Math.imul(h2 ^ ch, 1597334677);
  }
  h1 = Math.imul(h1 ^ (h1 >>> 16), 2246822507) ^ Math.imul(h2 ^ (h2 >>> 13), 3266489909);
  h2 = Math.imul(h2 ^ (h2 >>> 16), 2246822507) ^ Math.imul(h1 ^ (h1 >>> 13), 3266489909);
  return (4294967296 * (2097151 & h2) + (h1 >>> 0)).toString(16).padStart(16, '0');
}

let chartsInitialized = {};

function switchTab(tabId) {
  document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
  const panel = document.getElementById('tab-' + tabId);
  if (panel) panel.classList.add('active');
  document.querySelectorAll('.nav-btn').forEach(b => b.classList.toggle('active', b.dataset.tab === tabId));
  document.getElementById('mobileNav').classList.add('hidden');
  document.getElementById('menuIconOpen').classList.remove('hidden');
  document.getElementById('menuIconClose').classList.add('hidden');
  document.getElementById('mobileMenuBtn').setAttribute('aria-expanded', 'false');
  window.scrollTo({ top: 0, behavior: 'smooth' });
  if (tabId === 'evolution') { initEvolutionCharts(); replayActiveTabCharts('evolution'); }
  if (tabId === 'benchmarks') { initBenchmarkCharts(); replayActiveTabCharts('benchmarks'); }
}

document.querySelectorAll('.nav-btn').forEach(btn => btn.addEventListener('click', () => switchTab(btn.dataset.tab)));
document.querySelectorAll('.cta-nav').forEach(btn => btn.addEventListener('click', () => switchTab(btn.dataset.target)));
document.getElementById('brandHomeBtn').addEventListener('click', () => switchTab('overview'));
document.getElementById('mobileMenuBtn').addEventListener('click', () => {
  const m = document.getElementById('mobileNav');
  const open = m.classList.toggle('hidden');
  document.getElementById('menuIconOpen').classList.toggle('hidden', !open);
  document.getElementById('menuIconClose').classList.toggle('hidden', open);
  document.getElementById('mobileMenuBtn').setAttribute('aria-expanded', String(!open));
});

function renderVersion(verKey) {
  const v = VERSIONS[verKey];
  document.querySelectorAll('#versionBtns .pill-btn').forEach(b => b.classList.toggle('active', b.dataset.ver === verKey));
  document.getElementById('versionDetailCard').innerHTML = `
    <div class="flex flex-wrap items-center justify-between gap-2 mb-2">
      <span class="px-2.5 py-0.5 rounded-full text-[11px] font-bold bg-[#0f4c5c] text-white">${v.badge}</span>
      <span class="text-xs font-mono font-semibold text-emerald-700 bg-emerald-50 border border-emerald-200 px-2.5 py-0.5 rounded-lg">${v.metrics}</span>
    </div>
    <h3 class="text-lg font-extrabold text-slate-900">${v.title}</h3>
    <p class="text-xs sm:text-sm text-slate-700 mt-1.5"><strong>Why This Upgrade Happened:</strong> ${v.why}</p>
    <ul class="mt-3 space-y-1.5 list-disc pl-5">${v.built.map(i => `<li class="text-xs sm:text-sm text-slate-700">${i}</li>`).join('')}</ul>`;
}
document.querySelectorAll('#versionBtns .pill-btn').forEach(btn => btn.addEventListener('click', () => renderVersion(btn.dataset.ver)));

function renderAlgos(filter) {
  document.querySelectorAll('#algoFilters .pill-btn').forEach(b => b.classList.toggle('active', b.dataset.filter === filter));
  const badgeMap = { ADOPTED: 'bg-emerald-100 text-emerald-800 border-emerald-300', UPGRADED: 'bg-sky-100 text-sky-800 border-sky-300', REJECTED: 'bg-rose-100 text-rose-800 border-rose-300' };
  document.getElementById('algoTableBody').innerHTML = ALGOS.filter(a => filter === 'ALL' || a.status === filter).map((a, idx) => `
    <tr class="${idx % 2 === 1 ? 'bg-slate-50' : 'bg-white'}">
      <td class="p-3 font-semibold text-slate-900">${a.name}</td>
      <td class="p-3 font-mono text-xs text-slate-600 whitespace-nowrap">${a.stage}</td>
      <td class="p-3"><span class="px-2 py-0.5 rounded-full text-[11px] font-bold border ${badgeMap[a.status]}">${a.status}</span></td>
      <td class="p-3 text-slate-700">${a.why}</td>
    </tr>`).join('');
}
document.querySelectorAll('#algoFilters .pill-btn').forEach(btn => btn.addEventListener('click', () => renderAlgos(btn.dataset.filter)));

function renderBugs() {
  document.getElementById('bugTableBody').innerHTML = BUGS.map((b, idx) => `
    <tr class="${idx % 2 === 1 ? 'bg-slate-50' : 'bg-white'}">
      <td class="p-3 font-mono font-bold text-slate-500">#${b.id}</td>
      <td class="p-3 font-mono text-xs font-semibold text-[#0f4c5c] whitespace-nowrap">${b.ver}</td>
      <td class="p-3 font-semibold text-slate-900">${b.issue}</td>
      <td class="p-3 text-rose-700">${b.sym}</td>
      <td class="p-3 text-slate-700">${b.fix}</td>
    </tr>`).join('');
}

function renderStaticTables() {
  document.getElementById('benchTableBody').innerHTML = BENCH_V3.map((r, i) => `
    <tr class="${i === 0 ? 'bg-teal-50/70 font-semibold' : (i % 2 === 1 ? 'bg-slate-50' : 'bg-white')}">
      <td class="p-3 text-slate-900">${r.cohort}</td><td class="p-3 font-mono">${r.n}</td><td class="p-3 font-mono text-rose-700">${r.raw}</td>
      <td class="p-3 font-mono font-bold text-emerald-700">${r.cal}</td><td class="p-3 font-mono font-bold text-emerald-700">${r.cut}</td>
      <td class="p-3 font-mono">${r.cov}</td><td class="p-3 font-mono">${r.sing}</td><td class="p-3 font-mono text-center">${r.empty}</td>
    </tr>`).join('');
  document.getElementById('pgxTableBody').innerHTML = PGX_PANEL.map((p, i) => `
    <tr class="${i % 2 === 1 ? 'bg-slate-50' : 'bg-white'}">
      <td class="p-3 font-bold text-[#0f4c5c]">${p.gene}</td><td class="p-3 font-mono text-xs text-slate-700">${p.vars}</td>
      <td class="p-3 font-semibold text-amber-700">${p.anc}</td><td class="p-3 text-slate-800">${p.action}</td>
    </tr>`).join('');
  document.getElementById('shaTableBody').innerHTML = SHA_LEDGER.map((s, i) => `
    <tr class="${i % 2 === 1 ? 'bg-slate-50' : 'bg-white'}">
      <td class="p-3 font-mono font-semibold text-[#0f4c5c]">${s.file}</td><td class="p-3 text-slate-700">${s.out}</td>
      <td class="p-3 font-mono text-[11px] text-slate-600 break-all">${s.sha}</td>
    </tr>`).join('');
}

function renderDiscoveries() {
  const q = document.getElementById('discSearch').value.trim().toLowerCase();
  const anc = document.getElementById('discAncFilter').value;
  const filtered = DISCOVERIES.filter(d => (anc === 'ALL' || d.anc === anc) && (!q || d.id.toLowerCase().includes(q) || d.sig.toLowerCase().includes(q)));
  document.getElementById('discTableBody').innerHTML = filtered.length ? filtered.map((d, i) => `
    <tr class="${i % 2 === 1 ? 'bg-slate-50' : 'bg-white'}">
      <td class="p-3 font-mono font-bold text-slate-500">#${d.rank}</td><td class="p-3 font-mono font-bold text-[#0f4c5c]">${d.id}</td>
      <td class="p-3"><span class="px-2 py-0.5 rounded-full text-[11px] font-bold bg-amber-100 text-amber-800">${d.anc}</span></td>
      <td class="p-3 font-mono">${d.af}</td><td class="p-3 font-mono font-semibold text-emerald-700">${d.am}</td>
      <td class="p-3 font-mono">${d.plddt}</td><td class="p-3 text-slate-700">${d.sig}</td>
    </tr>`).join('') : `<tr><td colspan="7" class="p-4 text-center text-slate-500">No matching founder variants found.</td></tr>`;
}
document.getElementById('discSearch').addEventListener('input', renderDiscoveries);
document.getElementById('discAncFilter').addEventListener('change', renderDiscoveries);

function initLab() {
  document.getElementById('labPreset').innerHTML = LAB_PRESETS.map(p => `<option value="${p.id}">${p.label}</option>`).join('');
  runAudit(); runSimulator();
}

function runAudit() {
  const id = document.getElementById('labPreset').value, gt = document.getElementById('labGt').value, anc = document.getElementById('labAnc').value;
  const p = LAB_PRESETS.find(x => x.id === id) || LAB_PRESETS[0];
  let tier = p.prob >= 0.85 ? "Pathogenic_90pct" : (p.prob <= 0.15 ? "Benign_90pct" : "VUS_Review_Required"), safeguard = null;
  if (p.ba1Exc) { if (tier === "Benign_90pct") tier = `VUS_Review_Required (Layer-1 BA1 Exception: ${p.ba1Exc})`; safeguard = `Layer-1 BA1 Exception (${p.ba1Exc})`; }
  else if (p.af >= 0.05 && p.llr > -2.0) { tier = "Benign_90pct (Layer-2 Stand-Alone BA1 >= 5%)"; safeguard = "Layer-2 Stand-Alone BA1"; }
  else if (tier === "Benign_90pct" && p.af >= 0.005 && p.llr <= -2.0) { tier = "VUS_Review_Required (Layer-3 High-AF Functional Conflict)"; safeguard = `Layer-3 Functional Conflict (AF=${p.af}, LLR=${p.llr})`; }
  let zygosity = "Wildtype_Reference";
  if (gt === "0/1") zygosity = p.moi === "AR" ? "AR_Carrier_Only (Heterozygous)" : (p.moi === "XL" ? "XL_Heterozygous_Carrier_or_Risk" : "AD_Heterozygous_Affected_Risk");
  else if (gt === "1/1") zygosity = p.moi === "AR" ? "AR_Biallelic_Affected" : (p.moi === "XL" ? "XL_Hemizygous_or_Homozygous_Affected" : "AD_Homozygous_Affected");
  const payload = { variant_id: p.id, gene: p.gene, moi: p.moi, patient_ancestry: anc, sample_gt: gt, clinical_zygosity_status: zygosity, gnomad_grpmax_af: p.af, esm2_llr: p.llr, platt_calibrated_prob: p.prob, final_tier: tier, safeguard_triggered: safeguard, gec_pgx_alert: p.pgx || "None", empty_prediction_sets: 0, offline_db_mode: true, framing: "Computational research & clinical-pipeline QA audit only — not standalone SaMD." };
  payload.audit_digest = simpleHash(JSON.stringify(payload));
  document.getElementById('auditJsonOut').textContent = JSON.stringify(payload, null, 2);
  document.getElementById('auditSummaryBadges').innerHTML = `<span class="px-2.5 py-1 rounded-lg text-xs font-bold bg-teal-50 text-[#0f4c5c] border border-teal-200">Zygosity: ${zygosity}</span><span class="px-2.5 py-1 rounded-lg text-xs font-bold ${safeguard ? 'bg-amber-100 text-amber-900 border-amber-300' : 'bg-slate-100 text-slate-700 border-slate-200'} border">Safeguard: ${safeguard ? 'TRIGGERED' : 'None'}</span>${p.pgx ? '<span class="px-2.5 py-1 rounded-lg text-xs font-bold bg-rose-100 text-rose-800 border border-rose-300">CPIC PGx Alert Active</span>' : ''}`;
}

function runSimulator() {
  const prob = parseFloat(document.getElementById('simProb').value) || 0.5, llr = parseFloat(document.getElementById('simLlr').value) || 0, af = Math.max(0, parseFloat(document.getElementById('simAf').value) || 0), moi = document.getElementById('simMoi').value;
  document.getElementById('simProbVal').textContent = prob.toFixed(2); document.getElementById('simLlrVal').textContent = llr.toFixed(2);
  const bs1 = { AD: 0.0003, AR: 0.0050, XL: 0.0100 }[moi];
  let pejaver = ["Indeterminate", 0];
  if (prob >= 0.97104) pejaver = ["PP3_Strong", "+4"]; else if (prob >= 0.85273) pejaver = ["PP3_Moderate", "+2"]; else if (prob >= 0.70642) pejaver = ["PP3_Supporting", "+1"]; else if (prob <= 0.02896) pejaver = ["BP4_Strong", "-4"]; else if (prob <= 0.14727) pejaver = ["BP4_Moderate", "-2"]; else if (prob <= 0.29358) pejaver = ["BP4_Supporting", "-1"];
  const q = 0.18; let set = [];
  if (prob <= q) set.push("Benign"); if ((1 - prob) <= q) set.push("Pathogenic"); if (!set.length) set.push(prob >= 0.5 ? "Pathogenic" : "Benign"); if (prob > 0.20 && prob < 0.80) set = ["Benign", "Pathogenic"];
  let tier = set.length === 1 ? (set[0] === "Pathogenic" ? "Pathogenic_90pct" : "Benign_90pct") : "VUS_Review_Required", safeguard = "None";
  if (af >= 0.05 && llr > -2.0) { tier = "Benign_90pct (Layer-2 Stand-Alone BA1)"; safeguard = "Layer-2 BA1 (AF >= 0.05)"; }
  else if (tier === "Benign_90pct" && af >= bs1 && llr <= -2.0) { tier = `VUS_Review_Required (Layer-3 High-AF Functional Conflict; ${moi} BS1=${bs1})`; safeguard = `Layer-3 Conflict (AF=${af} >= BS1=${bs1}, LLR=${llr})`; }
  document.getElementById('simJsonOut').textContent = JSON.stringify({ input: { platt_prob: prob, grpmax_af: af, esm2_llr: llr, moi: moi, bs1_threshold: bs1 }, pejaver_acmg_tier: pejaver[0], pejaver_points: pejaver[1], conformal_90pct_set: set, final_tier: tier, safeguard: safeguard, empty_prediction_sets: 0 }, null, 2);
  document.getElementById('simBadges').innerHTML = `<span class="px-2.5 py-1 rounded-lg text-xs font-bold bg-sky-50 text-sky-900 border border-sky-200">Pejaver: ${pejaver[0]} (${pejaver[1]})</span><span class="px-2.5 py-1 rounded-lg text-xs font-bold bg-emerald-50 text-emerald-900 border border-emerald-200">Set: [${set.join(', ')}]</span>`;
}

['runAuditBtn','labPreset','labGt','labAnc'].forEach(id => document.getElementById(id).addEventListener(id === 'runAuditBtn' ? 'click' : 'change', runAudit));
['runSimBtn','simProb','simLlr','simAf','simMoi'].forEach(id => document.getElementById(id).addEventListener(id === 'runSimBtn' ? 'click' : (id === 'simMoi' ? 'change' : 'input'), runSimulator));
document.getElementById('copyAuditBtn').addEventListener('click', () => copyText(document.getElementById('auditJsonOut').textContent, 'Audit JSON'));
document.getElementById('copySimBtn').addEventListener('click', () => copyText(document.getElementById('simJsonOut').textContent, 'Simulator JSON'));
document.getElementById('copyQuickstartBtn').addEventListener('click', () => copyText(document.getElementById('quickstartCode').textContent, 'Quickstart commands'));
document.getElementById('copyBibtexBtn').addEventListener('click', () => copyText(document.getElementById('bibtexCode').textContent, 'BibTeX'));

renderVersion('genesis'); renderAlgos('ALL'); renderBugs(); renderStaticTables(); renderDiscoveries(); initLab();
['rubric','auroc','ece','mave','mondrian','ram'].forEach((k, i) => bindHoverPause(k, ['chartRubric','chartAuroc','chartEce','chartMave','chartMondrian','chartRam'][i]));
