let chartInstances = {};
let chartLoopEnabled = true;
const hoveredCharts = new Set();

function getStaggeredAnimation() {
  return {
    duration: 1400,
    easing: 'easeOutQuart',
    delay: (ctx) => (ctx.type === 'data' && ctx.mode === 'default') ? (ctx.dataIndex * 75 + ctx.datasetIndex * 120) : 0
  };
}

function bindHoverPause(key, canvasId) {
  const el = document.getElementById(canvasId);
  if (!el || el.dataset.hoverBound) return;
  el.dataset.hoverBound = '1';
  el.addEventListener('mouseenter', () => hoveredCharts.add(key));
  el.addEventListener('mouseleave', () => hoveredCharts.delete(key));
}

function replayActiveTabCharts(tabId) {
  const keys = tabId === 'evolution' ? ['rubric', 'auroc'] : (tabId === 'benchmarks' ? ['ece', 'mave', 'mondrian', 'ram'] : []);
  keys.forEach((k, idx) => {
    const ch = chartInstances[k];
    if (!ch || hoveredCharts.has(k)) return;
    setTimeout(() => {
      if (!chartLoopEnabled || hoveredCharts.has(k) || !chartInstances[k]) return;
      try { chartInstances[k].reset(); chartInstances[k].update(); } catch (e) {}
    }, idx * 220);
  });
}

function getCurrentActiveTab() {
  const active = document.querySelector('.tab-panel.active');
  return active ? active.id.replace('tab-', '') : 'overview';
}

setInterval(() => {
  if (!chartLoopEnabled || document.hidden) return;
  const tab = getCurrentActiveTab();
  if (tab === 'evolution' || tab === 'benchmarks') replayActiveTabCharts(tab);
}, 4500);

const loopBtn = document.getElementById('toggleChartLoopBtn');
if (loopBtn) {
  loopBtn.addEventListener('click', () => {
    chartLoopEnabled = !chartLoopEnabled;
    loopBtn.textContent = chartLoopEnabled ? 'Loop Charts: ON' : 'Loop Charts: PAUSED';
    loopBtn.className = chartLoopEnabled
      ? 'px-2.5 py-0.5 rounded-full bg-sky-500/20 hover:bg-sky-500/30 text-sky-200 border border-sky-400/40 font-semibold transition cursor-pointer'
      : 'px-2.5 py-0.5 rounded-full bg-slate-600/40 hover:bg-slate-600/60 text-slate-200 border border-slate-400/40 font-semibold transition cursor-pointer';
    showToast(chartLoopEnabled ? 'Chart loop animation resumed' : 'Chart loop animation paused');
    if (chartLoopEnabled) replayActiveTabCharts(getCurrentActiveTab());
  });
}

function initEvolutionCharts() {
  if (!chartsInitialized.rubric) {
    chartsInitialized.rubric = true;
    chartInstances.rubric = new Chart(document.getElementById('chartRubric'), {
      type: 'bar',
      data: {
        labels: ['Track 7: Genomic Equity (BUILT)', 'Track 2: AMR Gram-Neg (Backup)', 'Track 8: Drug Synergy', 'Track 1: Audited Fungal Pma1', 'Track 3: WHO Fungi', 'Track 4: MDR/XDR TB', 'Track 5: NTD Proteasome', 'Track 6: Viral RNA'],
        datasets: [{ label: 'Composite Score (N × F × Fu × I, max 625)', data: [500, 500, 400, 400, 320, 240, 192, 192], backgroundColor: ['#059669', '#0f4c5c', '#64748b', '#64748b', '#94a3b8', '#94a3b8', '#cbd5e1', '#cbd5e1'], borderRadius: 6 }]
      },
      options: { indexAxis: 'y', responsive: true, animation: getStaggeredAnimation(), maintainAspectRatio: false, plugins: { legend: { display: false } }, scales: { x: { max: 625, title: { display: true, text: 'Rubric Composite Score (Max 625)' } } } }
    });
  }
  if (!chartsInitialized.auroc) {
    chartsInitialized.auroc = true;
    chartInstances.auroc = new Chart(document.getElementById('chartAuroc'), {
      type: 'bar',
      data: {
        labels: ['Raw ESM-2 (8M)', 'v0.1.0 V1 Split', 'v0.1.0 V1 5-Fold', 'v0.1.0 V2 Platt', 'v1.0-clean (20g)', 'AlphaMissense Only', 'v1.0.0 Track B', 'Post-2023 Holdout', 'v2/v3 Pan-Global'],
        datasets: [{ label: 'OOF / Prospective AUROC', data: [0.6673, 0.7322, 0.7896, 0.8038, 0.8400, 0.9583, 0.9839, 0.9846, 0.9844], backgroundColor: ['#94a3b8', '#64748b', '#475569', '#0284c7', '#0f4c5c', '#d97706', '#059669', '#047857', '#065f46'], borderRadius: 6 }]
      },
      options: { responsive: true, maintainAspectRatio: false, animation: getStaggeredAnimation(), plugins: { legend: { display: false } }, scales: { y: { min: 0.55, max: 1.0, title: { display: true, text: 'AUROC (Zero Leakage)' } } } }
    });
  }
}

function renderEceChart(mode) {
  document.querySelectorAll('#eceModeBtns .pill-btn').forEach(b => b.classList.toggle('active', b.dataset.ecemode === mode));
  const cfg = {
    v3: { labels: ['Overall', 'SAS', 'AFR', 'EAS', 'AMR', 'MID*', 'NFE', 'FOUNDER', 'ULTRA_RARE'], d1Label: 'Raw AlphaMissense ECE', d1: [0.0686, 0.0530, 0.0543, 0.0685, 0.0480, 0.1127, 0.0941, 0.0669, 0.1333], d2Label: 'v2/v3 Nested OOF Platt ECE (−90.4%)', d2: [0.0066, 0.0229, 0.0141, 0.0335, 0.0288, 0.0566, 0.0180, 0.0179, 0.0054] },
    v1: { labels: ['Overall (9,142)', 'EUR (2,209)', 'SAS (1,098)', 'AFR (1,412)', 'ULTRA_RARE (4,423)'], d1Label: 'Raw AlphaMissense ECE', d1: [0.0686, 0.0622, 0.0527, 0.0637, 0.1192], d2Label: 'v1.0.0 Nested OOF Platt ECE (−82.2%)', d2: [0.0122, 0.0343, 0.0298, 0.0233, 0.0072] },
    v01: { labels: ['EUR (n=311)', 'SAS (n=95)', 'AFR (n=150)', 'ULTRA_RARE (n=650)'], d1Label: 'V1 Isotonic Regression (Overfit — Rejected)', d1: [0.196, 0.183, 0.210, 0.055], d2Label: 'V2 Regularized MLP + Platt (Adopted)', d2: [0.0879, 0.1324, 0.1136, 0.0211] }
  }[mode];
  if (chartInstances.ece) chartInstances.ece.destroy();
  chartInstances.ece = new Chart(document.getElementById('chartEce'), {
    type: 'bar',
    data: { labels: cfg.labels, datasets: [{ label: cfg.d1Label, data: cfg.d1, backgroundColor: '#dc2626', borderRadius: 5 }, { label: cfg.d2Label, data: cfg.d2, backgroundColor: '#0f4c5c', borderRadius: 5 }] },
    options: { responsive: true, maintainAspectRatio: false, animation: getStaggeredAnimation(), scales: { y: { beginAtZero: true, title: { display: true, text: 'Expected Calibration Error (Lower is Better)' } } } }
  });
}

function renderMaveChart(mode) {
  document.querySelectorAll('#maveModeBtns .pill-btn').forEach(b => b.classList.toggle('active', b.dataset.mavemode === mode));
  const cfg = mode === 'v2' ? {
    labels: ['TP53 (WIN)', 'BRCA1 (WIN)', 'PTEN (WIN)', 'MSH2 (WIN)', 'LDLR (MISS)'],
    d1Label: 'Raw ESM2-35M |ρ|', d1: [0.2485, 0.0085, 0.1669, 0.1635, 0.1532],
    d2Label: 'v2.0.0 ESM2-35M + AFDB v6 C-Beta RSA', d2: [0.2488, 0.0124, 0.1698, 0.1665, 0.1530],
    note: 'v2.0.0 Upgrade Proof: Replacing aggressive pLDDT scaling with pure-NumPy 10Å C-Beta contact density flipped TP53 from miss to WIN (+0.0003). LDLR (−0.0002) is transparently reported as an honest miss (4/5 wins).'
  } : {
    labels: ['PTEN (WIN)', 'MSH2 (WIN)', 'LDLR (WIN)', 'BRCA1 (WIN)', 'TP53 (MISS)'],
    d1Label: 'Raw ESM2-35M |ρ|', d1: [0.4276, 0.3727, 0.2902, 0.0101, 0.1649],
    d2Label: 'v1.0.0 ESM2-35M + AFDB v6 pLDDT', d2: [0.4497, 0.3767, 0.2920, 0.0174, 0.1021],
    note: 'v1.0.0 Baseline: Global pLDDT modulation won 4/5 genes (PTEN +0.0221, MSH2 +0.0040, LDLR +0.0018, BRCA1 +0.0073) but underperformed on TP53 (−0.0628), motivating the v2.0.0 C-Beta RSA burial upgrade.'
  };
  document.getElementById('maveNote').textContent = cfg.note;
  if (chartInstances.mave) chartInstances.mave.destroy();
  chartInstances.mave = new Chart(document.getElementById('chartMave'), {
    type: 'bar',
    data: { labels: cfg.labels, datasets: [{ label: cfg.d1Label, data: cfg.d1, backgroundColor: '#94a3b8', borderRadius: 5 }, { label: cfg.d2Label, data: cfg.d2, backgroundColor: '#059669', borderRadius: 5 }] },
    options: { responsive: true, maintainAspectRatio: false, animation: getStaggeredAnimation(), scales: { y: { beginAtZero: true, title: { display: true, text: 'Spearman |ρ| (Zero ClinVar Labels)' } } } }
  });
}

function initBenchmarkCharts() {
  if (!chartsInitialized.bench) {
    chartsInitialized.bench = true;
    renderEceChart('v3');
    renderMaveChart('v2');
    chartInstances.mondrian = new Chart(document.getElementById('chartMondrian'), {
      type: 'bar',
      data: {
        labels: ['SAS (700)', 'AFR (1,221)', 'EAS (740)', 'AMR (717)', 'MID* (122)', 'NFE (1,192)', 'FOUNDER (622)', 'ULTRA_RARE (3,828)'],
        datasets: [
          { label: 'Mondrian Coverage (Target >= 0.90)', data: [0.9114, 0.9328, 0.9149, 0.9121, 0.8934, 0.9077, 0.9277, 0.9875], backgroundColor: '#0f4c5c', borderRadius: 5 },
          { label: 'Decisive Singleton Rate', data: [0.979, 1.000, 0.958, 0.943, 1.000, 0.905, 1.000, 1.000], backgroundColor: '#d97706', borderRadius: 5 }
        ]
      },
      options: { responsive: true, maintainAspectRatio: false, animation: getStaggeredAnimation(), scales: { y: { min: 0.80, max: 1.02, title: { display: true, text: 'Proportion (0 Empty Sets)' } } } }
    });
    chartInstances.ram = new Chart(document.getElementById('chartRam'), {
      type: 'bar',
      data: {
        labels: ['ClinVar Stream', 'AM 71M Stream', '10K WES v3.0', 'ESM2-8M CPU', 'ESM2-35M Cache', 'Free Tier Cap'],
        datasets: [{ label: 'Peak RAM (MB)', data: [34.0, 55.25, 181.1, 507.06, 629.0, 8192], backgroundColor: ['#059669', '#059669', '#0284c7', '#0f4c5c', '#0f4c5c', '#cbd5e1'], borderRadius: 6 }]
      },
      options: { responsive: true, maintainAspectRatio: false, animation: getStaggeredAnimation(), plugins: { legend: { display: false } }, scales: { y: { type: 'logarithmic', title: { display: true, text: 'Peak RAM (MB, log scale)' } } } }
    });
  }
}

document.querySelectorAll('#eceModeBtns .pill-btn').forEach(btn => btn.addEventListener('click', () => renderEceChart(btn.dataset.ecemode)));
document.querySelectorAll('#maveModeBtns .pill-btn').forEach(btn => btn.addEventListener('click', () => renderMaveChart(btn.dataset.mavemode)));
