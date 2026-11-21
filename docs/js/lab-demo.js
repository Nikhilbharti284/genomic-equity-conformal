(function(){
  // Paste your YouTube embed URL (e.g. https://www.youtube.com/embed/VIDEO_ID) or Loom embed URL here once recorded:
  const DEMO_VIDEO_URL = "";
  const mount = document.getElementById("vcf-demo-mount");
  if (!mount) return;
  if (DEMO_VIDEO_URL) {
    mount.innerHTML = '<div style="position:relative;padding-bottom:56.25%;height:0;overflow:hidden;border-radius:10px"><iframe src="'+DEMO_VIDEO_URL+'" style="position:absolute;top:0;left:0;width:100%;height:100%;border:0" allowfullscreen></iframe></div>';
    return;
  }
  const steps = [
    {t:"0:00 / 0:60", title:"Stage 1 — Raw Hospital VCF 4.2 Uploaded", log:"$ python3 src/audit_vcf.py --vcf hospital_sample.vcf\nLoaded 3 GRCh38 variants: chr17:7675088 C>T | chr11:5227002 T>A | chr17:7676154 G>C\nStreamed AlphaMissense (CC BY 4.0) + gnomAD v4 SAS/AFR AF in 55.25 MB RAM.", pct:20},
    {t:"0:15 / 0:60", title:"Stage 2 — Homology-Safe Feature Resolution", log:"Resolved HGVS & Protein Changes:\n  1. TP53_R175H  (SAS AF: 0.000000, AFR AF: 0.000000)\n  2. HBB_E7V     (SAS AF: 0.000621, AFR AF: 0.046998 — Sickle Cell!)\n  3. TP53_P72R   (SAS AF: 0.500196, AFR AF: 0.379906)", pct:40},
    {t:"0:30 / 0:60", title:"Stage 3 — Nested OOF Platt Calibration (-82.2% ECE)", log:"Calibrated Probabilities (ECE 0.0686 -> 0.0122):\n  TP53_R175H : Platt p = 0.9202\n  HBB_E7V    : Platt p = 0.0195 (raw model suppressed by high AFR AF!)\n  TP53_P72R  : Platt p = 0.0015", pct:60},
    {t:"0:45 / 0:60", title:"Stage 4 & 5 — Mondrian Conformal + 3-Layer BA1 Safeguard", log:"[SAFEGUARD TRIGGERED] HBB_E7V matches ClinGen BA1 Exception Registry!\n  -> Blocked false Benign override; routed to VUS (Layer-1 BA1 Exception).\n  -> TP53_R175H => Pathogenic_90pct (Singleton)\n  -> TP53_P72R  => Benign_90pct (Singleton)", pct:85},
    {t:"1:00 / 0:60", title:"Complete — W3C RO-Crate SHA-256 Audit Report Ready", log:"[PASS] 3/3 Hospital VCF variants audited | 0 Empty Sets | SHA-256 Verified\nExported: raw_vcf_audit_output.json (e32c23e2a93c...4021)\nScroll below to test your own variants live in the Interactive Lab!", pct:100}
  ];
  let idx = 0, timer = null;
  mount.innerHTML = '<div style="background:#0b1f26;color:#e2e8f0;border-radius:10px;padding:14px;border:1px solid #1e4b58">'
    + '<div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:8px;margin-bottom:8px">'
    + '<strong id="demo-title" style="color:#f6c453;font-size:14px"></strong>'
    + '<span id="demo-time" style="font-size:12px;color:#94a3b8;font-family:monospace"></span></div>'
    + '<div style="height:6px;background:#1e3a44;border-radius:4px;overflow:hidden;margin-bottom:10px"><div id="demo-bar" style="height:100%;width:0%;background:#c5942f;transition:width .4s"></div></div>'
    + '<pre id="demo-log" style="margin:0;background:#071419;color:#a7f3d0;padding:12px;border-radius:8px;font-size:12px;overflow-x:auto;white-space:pre-wrap;min-height:92px"></pre>'
    + '<div style="display:flex;gap:8px;margin-top:10px;flex-wrap:wrap">'
    + '<button id="demo-play" class="btn-gold" style="padding:6px 14px;cursor:pointer">&#9654; Play 60s Walkthrough</button>'
    + '<button id="demo-next" style="padding:6px 12px;border-radius:6px;border:1px solid #475569;background:#132f38;color:#fff;cursor:pointer">Next Step &#8594;</button>'
    + '</div></div>';
  function render(){
    const s = steps[idx];
    document.getElementById("demo-title").textContent = s.title;
    document.getElementById("demo-time").textContent = s.t;
    document.getElementById("demo-bar").style.width = s.pct + "%";
    document.getElementById("demo-log").textContent = s.log;
  }
  function stop(){ if(timer){clearInterval(timer);timer=null;} document.getElementById("demo-play").innerHTML = "&#9654; Play 60s Walkthrough"; }
  document.getElementById("demo-play").onclick = function(){
    if(timer){ stop(); return; }
    if(idx === steps.length - 1) idx = 0;
    render();
    this.innerHTML = "&#10074;&#10074; Pause";
    timer = setInterval(function(){
      if(idx < steps.length - 1){ idx++; render(); } else { stop(); }
    }, 3200);
  };
  document.getElementById("demo-next").onclick = function(){ stop(); idx = (idx + 1) % steps.length; render(); };
  render();
})();