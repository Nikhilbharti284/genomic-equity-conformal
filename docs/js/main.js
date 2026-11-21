(function(){
const r=document.documentElement,s=localStorage.getItem("gec-theme"),d=window.matchMedia&&window.matchMedia("(prefers-color-scheme: dark)").matches;
r.setAttribute("data-theme",s||(d?"dark":"light"));
const tb=document.getElementById("theme-toggle");
function sync(){if(tb)tb.textContent=r.getAttribute("data-theme")==="dark"?"☀ Light":"☾ Dark";}
sync();
if(tb)tb.addEventListener("click",()=>{const n=r.getAttribute("data-theme")==="dark"?"light":"dark";r.setAttribute("data-theme",n);localStorage.setItem("gec-theme",n);sync();window.dispatchEvent(new CustomEvent("gec-theme-change"));});
const dr=document.getElementById("who-drawer"),op=document.getElementById("menu-open"),cl=document.getElementById("menu-close");
if(op&&dr)op.addEventListener("click",()=>{dr.classList.add("open");document.body.style.overflow="hidden";});
if(cl&&dr)cl.addEventListener("click",()=>{dr.classList.remove("open");document.body.style.overflow="";});
document.addEventListener("keydown",e=>{if(e.key==="Escape"&&dr){dr.classList.remove("open");document.body.style.overflow="";}});
const bar=document.getElementById("progress-bar"),top=document.getElementById("back-to-top");
window.addEventListener("scroll",()=>{const h=document.documentElement.scrollHeight-window.innerHeight;if(bar&&h>0)bar.style.width=Math.min(100,(window.scrollY/h)*100)+"%";if(top)top.classList.toggle("show",window.scrollY>320);},{passive:true});
if(top)top.addEventListener("click",()=>window.scrollTo({top:0,behavior:"smooth"}));
document.addEventListener("click",e=>{const b=e.target.closest("[data-copy]");if(!b)return;const sel=b.getAttribute("data-copy");const txt=sel.startsWith("#")?(document.querySelector(sel)?.textContent||""):sel;navigator.clipboard.writeText(txt).then(()=>{const o=b.textContent;b.textContent="Copied!";setTimeout(()=>b.textContent=o,1400);});});
document.querySelectorAll(".acc-head").forEach(h=>h.addEventListener("click",()=>h.parentElement.classList.toggle("open")));
const af=document.getElementById("audit-filter");
if(af)af.addEventListener("input",()=>{const q=af.value.trim().toLowerCase();document.querySelectorAll(".acc-item").forEach(i=>i.style.display=i.textContent.toLowerCase().includes(q)?"":"none");});
const st=document.querySelectorAll(".pipe-stage"),pd=document.getElementById("pipe-detail");
if(st.length&&pd)st.forEach(x=>x.addEventListener("click",()=>{st.forEach(y=>y.classList.remove("active"));x.classList.add("active");pd.innerHTML=x.getAttribute("data-detail");}));
const si=document.getElementById("drawer-search"),sb=document.getElementById("drawer-results");let idx=null;
if(si&&sb){si.addEventListener("focus",async()=>{if(!idx){try{idx=await(await fetch("data/search-index.json")).json();}catch(_){idx=[];}}});
si.addEventListener("input",()=>{const q=si.value.trim().toLowerCase();if(!q||!idx){sb.innerHTML="";return;}const hits=idx.filter(p=>(p.title+" "+p.text).toLowerCase().includes(q)).slice(0,6);sb.innerHTML=hits.length?hits.map(h=>`<a href="${h.url}"><strong>${h.title}</strong><small>${h.snippet}</small></a>`).join(""):'<div style="padding:10px 0;color:#9ec4d4;font-size:.88rem">No matching pages.</div>';});}
})();