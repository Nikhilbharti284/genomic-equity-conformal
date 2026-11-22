function toggleDrawer(open){
  var d=document.getElementById("nav-drawer");
  if(d) d.classList.toggle("open", open);
}
function copyText(id, btn){
  var el=document.getElementById(id);
  if(!el) return;
  navigator.clipboard.writeText(el.textContent).then(function(){
    var old=btn.textContent; btn.textContent="Copied!";
    setTimeout(function(){ btn.textContent=old; }, 1600);
  });
}
