// static/adminshell.js
//
// The admin shell (base.html, ADM.92): the density switch. Compact is the default (the CSS styles `.app` unless
// `data-density="roomy"`); the choice is kept in this browser. It is a classic script, loaded right after the shell's
// markup, so the page is drawn at the chosen density without a flash. Without script the page is compact and the
// switch does nothing.

(function () {
  var KEY = "planetgen.admin.density";
  var box = document.querySelector("[data-admin-density]");
  var stored = null;
  try {
    stored = window.localStorage.getItem(KEY);
  } catch (e) {
    // private window or blocked storage: compact, nothing remembered
  }
  function apply(compact) {
    if (compact) document.body.removeAttribute("data-density");
    else document.body.setAttribute("data-density", "roomy");
    if (box) box.checked = compact;
  }
  apply(stored !== "roomy");
  if (!box) return;
  box.addEventListener("change", function () {
    apply(box.checked);
    try {
      window.localStorage.setItem(KEY, box.checked ? "compact" : "roomy");
    } catch (e) {
      // not remembered
    }
  });
})();
