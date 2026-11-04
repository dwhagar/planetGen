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

  // UX.142: the shell (rail, search, density switch) stays in view under the held site header while a long admin page
  // scrolls. The header's height is not a constant, so it is measured and handed to the CSS as --site-header-h.
  var header = document.querySelector(".site-header");
  if (!header) return;
  function holdUnderHeader() {
    document.documentElement.style.setProperty("--site-header-h", header.getBoundingClientRect().height + "px");
  }
  holdUnderHeader();
  if (typeof ResizeObserver === "function") new ResizeObserver(holdUnderHeader).observe(header);
  else window.addEventListener("resize", holdUnderHeader);
})();
