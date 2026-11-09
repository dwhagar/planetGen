// static/adminshell.js
//
// The admin shell (base.html, ADM.92): keeps the shell under the held site header.

(function () {
  // UX.142: the shell (rail, search) stays in view under the held site header while a long admin page
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
