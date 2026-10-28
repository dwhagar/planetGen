// UX.204: the collapsible sections of a page (`<details>`) start closed. A link to one (`#habitability`, a contents
// entry, a symbol that points at a body's Habitability group) opens it, and so opens any collapsed section around the
// target, then brings the target into view. A section the visitor opened is left alone.
(function () {
  function reveal() {
    const id = decodeURIComponent(location.hash.slice(1));
    if (!id) return;
    const target = document.getElementById(id);
    if (!target) return;
    let opened = false;
    for (let el = target; el; el = el.parentElement && el.parentElement.closest("details")) {
      if (el.tagName === "DETAILS" && !el.open) {
        el.open = true;
        opened = true;
      }
    }
    if (opened) target.scrollIntoView();
  }
  window.addEventListener("hashchange", reveal);
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", reveal);
  else reveal();
})();
