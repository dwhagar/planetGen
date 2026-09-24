// html/static/copycode.js
//
// The system page's Copy button (html/system.py's `_code_html`): copies
// the generated Wikitext/Markdown out of the read-only code box named by
// the button's `data-copy-target`. A separate file because the site's
// Content-Security-Policy (`default-src 'self'`) allows no inline script.
// The textarea stays selectable by hand, so the page still works without
// JavaScript or clipboard access.

function copyFrom(button) {
  const box = document.getElementById(button.dataset.copyTarget);
  if (!box) {
    return;
  }
  const done = (label) => {
    button.textContent = label;
    setTimeout(() => { button.textContent = "Copy"; }, 1500);
  };
  box.select();
  if (navigator.clipboard && window.isSecureContext) {
    navigator.clipboard.writeText(box.value).then(() => done("Copied"), () => done("Press Ctrl+C"));
  } else {
    // Plain-HTTP deployments have no async clipboard API.
    done(document.execCommand("copy") ? "Copied" : "Press Ctrl+C");
  }
}

for (const button of document.querySelectorAll("[data-copy-target]")) {
  button.addEventListener("click", () => copyFrom(button));
}
