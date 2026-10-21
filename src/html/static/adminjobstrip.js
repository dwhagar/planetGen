// static/adminjobstrip.js
//
// The job strip on every admin page (partials/admin_jobstrip.html, ADM.119): the running job with its stage, bar,
// elapsed time and a Cancel button, kept current from `/admin/generate/status` (the `strip` field). It polls every
// few seconds while the tab is visible and not at all while it is hidden. When the job ends the strip says so and
// links to the job; the next page load hides it. Without script the strip is the job as of the page load.

export const POLL_MS = 3000;

// Fill the strip's parts from the status JSON's `strip` object; `null` means nothing is running.
export function apply(strip, parts, previous) {
  if (!strip) {
    if (!previous) { parts.root.hidden = true; return null; }
    parts.text.textContent = "Finished.";
    parts.bar.hidden = true;
    parts.cancel.hidden = true;
    parts.root.hidden = false;
    return null;
  }
  parts.root.hidden = false;
  parts.title.textContent = strip.title;
  parts.title.setAttribute("href", strip.url);
  parts.text.textContent = strip.text;
  parts.elapsed.textContent = strip.elapsed_text;
  parts.jobId.value = strip.id;
  parts.cancel.hidden = false;
  parts.bar.hidden = false;
  if (strip.value === null || strip.value === undefined) parts.bar.removeAttribute("value");
  else parts.bar.value = strip.value;
  return strip;
}

export function init(root) {
  const parts = {
    root,
    title: root.querySelector("[data-strip-title]"),
    text: root.querySelector("[data-strip-text]"),
    bar: root.querySelector("[data-strip-bar]"),
    elapsed: root.querySelector("[data-strip-elapsed]"),
    cancel: root.querySelector("[data-strip-cancel]"),
    jobId: root.querySelector("[data-strip-job]"),
  };
  let previous = root.hidden ? null : { id: parts.jobId.value };
  let timer = null;

  async function poll() {
    timer = null;
    try {
      const response = await fetch(root.dataset.statusUrl, { credentials: "same-origin", headers: { Accept: "application/json" } });
      if (response.ok) previous = apply((await response.json()).strip, parts, previous);
    } catch {
      // a dropped connection: try again at the next tick
    }
    schedule();
  }

  function schedule() {
    if (timer === null && !document.hidden) timer = setTimeout(poll, POLL_MS);
  }

  document.addEventListener("visibilitychange", () => {
    if (document.hidden) { clearTimeout(timer); timer = null; } else poll();
  });
  parts.cancel.addEventListener("submit", (event) => {
    if (!window.confirm(`Cancel "${parts.title.textContent}"?`)) event.preventDefault();
  });
  schedule();
}

if (typeof document !== "undefined") {
  const root = document.querySelector("[data-job-strip]");
  if (root) init(root);
}
