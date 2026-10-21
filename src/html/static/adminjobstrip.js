// static/adminjobstrip.js
//
// The job strip on every admin page (partials/admin_jobstrip.html, ADM.119, ADM.123): the running job with its stage,
// bar, elapsed time and a Cancel button, and under it the one job waiting behind it (Discard, and Start now once the
// job before it failed). Kept current from `/admin/generate/status` (the `strip` field). It polls every few seconds
// while the tab is visible and not at all while it is hidden. When the job ends the strip says so and links to the
// job; the next page load hides it. Without script the strip is the jobs as of the page load.

export const POLL_MS = 3000;

// Fill the strip's parts from the status JSON's `strip` (`{job, queued}`, or `null` when nothing runs or waits).
// `previous` is the job shown last time; returns the job shown now.
export function apply(strip, parts, previous) {
  const job = strip ? strip.job : null;
  const queued = strip ? strip.queued : null;
  parts.queued.hidden = !queued;
  if (queued) {
    parts.queuedTitle.textContent = queued.title;
    parts.queuedNote.textContent = queued.held || "Starts when the running job ends.";
    parts.start.hidden = !queued.held;
  }
  if (job) {
    parts.run.hidden = false;
    parts.title.textContent = job.title;
    parts.title.setAttribute("href", job.url);
    parts.text.textContent = job.text;
    parts.elapsed.textContent = job.elapsed_text;
    parts.jobId.value = job.id;
    parts.cancel.hidden = false;
    parts.bar.hidden = false;
    if (job.value === null || job.value === undefined) parts.bar.removeAttribute("value");
    else parts.bar.value = job.value;
    parts.root.hidden = false;
    return job;
  }
  if (previous) {
    parts.text.textContent = "Finished.";
    parts.bar.hidden = true;
    parts.cancel.hidden = true;
    parts.run.hidden = false;
    parts.root.hidden = false;
  } else {
    parts.run.hidden = true;
    parts.root.hidden = !queued;
  }
  return null;
}

export function init(root) {
  const find = (name) => root.querySelector(`[data-strip-${name}]`);
  const parts = {
    root, run: find("run"), title: find("title"), text: find("text"), bar: find("bar"), elapsed: find("elapsed"),
    cancel: find("cancel"), jobId: find("job"), queued: find("queued"), queuedTitle: find("queued-title"),
    queuedNote: find("queued-note"), start: find("start"),
  };
  let previous = parts.run.hidden ? null : { id: parts.jobId.value };
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
