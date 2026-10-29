// static/adminjobstrip.js
//
// The job strip on every admin page (partials/admin_jobstrip.html, ADM.119, ADM.123, PERF.154): the running job with its
// stage, bar, elapsed time and a Cancel button, and under it the waiting list in the order it arrived: each waiting job
// with Discard (and Start now once the job before it failed), each waiting run with its place. Kept current from `/admin/generate/status` (the `strip` field). It polls every few seconds
// while the tab is visible and not at all while it is hidden. When the job ends the strip says so and links to the
// job; the next page load hides it. Without script the strip is the jobs as of the page load.

export const POLL_MS = 3000;

// A button that posts `action` for the waiting entry `ticket` to the strip's form address.
function form(parts, action, ticket, label) {
  const element = document.createElement("form");
  element.setAttribute("method", "post");
  element.setAttribute("action", parts.root.dataset.actionUrl);
  element.setAttribute("class", "admin-jobstrip-form");
  for (const [name, value] of [[parts.root.dataset.csrfField, parts.root.dataset.csrfToken], ["action", action], ["ticket", ticket]]) {
    const input = document.createElement("input");
    input.setAttribute("type", "hidden");
    input.setAttribute("name", name);
    input.setAttribute("value", value);
    element.appendChild(input);
  }
  const button = document.createElement("button");
  button.setAttribute("type", "submit");
  button.textContent = label;
  element.appendChild(button);
  return element;
}

// One waiting entry of the list: its title, what it waits for, and its buttons.
function entryItem(parts, entry) {
  const item = document.createElement("li");
  const title = document.createElement("strong");
  title.textContent = entry.title;
  const note = document.createElement("span");
  note.setAttribute("class", "admin-jobstrip-text");
  note.textContent = entry.note;
  item.appendChild(title);
  item.appendChild(note);
  if (entry.startable) item.appendChild(form(parts, "start_queued", entry.ticket, "Start now"));
  if (entry.discardable) item.appendChild(form(parts, "discard_queued", entry.ticket, "Discard"));
  return item;
}

// Fill the strip's parts from the status JSON's `strip` (`{job, line}`, or `null` when nothing runs or waits).
// `previous` is the job shown last time; returns the job shown now.
export function apply(strip, parts, previous) {
  const job = strip ? strip.job : null;
  const line = strip && strip.line ? strip.line : [];
  const queued = line.length > 0;
  parts.queued.hidden = !queued;
  parts.line.replaceChildren(...line.map((entry) => entryItem(parts, entry)));
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
    cancel: find("cancel"), jobId: find("job"), queued: find("queued"), line: find("line"),
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
