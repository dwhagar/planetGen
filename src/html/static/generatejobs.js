// static/generatejobs.js
//
// Keeps the admin Generate page's "Current job" panel live
// (web/templates/generate.html and generate_job.html) from the job's
// Server-Sent Events stream (ADM.22, `/admin/generate/jobs/<id>/stream`):
// `state` events redraw the step, progress bar, elapsed time and error in
// place, `log` events feed the output terminal (static/joblog.js), and
// `done` ends it. The stream closes itself every so often (a request must
// stay under Apache's request timeout); the browser reconnects with the
// byte offset it reached, so no output is lost or repeated. When the job
// finishes, the Generate page reloads itself (unless it failed: the log then
// stays open until Continue is clicked) (so the galaxy summary, the
// job list and the forms catch up); a job's own page just stops.
// Without EventSource it polls /admin/generate/status instead, and
// without JavaScript the page still works: reload it to see progress.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { createJobLog } = await import(`./joblog.js${VERSION_QUERY}`);
const { quietReconnect } = await import(`./connectionnote.js${VERSION_QUERY}`);

const POLL_MS = 2000;
const RECONNECT_GRACE_MS = 8000;

const panel = document.getElementById("current-job");

function setText(selector, text) {
  const el = panel.querySelector(selector);
  if (el && el.textContent !== text) el.textContent = text;
}

function render(job) {
  failed = job.status === "failed";
  setText("[data-job-title]", job.title || "Job");
  const status = panel.querySelector("[data-job-status]");
  if (status) {
    status.textContent = job.status;
    status.className = `badge job-${job.status}`;
  }
  setText("[data-job-step]", job.step_label ? `Step ${job.step} of ${job.steps.length}: ${job.step_label}` : "");

  const bar = panel.querySelector("[data-job-progress]");
  if (bar) {
    if (job.progress_total) {
      bar.max = job.progress_total;
      bar.value = job.progress_completed || 0;
    } else if (job.finished) {
      bar.max = 1;
      bar.value = 1;
    } else {
      bar.removeAttribute("value"); // indeterminate
    }
  }
  setText("[data-job-progress-text]", job.progress_text || "");
  const detail = panel.querySelector("[data-job-progress-detail]");
  if (detail) {
    detail.textContent = job.progress_detail_text || "";
    detail.hidden = !job.progress_detail_text;
  }
  setText("[data-job-elapsed]", job.elapsed_text || "");
  setText("[data-job-remaining]", job.remaining_label || "");

  const error = panel.querySelector("[data-job-error]");
  if (error) {
    error.textContent = job.error || "";
    error.hidden = !job.error;
  }

  if (job.finished) {
    const cancel = panel.querySelector("[data-job-cancel]");
    if (cancel) cancel.remove();
  }
}

// A job that succeeded or was cancelled reloads the page (it is pinned to
// the job's own page, which just stops). One that failed stays as it is,
// with its error and log in view, until Continue is clicked (ADM.24).
function finish() {
  if ("jobPinned" in panel.dataset) return;
  if (failed) {
    const box = panel.querySelector("[data-job-continue]");
    const button = panel.querySelector("[data-job-continue-button]");
    if (box && button) {
      box.hidden = false;
      button.addEventListener("click", () => window.location.reload());
      button.focus();
      return;
    }
  }
  window.setTimeout(() => window.location.reload(), 1000);
}

let failed = false;

// Copy log (ADM.25): the job's whole output, fetched from the download
// route, onto the clipboard, so an error trace can be pasted into a report.
function wireCopy(url) {
  const button = panel.querySelector("[data-job-copy]");
  const note = panel.querySelector("[data-job-copy-note]");
  if (!button || !navigator.clipboard) return;
  button.hidden = false;
  const say = (text) => {
    if (note) {
      note.textContent = text;
      note.hidden = !text;
    }
  };
  button.addEventListener("click", async () => {
    try {
      const response = await fetch(url, { credentials: "same-origin", cache: "no-store" });
      if (!response.ok) throw new Error(`status ${response.status}`);
      await navigator.clipboard.writeText(await response.text());
      say("Copied.");
    } catch (err) {
      say("Could not copy the log.");
    }
  });
}

function connection(text) {
  const note = panel.querySelector("[data-job-connection]");
  if (note) {
    note.textContent = text;
    note.hidden = !text;
  }
}

function stream(url, log) {
  const source = new EventSource(url, { withCredentials: true });
  source.addEventListener("log", (event) => log.write(JSON.parse(event.data).text));
  source.addEventListener("state", (event) => render(JSON.parse(event.data)));
  source.addEventListener("done", () => {
    source.close();
    connection("");
    finish();
  });
  // EventSource reconnects by itself, resuming from the last event id; the
  // server ends every response on purpose (STREAM_SECONDS), which the browser
  // reports as an error too, so the note waits to see whether it reconnects.
  const link = quietReconnect(connection, RECONNECT_GRACE_MS);
  source.addEventListener("open", () => link.restored());
  source.addEventListener("error", () => link.lost());
}

// The polling fallback: the job's state and the tail of its output.
async function poll(pre) {
  const url = new URL(panel.dataset.statusUrl, window.location.href);
  url.searchParams.set("job", panel.dataset.job);
  let body;
  try {
    const response = await fetch(url, { credentials: "same-origin", cache: "no-store" });
    if (!response.ok) throw new Error(`status ${response.status}`);
    body = await response.json();
  } catch (err) {
    window.setTimeout(() => poll(pre), POLL_MS * 3);
    return;
  }
  if (!body.job) return;
  render(body.job);
  if (pre && typeof body.log === "string" && pre.textContent !== body.log) {
    const atBottom = pre.scrollHeight - pre.scrollTop - pre.clientHeight < 24;
    pre.textContent = body.log;
    if (atBottom) pre.scrollTop = pre.scrollHeight;
  }
  if (!body.job.finished) {
    window.setTimeout(() => poll(pre), POLL_MS);
  } else {
    finish();
  }
}

if (panel && panel.dataset.job && !("jobFinished" in panel.dataset)) {
  const pre = panel.querySelector("[data-job-log]");
  if (typeof EventSource !== "undefined" && panel.dataset.jobUrlTemplate) {
    const base = panel.dataset.jobUrlTemplate.replace("__JOB__", encodeURIComponent(panel.dataset.job));
    wireCopy(`${base}/log`);
    stream(`${base}/stream`, await createJobLog(pre));
  } else {
    if (pre) pre.scrollTop = pre.scrollHeight;
    if (panel.dataset.jobUrlTemplate) {
      wireCopy(panel.dataset.jobUrlTemplate.replace("__JOB__", encodeURIComponent(panel.dataset.job)) + "/log");
    }
    window.setTimeout(() => poll(pre), POLL_MS);
  }
}
