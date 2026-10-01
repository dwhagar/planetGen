// static/generatejobs.js
//
// Keeps the admin Generate page's "Current job" panel live
// (web/templates/generate.html and generate_job.html): every few seconds
// it asks /admin/generate/status for the job's state and the tail of its
// output, and redraws the step, progress bar, elapsed time, error and
// log in place. When the job finishes, the Generate page reloads itself
// (so the galaxy summary, the job list and the forms catch up); a job's
// own page just stops polling. Without JavaScript the page still works:
// reload it to see progress.

const POLL_MS = 2000;

const panel = document.getElementById("current-job");

function setText(selector, text) {
  const el = panel.querySelector(selector);
  if (el && el.textContent !== text) el.textContent = text;
}

function render(job, log) {
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
  setText("[data-job-remaining]", job.remaining_text ? `about ${job.remaining_text} left` : "");

  const error = panel.querySelector("[data-job-error]");
  if (error) {
    error.textContent = job.error || "";
    error.hidden = !job.error;
  }

  const pre = panel.querySelector("[data-job-log]");
  if (pre && typeof log === "string" && pre.textContent !== log) {
    const atBottom = pre.scrollHeight - pre.scrollTop - pre.clientHeight < 24;
    pre.textContent = log;
    if (atBottom) pre.scrollTop = pre.scrollHeight;
  }

  if (job.finished) {
    const cancel = panel.querySelector("[data-job-cancel]");
    if (cancel) cancel.remove();
  }
}

async function poll() {
  const url = new URL(panel.dataset.statusUrl, window.location.href);
  url.searchParams.set("job", panel.dataset.job);
  let body;
  try {
    const response = await fetch(url, { credentials: "same-origin", cache: "no-store" });
    if (!response.ok) throw new Error(`status ${response.status}`);
    body = await response.json();
  } catch (err) {
    window.setTimeout(poll, POLL_MS * 3);
    return;
  }
  if (!body.job) return;
  render(body.job, body.log);
  if (!body.job.finished) {
    window.setTimeout(poll, POLL_MS);
  } else if (!("jobPinned" in panel.dataset)) {
    window.setTimeout(() => window.location.reload(), 1000);
  }
}

if (panel && panel.dataset.job && !("jobFinished" in panel.dataset)) {
  const pre = panel.querySelector("[data-job-log]");
  if (pre) pre.scrollTop = pre.scrollHeight;
  window.setTimeout(poll, POLL_MS);
}
