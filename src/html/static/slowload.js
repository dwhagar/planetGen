// static/slowload.js
//
// The loading state of a slow page (web/lib/slowroute.py, UX.118). The server answers a page whose data is being
// built by a queued job (PERF.105) with the page's frame and a progress bar; this polls the job
// (`GET /api/page-jobs/<id>`), draws its estimated progress, and when the job has succeeded fetches the same address
// again and swaps the finished content in without reloading the page. If that page is waiting for another job (a
// search asks for several answers) it carries on with that one. A failed or cancelled job shows its reason and a
// Try again (the fetch carries `X-Slow-Retry`, so a new job replaces the failed one). Stop waiting only stops this
// page polling: the job is shared by every visitor who asked for the same thing, so it goes on, and Check again
// picks it up. Without JavaScript the page refreshes itself (a `<noscript>` meta refresh).
//
// The progress is an estimate, not a count (the answer is one query), and the text says so.

const VERSION_QUERY = new URL(import.meta.url).search;
// The formatters load by the script's own version (a static import cannot carry it). The buttons are wired at once
// and the polling starts when they have arrived: a top-level `await` here left Stop waiting on the screen, clicked
// and doing nothing, until both had loaded (TEST.248).
let formatNumber = null;
let formatDurationSeconds = null;
const formatters = Promise.all([
  import(`./numberformat.js${VERSION_QUERY}`),
  import(`./unitladder.js${VERSION_QUERY}`),
]).then(function ([numbers, ladder]) {
  formatNumber = numbers.formatNumber;
  formatDurationSeconds = ladder.formatDurationSeconds;
});

const POLL_MS = 1000;
const RETRY_HEADER = "X-Slow-Retry";
const STATE_TEXT = {
  queued: "Waiting for a worker to start building this page.",
  running: "Building this page.",
  failed: "This page could not be built.",
  cancelled: "An administrator stopped building this page.",
  succeeded: "Built. Showing it now.",
};

let panel = document.querySelector("[data-slow-load]");
let timer = null;
let stopped = false;
let busy = false;

function part(name) {
  return panel.querySelector(`[data-slow-${name}]`);
}

function setText(name, text) {
  const el = part(name);
  if (el && el.textContent !== text) el.textContent = text;
}

function show(name, visible) {
  const el = part(name);
  if (el) el.hidden = !visible;
}

// The bar's value, 0 to 100, or null when it can't say (queued, or the size is not known).
function percentOf(job) {
  if (job.state === "succeeded") return 100;
  if (job.state !== "running" || !job.total || job.done == null) return null;
  return Math.min(99, Math.floor((100 * job.done) / job.total));
}

function statusText(job) {
  const percent = percentOf(job);
  if (job.state !== "running" || percent === null) return STATE_TEXT[job.state] || "";
  const left = job.eta_seconds == null ? "still working" : `about ${formatDurationSeconds(job.eta_seconds)} left`;
  return `About ${percent}% (an estimate: ${formatNumber(job.done)} of ${formatNumber(job.total)} rows), ${left}.`;
}

function render(job) {
  panel.dataset.state = job.state;
  const bar = part("bar");
  const percent = percentOf(job);
  if (bar) {
    if (percent === null) bar.removeAttribute("value"); // indeterminate
    else bar.value = percent;
  }
  setText("percent", percent === null ? "" : `${percent}%`);
  const failed = job.state === "failed" || job.state === "cancelled";
  setText("status", failed ? "" : statusText(job));
  const error = part("error");
  if (error) {
    error.textContent = failed ? (job.state === "failed" ? job.error || STATE_TEXT.failed : STATE_TEXT.cancelled) : "";
    error.hidden = !failed;
  }
  show("stop", !failed && !stopped);
  show("retry", failed);
}

function schedule() {
  window.clearTimeout(timer);
  if (!stopped) timer = window.setTimeout(poll, POLL_MS);
}

async function poll() {
  if (stopped || busy) return;
  busy = true;
  try {
    const response = await fetch(panel.dataset.statusUrl, { cache: "no-store", headers: { Accept: "application/json" } });
    if (stopped) {
      busy = false; // "Stop waiting" was pressed while this answer was on its way: it must not undo the stop
      return;
    }
    if (response.status === 404) {
      busy = false;
      await bring(false); // the job is gone (finished long ago, or expired): ask for the page itself
      return;
    }
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    const job = await response.json();
    if (stopped) {
      busy = false;
      return;
    }
    render(job);
    if (job.state === "succeeded") {
      busy = false;
      await bring(false);
      return;
    }
    if (job.state === "failed" || job.state === "cancelled") {
      busy = false;
      return;
    }
    setText("status", statusText(job));
  } catch (error) {
    setText("status", "Can't reach the server; trying again.");
  }
  busy = false;
  schedule();
}

// Fetches the page at this address and puts what it answers in place: the finished content, or, when it is
// waiting for another job, that job's progress.
async function bring(retry) {
  busy = true;
  try {
    const headers = { Accept: "text/html" };
    if (retry) headers[RETRY_HEADER] = "1";
    const response = await fetch(window.location.href, { cache: "no-store", headers });
    const text = await response.text();
    if (response.status !== 200 && response.status !== 202) {
      fail(`The page could not be built (HTTP ${response.status}).`);
      return;
    }
    const next = new DOMParser().parseFromString(text, "text/html");
    const waiting = next.querySelector("[data-slow-load]");
    if (response.status === 202 && waiting) {
      panel.dataset.job = waiting.dataset.job;
      panel.dataset.statusUrl = waiting.dataset.statusUrl;
      busy = false;
      stopped = false;
      await poll();
      return;
    }
    swapIn(next);
  } catch (error) {
    fail("Can't reach the server.");
  }
  busy = false;
}

function fail(message) {
  stopped = true;
  window.clearTimeout(timer);
  panel.dataset.state = "failed";
  setText("status", "");
  const error = part("error");
  if (error) {
    error.textContent = message;
    error.hidden = false;
  }
  show("stop", false);
  show("resume", false);
  show("retry", true);
}

function swapIn(next) {
  const main = document.getElementById("main");
  const fresh = next.getElementById("main");
  if (!main || !fresh) {
    window.location.reload();
    return;
  }
  document.title = next.title || document.title;
  // The finished page's own scripts (its table, say) run once, when they load: bring in the ones this page lacks.
  const have = new Set(Array.from(document.querySelectorAll("script[src]")).map((el) => el.getAttribute("src")));
  for (const script of next.querySelectorAll("script[type=module][src]")) {
    const src = script.getAttribute("src");
    if (have.has(src)) continue;
    const copy = document.createElement("script");
    copy.type = "module";
    copy.src = src;
    document.head.appendChild(copy);
  }
  main.replaceChildren(...Array.from(fresh.childNodes));
  main.focus();
}

function stop() {
  stopped = true;
  window.clearTimeout(timer);
  setText("status", "Stopped waiting. The page is still being built for anyone else who asked; check again when you like.");
  show("stop", false);
  show("resume", true);
}

function resume() {
  stopped = false;
  show("resume", false);
  show("stop", true);
  formatters.then(poll);
}

if (panel) {
  panel.querySelector("[data-slow-stop]")?.addEventListener("click", stop);
  panel.querySelector("[data-slow-resume]")?.addEventListener("click", resume);
  panel.querySelector("[data-slow-retry]")?.addEventListener("click", () => {
    stopped = false;
    show("retry", false);
    setText("status", STATE_TEXT.queued);
    formatters.then(function () { bring(true); });
  });
  const failedAtStart = panel.dataset.state === "failed" || panel.dataset.state === "cancelled";
  if (!failedAtStart) formatters.then(schedule);
}
