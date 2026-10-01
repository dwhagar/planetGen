// tests/js/generatejobs.test.mjs -- static/generatejobs.js, the Generate
// page's live "Current job" panel (TEST.58): what it asks for, how it
// redraws, and when it polls again, stops or reloads.

import { test } from "node:test";
import assert from "node:assert/strict";

import { h, importPage, installDom, jsonResponse, settle } from "./fakedom.mjs";

function panel(attrs) {
  const el = h("section", Object.assign({
    id: "current-job", "data-job": "20261001-abc", "data-status-url": "/admin/generate/status",
  }, attrs || {}), [
    h("h2", { "data-job-title": "" }, "Job"),
    h("span", { "data-job-status": "", class: "badge job-running" }, "running"),
    h("p", { "data-job-step": "" }),
    h("progress", { "data-job-progress": "" }),
    h("span", { "data-job-progress-text": "" }),
    h("span", { "data-job-progress-detail": "", hidden: true }),
    h("span", { "data-job-elapsed": "" }),
    h("span", { "data-job-remaining": "" }),
    h("p", { "data-job-error": "", hidden: true }),
    h("pre", { "data-job-log": "" }, "old log"),
    h("button", { "data-job-cancel": "" }, "Cancel"),
  ]);
  document.body.appendChild(el);
  return el;
}

function job(fields) {
  return Object.assign({
    title: "Generate 4 sectors", status: "running", step: 1, steps: [{}, {}], step_label: "Sectors",
    progress_total: 4, progress_completed: 1, progress_text: "1 of 4 sectors", progress_detail_text: "",
    elapsed_text: "12 s", remaining_text: "36 s", error: null, finished: false,
  }, fields || {});
}

async function start(win, attrs) {
  const el = panel(attrs);
  await importPage("generatejobs.js");
  return el;
}

test("asks for the job's status two seconds in and draws it", async () => {
  const win = installDom("http://localhost/admin/generate");
  win.fetchHandler = () => ({ job: job(), log: "line 1\nline 2" });
  const el = await start(win);
  assert.equal(win.fetchCalls.length, 0, "nothing is fetched before the first interval");
  win.timers.advance(2000);
  await settle();
  assert.equal(win.fetchCalls.length, 1);
  const url = new URL(win.fetchCalls[0].url);
  assert.equal(url.pathname, "/admin/generate/status");
  assert.equal(url.searchParams.get("job"), "20261001-abc");
  assert.equal(win.fetchCalls[0].options.cache, "no-store");

  assert.equal(el.querySelector("[data-job-title]").textContent, "Generate 4 sectors");
  const status = el.querySelector("[data-job-status]");
  assert.equal(status.textContent, "running");
  assert.equal(status.className, "badge job-running");
  assert.equal(el.querySelector("[data-job-step]").textContent, "Step 1 of 2: Sectors");
  const bar = el.querySelector("[data-job-progress]");
  assert.equal(Number(bar.max), 4);
  assert.equal(bar.value, 1);
  assert.equal(el.querySelector("[data-job-progress-text]").textContent, "1 of 4 sectors");
  assert.equal(el.querySelector("[data-job-progress-detail]").hidden, true);
  assert.equal(el.querySelector("[data-job-elapsed]").textContent, "12 s");
  assert.equal(el.querySelector("[data-job-remaining]").textContent, "about 36 s left");
  assert.equal(el.querySelector("[data-job-error]").hidden, true);
  assert.equal(el.querySelector("[data-job-log]").textContent, "line 1\nline 2");
  assert.ok(el.querySelector("[data-job-cancel]"), "a running job keeps its Cancel button");
});

test("keeps polling every two seconds while the job runs", async () => {
  const win = installDom("http://localhost/admin/generate");
  let completed = 0;
  win.fetchHandler = () => ({ job: job({ progress_completed: ++completed }), log: "" });
  const el = await start(win);
  for (let n = 1; n <= 3; n++) {
    win.timers.advance(2000);
    await settle();
    assert.equal(win.fetchCalls.length, n);
    assert.equal(el.querySelector("[data-job-progress]").value, n);
  }
  assert.deepEqual(win.timers.pending(), [2000]);
});

test("a failed request waits three intervals before trying again", async () => {
  const win = installDom("http://localhost/admin/generate");
  let calls = 0;
  win.fetchHandler = () => (++calls === 1 ? jsonResponse({}, 503) : { job: job(), log: "" });
  await start(win);
  win.timers.advance(2000);
  await settle();
  assert.deepEqual(win.timers.pending(), [6000]);
  win.timers.advance(5999);
  await settle();
  assert.equal(calls, 1);
  win.timers.advance(1);
  await settle();
  assert.equal(calls, 2);
});

test("a network error also backs off instead of stopping", async () => {
  const win = installDom("http://localhost/admin/generate");
  win.fetchHandler = () => { throw new TypeError("network down"); };
  await start(win);
  win.timers.advance(2000);
  await settle();
  assert.deepEqual(win.timers.pending(), [6000]);
});

test("a finished job drops Cancel, fills the bar and reloads the Generate page", async () => {
  const win = installDom("http://localhost/admin/generate");
  win.fetchHandler = () => ({
    job: job({ status: "succeeded", finished: true, progress_total: 0, remaining_text: "", error: null }), log: "done",
  });
  const el = await start(win);
  win.timers.advance(2000);
  await settle();
  assert.equal(el.querySelector("[data-job-cancel]"), null);
  const bar = el.querySelector("[data-job-progress]");
  assert.equal(Number(bar.max), 1);
  assert.equal(bar.value, 1);
  assert.equal(el.querySelector("[data-job-remaining]").textContent, "");
  assert.equal(win.location.reloads, 0);
  win.timers.advance(1000);
  assert.equal(win.location.reloads, 1);
  assert.deepEqual(win.timers.pending(), [], "no more polling after the job finished");
});

test("a job's own page (pinned) stops polling without reloading", async () => {
  const win = installDom("http://localhost/admin/generate/jobs/20261001-abc");
  win.fetchHandler = () => ({ job: job({ status: "failed", finished: true, error: "Step 2 failed" }), log: "" });
  const el = await start(win, { "data-job-pinned": "" });
  win.timers.advance(2000);
  await settle();
  const error = el.querySelector("[data-job-error]");
  assert.equal(error.hidden, false);
  assert.equal(error.textContent, "Step 2 failed");
  assert.equal(el.querySelector("[data-job-status]").className, "badge job-failed");
  win.timers.advance(10000);
  assert.equal(win.location.reloads, 0);
  assert.equal(win.fetchCalls.length, 1);
});

test("an answer without a job stops polling", async () => {
  const win = installDom("http://localhost/admin/generate");
  win.fetchHandler = () => ({ job: null });
  await start(win);
  win.timers.advance(2000);
  await settle();
  assert.deepEqual(win.timers.pending(), []);
});

test("an already finished job, or no job at all, is never polled", async () => {
  let win = installDom("http://localhost/admin/generate");
  await start(win, { "data-job-finished": "" });
  assert.deepEqual(win.timers.pending(), []);

  win = installDom("http://localhost/admin/generate");
  await start(win, { "data-job": null });
  assert.deepEqual(win.timers.pending(), []);

  win = installDom("http://localhost/admin/generate");
  await importPage("generatejobs.js");
  assert.deepEqual(win.timers.pending(), []);
});

test("unknown progress shows an indeterminate bar; detail text shows when given", async () => {
  const win = installDom("http://localhost/admin/generate");
  win.fetchHandler = () => ({ job: job({ progress_total: 0, progress_detail_text: "sector 3: planets" }), log: "" });
  const el = await start(win);
  el.querySelector("[data-job-progress]").value = 0.5;
  win.timers.advance(2000);
  await settle();
  assert.equal(el.querySelector("[data-job-progress]").hasAttribute("value"), false);
  const detail = el.querySelector("[data-job-progress-detail]");
  assert.equal(detail.hidden, false);
  assert.equal(detail.textContent, "sector 3: planets");
  assert.equal(el.querySelector("[data-job-step]").textContent, "Step 1 of 2: Sectors");
});

test("the log follows new output only while it is scrolled to the bottom", async () => {
  const win = installDom("http://localhost/admin/generate");
  let n = 0;
  win.fetchHandler = () => ({ job: job(), log: "line " + ++n });
  const el = await start(win);
  const pre = el.querySelector("[data-job-log]");
  // At the bottom: the new text is followed.
  pre.scrollHeight = 500;
  pre.clientHeight = 100;
  pre.scrollTop = 400;
  win.timers.advance(2000);
  await settle();
  assert.equal(pre.textContent, "line 1");
  assert.equal(pre.scrollTop, 500);
  // Scrolled up to read: left where it is.
  pre.scrollTop = 100;
  win.timers.advance(2000);
  await settle();
  assert.equal(pre.textContent, "line 2");
  assert.equal(pre.scrollTop, 100);
});
