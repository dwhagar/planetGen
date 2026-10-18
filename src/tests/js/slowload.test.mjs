// tests/js/slowload.test.mjs -- static/slowload.js, the loading state of a slow page (UX.118): what it asks for,
// how it draws the job's estimated progress, and what it does when the job ends, fails, is gone or is stopped.

import { test } from "node:test";
import assert from "node:assert/strict";

import { h, importPage, installDom, jsonResponse, settle } from "./fakedom.mjs";

function panel(attrs) {
  const el = h("section", Object.assign({
    "data-slow-load": "", "data-job": "0123456789abcdef", "data-status-url": "/api/page-jobs/0123456789abcdef",
    "data-state": "queued",
  }, attrs || {}), [
    h("progress", { "data-slow-bar": "", max: "100" }),
    h("span", { "data-slow-percent": "" }),
    h("p", { "data-slow-status": "" }, "Waiting for a worker to start building this page."),
    h("p", { "data-slow-error": "", hidden: true }),
    h("button", { "data-slow-stop": "" }, "Stop waiting"),
    h("button", { "data-slow-resume": "", hidden: true }, "Check again"),
    h("button", { "data-slow-retry": "", hidden: true }, "Try again"),
  ]);
  document.body.appendChild(el);
  return el;
}

function job(fields) {
  return Object.assign({
    id: "0123456789abcdef", state: "running", done: 300, total: 1000, unit: "rows", estimated: true,
    eta_seconds: 5, seconds: 3, error: null, cancelling: false,
  }, fields || {});
}

async function start(win, attrs) {
  const el = panel(attrs);
  await importPage("slowload.js");
  return el;
}

function text(el, name) {
  return el.querySelector(`[data-slow-${name}]`).textContent;
}

test("asks for the job's status a second in and draws an estimate that says so", async () => {
  const win = installDom("http://localhost/phenomena");
  win.fetchHandler = () => job();
  const el = await start(win);
  assert.equal(win.fetchCalls.length, 0, "nothing is fetched before the first interval");
  win.timers.advance(1000);
  await settle();
  assert.equal(win.fetchCalls.length, 1);
  assert.equal(win.fetchCalls[0].url, "/api/page-jobs/0123456789abcdef");
  assert.equal(win.fetchCalls[0].options.cache, "no-store");
  assert.equal(el.querySelector("[data-slow-bar]").value, 30);
  assert.equal(text(el, "percent"), "30%");
  assert.equal(text(el, "status"), "About 30% (an estimate: 300 of 1,000 rows), about 5 s left.");
  assert.equal(el.dataset.state, "running");
});

test("a job still queued leaves the bar with no value, and one past its expected time says it is still working", async () => {
  const win = installDom("http://localhost/phenomena");
  let answer = job({ state: "queued", done: 0, eta_seconds: null });
  win.fetchHandler = () => answer;
  const el = await start(win);
  win.timers.advance(1000);
  await settle();
  const bar = el.querySelector("[data-slow-bar]");
  assert.equal(bar.hasAttribute("value"), false);
  assert.equal(text(el, "percent"), "");
  assert.equal(text(el, "status"), "Waiting for a worker to start building this page.");
  answer = job({ done: 950, eta_seconds: null });
  win.timers.advance(1000);
  await settle();
  assert.equal(bar.value, 95);
  assert.match(text(el, "status"), /still working\.$/);
});

test("keeps polling every second while the job runs", async () => {
  const win = installDom("http://localhost/phenomena");
  let done = 100;
  win.fetchHandler = () => job({ done: (done += 100) });
  const el = await start(win);
  for (let n = 1; n <= 3; n++) {
    win.timers.advance(1000);
    await settle();
    assert.equal(win.fetchCalls.length, n);
  }
  assert.equal(el.querySelector("[data-slow-bar]").value, 40);
});

test("a failed job shows why and offers Try again, and stops polling", async () => {
  const win = installDom("http://localhost/phenomena");
  win.fetchHandler = () => job({ state: "failed", done: 0, error: "Lost connection to MySQL" });
  const el = await start(win);
  win.timers.advance(1000);
  await settle();
  const error = el.querySelector("[data-slow-error]");
  assert.equal(error.hidden, false);
  assert.equal(error.textContent, "Lost connection to MySQL");
  assert.equal(el.querySelector("[data-slow-retry]").hidden, false);
  assert.equal(el.querySelector("[data-slow-stop]").hidden, true);
  win.timers.advance(5000);
  await settle();
  assert.equal(win.fetchCalls.length, 1);
});

test("a cancelled job says an administrator stopped it", async () => {
  const win = installDom("http://localhost/phenomena");
  win.fetchHandler = () => job({ state: "cancelled", done: 0 });
  const el = await start(win);
  win.timers.advance(1000);
  await settle();
  assert.equal(text(el, "error"), "An administrator stopped building this page.");
  assert.equal(el.querySelector("[data-slow-retry]").hidden, false);
});

test("a page that starts out failed does not poll", async () => {
  const win = installDom("http://localhost/phenomena");
  win.fetchHandler = () => job();
  await start(win, { "data-state": "failed" });
  win.timers.advance(5000);
  await settle();
  assert.equal(win.fetchCalls.length, 0);
});

test("Stop waiting stops polling and Check again resumes it", async () => {
  const win = installDom("http://localhost/phenomena");
  win.fetchHandler = () => job();
  const el = await start(win);
  win.timers.advance(1000);
  await settle();
  el.querySelector("[data-slow-stop]").dispatchEvent(new Event("click"));
  assert.match(text(el, "status"), /^Stopped waiting\./);
  assert.equal(el.querySelector("[data-slow-resume]").hidden, false);
  win.timers.advance(5000);
  await settle();
  assert.equal(win.fetchCalls.length, 1, "no more polls once stopped");
  el.querySelector("[data-slow-resume]").dispatchEvent(new Event("click"));
  await settle();
  assert.equal(win.fetchCalls.length, 2);
  assert.equal(el.querySelector("[data-slow-stop]").hidden, false);
});

test("an answer that arrives after Stop waiting does not undo the stop (TEST.235)", async () => {
  const win = installDom("http://localhost/phenomena");
  let answer;
  win.fetchHandler = () => new Promise((resolve) => { answer = resolve; });
  const el = await start(win);
  win.timers.advance(1000);
  await settle();
  assert.equal(win.fetchCalls.length, 1, "the poll is on its way");
  el.querySelector("[data-slow-stop]").dispatchEvent(new Event("click"));
  answer(job());
  await settle();
  assert.match(text(el, "status"), /^Stopped waiting\./);
  assert.equal(el.querySelector("[data-slow-resume]").hidden, false);
  assert.equal(el.querySelector("[data-slow-stop]").hidden, true);
  win.timers.advance(5000);
  await settle();
  assert.equal(win.fetchCalls.length, 1, "no more polls once stopped");
  win.fetchHandler = () => job({ done: 600 });
  el.querySelector("[data-slow-resume]").dispatchEvent(new Event("click"));
  await settle();
  assert.equal(text(el, "percent"), "60%", "Check again still works after a stop that raced a poll");
});

test("a finished job asks for the page at the same address, and a page waiting for another job carries on with that one", async () => {
  const win = installDom("http://localhost/search?q=sol");
  const next = h("section", { "data-slow-load": "", "data-job": "fedcba9876543210", "data-state": "queued",
    "data-status-url": "/api/page-jobs/fedcba9876543210" });
  globalThis.DOMParser = class { parseFromString() { return { querySelector: () => next }; } };
  win.fetchHandler = (url) => {
    if (url.startsWith("/api/page-jobs/0123")) return job({ state: "succeeded", done: 1000, eta_seconds: 0 });
    if (url.startsWith("/api/page-jobs/fedc")) return job({ id: "fedcba9876543210", done: 500 });
    return { ok: true, status: 202, json: async () => ({}), text: async () => "<html>waiting</html>" };
  };
  const el = await start(win);
  win.timers.advance(1000);
  await settle(30);
  const pages = win.fetchCalls.filter((call) => !call.url.startsWith("/api/"));
  assert.equal(pages.length, 1);
  assert.equal(pages[0].url, "http://localhost/search?q=sol");
  assert.equal(pages[0].options.headers["X-Slow-Retry"], undefined);
  assert.equal(el.dataset.job, "fedcba9876543210");
  assert.equal(el.dataset.statusUrl, "/api/page-jobs/fedcba9876543210");
  assert.equal(el.querySelector("[data-slow-bar]").value, 50);
});

test("a job that has expired (404) sends the visitor to the page itself", async () => {
  const win = installDom("http://localhost/phenomena");
  globalThis.DOMParser = class { parseFromString() { return { querySelector: () => null, getElementById: () => null }; } };
  win.fetchHandler = (url) => url.startsWith("/api/page-jobs/")
    ? jsonResponse({ error: "no such page job" }, 404)
    : { ok: true, status: 200, json: async () => ({}), text: async () => "<html></html>" };
  await start(win);
  win.timers.advance(1000);
  await settle(30);
  assert.equal(win.fetchCalls.filter((call) => !call.url.startsWith("/api/")).length, 1);
  assert.equal(win.location.reloads, 1, "a page it cannot swap into is simply reloaded");
});

test("Try again asks for the page with the retry header", async () => {
  const win = installDom("http://localhost/phenomena");
  const next = h("section", { "data-slow-load": "", "data-job": "aaaaaaaaaaaaaaaa", "data-state": "queued",
    "data-status-url": "/api/page-jobs/aaaaaaaaaaaaaaaa" });
  globalThis.DOMParser = class { parseFromString() { return { querySelector: () => next }; } };
  win.fetchHandler = (url) => url.startsWith("/api/page-jobs/")
    ? job({ id: "aaaaaaaaaaaaaaaa", state: "queued", done: 0 })
    : { ok: true, status: 202, json: async () => ({}), text: async () => "<html>waiting</html>" };
  const el = await start(win, { "data-state": "failed" });
  el.querySelector("[data-slow-retry]").dispatchEvent(new Event("click"));
  await settle(30);
  const pages = win.fetchCalls.filter((call) => !call.url.startsWith("/api/"));
  assert.equal(pages.length, 1);
  assert.equal(pages[0].options.headers["X-Slow-Retry"], "1");
  assert.equal(el.dataset.job, "aaaaaaaaaaaaaaaa");
});
