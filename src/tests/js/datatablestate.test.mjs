// tests/js/datatablestate.test.mjs -- static/datatablestate.js, the data
// table's logic (UX.41): the query a state turns into, the pages the
// visible rows need and the store of fetched pages.

import { test } from "node:test";
import assert from "node:assert/strict";

import { PAGE_SIZE, createPageStore, dataQuery, pageAddress, pagesFor, stateParams } from "../../html/static/datatablestate.js";

const STATE = { sort: "type", descending: true, filters: { type: ["nebula", "rogue_planet"], descriptor: [] } };

test("a state becomes sort, order and one pair per chosen filter value", () => {
  assert.deepEqual(stateParams(STATE, "name"), [["sort", "type"], ["order", "desc"], ["type", "nebula"], ["type", "rogue_planet"]]);
});

test("the default sort, ascending, is left out", () => {
  assert.deepEqual(stateParams({ sort: "name", descending: false, filters: {} }, "name"), []);
  assert.deepEqual(stateParams({ sort: "name", descending: true, filters: {} }, "name"), [["sort", "name"], ["order", "desc"]]);
});

test("the data query always names the sort and carries the window", () => {
  assert.equal(dataQuery({ sort: "name", descending: false, filters: {} }, 100, 50, true), "sort=name&offset=100&limit=50&facets=1");
  assert.equal(dataQuery(STATE, 0, 50, false), "sort=type&order=desc&type=nebula&type=rogue_planet&offset=0&limit=50");
});

test("values are percent-encoded", () => {
  const state = { sort: "name", descending: false, filters: { descriptor: ["a b&c"] } };
  assert.equal(dataQuery(state, 0, 50, false), "sort=name&descriptor=a%20b%26c&offset=0&limit=50");
});

test("the jump-by-value start rides in the data query as from, and in a prefixed address as the table's own from (PERF.75)", () => {
  const state = { sort: "name", descending: false, filters: {}, start: "Be a&b" };
  assert.equal(dataQuery(state, 0, 50, false), "sort=name&from=Be%20a%26b&offset=0&limit=50");
  assert.deepEqual(stateParams(state, "name", "systems_"), [["systems_from", "Be a&b"]]);
  assert.deepEqual(stateParams({ ...state, start: "" }, "name", "systems_"), []);
  assert.equal(pageAddress("/systems", "?x=1&systems_from=old", state, "name", "systems_", ["systems_sort", "systems_order", "systems_page", "systems_from"]),
    "/systems?x=1&systems_from=Be%20a%26b");
});

const OWNED = ["sort", "order", "page", "type", "descriptor"];

test("the page address is the path with the same pairs the no-script links use", () => {
  assert.equal(pageAddress("/phenomena", "", STATE, "name", "", OWNED), "/phenomena?sort=type&order=desc&type=nebula&type=rogue_planet");
  assert.equal(pageAddress("/phenomena", "?page=3&sort=x", { sort: "name", descending: false, filters: {} }, "name", "", OWNED), "/phenomena");
});

test("a table's prefix names its sort and order, and the page's other parameters stay", () => {
  const state = { sort: "systems", descending: true, filters: {} };
  assert.deepEqual(stateParams(state, "distance", "sectors_"), [["sectors_sort", "systems"], ["sectors_order", "desc"]]);
  const owned = ["sectors_sort", "sectors_order", "sectors_page"];
  assert.equal(
    pageAddress("/", "?standalone_sort=binary&sectors_page=4&sectors_sort=name", state, "distance", "sectors_", owned),
    "/?standalone_sort=binary&sectors_sort=systems&sectors_order=desc",
  );
});

test("rows map to the zero-based pages that hold them", () => {
  assert.deepEqual(pagesFor(0, 49, 50), [0]);
  assert.deepEqual(pagesFor(48, 52, 50), [0, 1]);
  assert.deepEqual(pagesFor(120, 260, 50), [2, 3, 4, 5]);
  assert.equal(PAGE_SIZE, 50);
});

function fakeServer(total) {
  const asked = [];
  const fetchPage = (page) => {
    asked.push(page);
    const rows = [];
    for (let i = page * 2; i < Math.min(total, page * 2 + 2); i += 1) {
      rows.push({ text: "row " + i });
    }
    return Promise.resolve({ rows, total, facets: page === 0 ? { type: [] } : null });
  };
  return { asked, fetchPage };
}

test("the store fetches a page once and answers rows from it", async () => {
  const server = fakeServer(5);
  const store = createPageStore(2, server.fetchPage);
  assert.equal(store.row(3), undefined);
  assert.equal(await store.ensure(1), true);
  assert.equal(await store.ensure(1), true);
  assert.deepEqual(server.asked, [1]);
  assert.deepEqual(store.row(3), { text: "row 3" });
  assert.equal(store.total, 5);
  assert.equal(store.has(1), true);
  assert.equal(store.has(0), false);
});

test("two callers asking for one page share its request", async () => {
  const server = fakeServer(5);
  const store = createPageStore(2, server.fetchPage);
  await Promise.all([store.ensure(0), store.ensure(0)]);
  assert.deepEqual(server.asked, [0]);
  assert.deepEqual(store.facets, { type: [] });
});

test("a reset forgets the pages and drops answers still on their way", async () => {
  let release;
  const slow = new Promise((resolve) => { release = resolve; });
  const store = createPageStore(2, () => slow);
  const stale = store.ensure(0);
  store.reset();
  release({ rows: [{ text: "old" }], total: 9 });
  assert.equal(await stale, false);
  assert.equal(store.row(0), undefined);
  assert.equal(store.total, 0);
});

test("a failed page is not kept, so the next call asks again", async () => {
  let calls = 0;
  const store = createPageStore(2, () => (++calls === 1 ? Promise.reject(new Error("down")) : Promise.resolve({ rows: [{ text: "x" }], total: 1 })));
  assert.equal(await store.ensure(0), false);
  assert.equal(await store.ensure(0), true);
  assert.equal(calls, 2);
});

import { busyPause, fetchWhileBusy } from "../../html/static/datatablestate.js";

test("the pause after a busy answer grows 2 s, 4 s ... and stops at 30 s", () => {
  assert.deepEqual([1, 2, 3, 4, 5, 6, 20].map(busyPause), [2000, 4000, 8000, 16000, 30000, 30000, 30000]);
});

function answer(status, body) {
  return { status, ok: status >= 200 && status < 300, json: async () => body };
}

test("a busy answer is asked again after the pause and the busy line shows meanwhile", async () => {
  const answers = [answer(504), answer(504), answer(200, { rows: [] })];
  const pauses = [];
  const busy = [];
  const body = await fetchWhileBusy(async () => answers.shift(), () => true, (state) => busy.push(state),
    async (ms) => { pauses.push(ms); });
  assert.deepEqual(body, { rows: [] });
  assert.deepEqual(pauses, [2000, 4000]);
  assert.deepEqual(busy, [true, true, false]);
});

test("a page that is no longer wanted stops asking", async () => {
  let wanted = true;
  const busy = [];
  await assert.rejects(fetchWhileBusy(async () => answer(504), () => wanted, (state) => busy.push(state),
    async () => { wanted = false; }), /no longer wanted/);
  assert.deepEqual(busy, [true, false]);
});

test("any other failing status throws and clears the busy line", async () => {
  const busy = [];
  await assert.rejects(fetchWhileBusy(async () => answer(502), () => true, (state) => busy.push(state), async () => {}),
    /table 502/);
  assert.deepEqual(busy, [false]);
});

test("the store knows an answer for an older generation is stale", () => {
  const store = createPageStore(50, async () => ({ rows: [], total: 0 }));
  assert.equal(store.stale(0), false);
  store.reset();
  assert.equal(store.stale(0), true);
});

test("the store keeps the note of the page it fetched and forgets it on reset (PERF.135)", async () => {
  const store = createPageStore(2, () => Promise.resolve({ rows: [{ text: "x" }], total: 1, note: "Nothing from “9” on; showing the last 1." }));
  assert.equal(store.note, "");
  await store.ensure(0);
  assert.equal(store.note, "Nothing from “9” on; showing the last 1.");
  store.reset();
  assert.equal(store.note, "");
  const plain = createPageStore(2, () => Promise.resolve({ rows: [{ text: "x" }], total: 1 }));
  await plain.ensure(0);
  assert.equal(plain.note, "");
});

import { JOB_POLL_MS, percentOf } from "../../html/static/datatablestate.js";

const pending = (state, extra = {}) => answer(202, { pending: true, percent: null, status_url: "/api/page-jobs/j1", job: { id: "j1", state, ...extra } });

function jobsFor(statuses, log) {
  return {
    shown: [], stopped: () => false, retry: false,
    show(job) { log.push(job && job.state); },
    poll: async (url) => { log.push(url); return statuses.length > 1 ? statuses.shift() : statuses[0]; },
  };
}

test("a 202 waits on the job, polling its status address, and asks for the rows again when it has succeeded", async () => {
  const log = [];
  const sleeps = [];
  const answers = [pending("queued"), answer(200, { rows: [1], total: 1 })];
  const jobs = jobsFor([answer(200, { id: "j1", state: "running", done: 50, total: 100 }), answer(200, { id: "j1", state: "succeeded" })], log);
  const body = await fetchWhileBusy(async () => answers.shift(), () => true, () => {}, async (ms) => { sleeps.push(ms); }, jobs);
  assert.deepEqual(body, { rows: [1], total: 1 });
  assert.deepEqual(log, ["queued", "/api/page-jobs/j1", "running", "/api/page-jobs/j1", "succeeded", null]);
  assert.deepEqual(sleeps, [JOB_POLL_MS, JOB_POLL_MS]);
});

test("a failed or cancelled job stays on the screen and the fetch fails", async () => {
  for (const state of ["failed", "cancelled"]) {
    const log = [];
    await assert.rejects(fetchWhileBusy(async () => pending(state, { error: "Disk full" }), () => true, () => {}, async () => {}, jobsFor([], log)),
      new RegExp(state));
    assert.deepEqual(log, [state]);
  }
});

test("Stop waiting ends the wait and a table no longer wanted clears the display", async () => {
  const log = [];
  const stopping = { ...jobsFor([], log), stopped: () => true };
  await assert.rejects(fetchWhileBusy(async () => pending("running"), () => true, () => {}, async () => {}, stopping), /stopped/);
  const logged = [];
  await assert.rejects(fetchWhileBusy(async () => pending("running"), () => false, () => {}, async () => {}, jobsFor([], logged)), /no longer wanted/);
  assert.deepEqual(logged, ["running", null]);
});

test("a job that has expired (404) is skipped: the rows are asked for again", async () => {
  const answers = [pending("running"), answer(200, { rows: [], total: 0 })];
  const body = await fetchWhileBusy(async () => answers.shift(), () => true, () => {}, async () => {}, jobsFor([answer(404, {})], []));
  assert.deepEqual(body, { rows: [], total: 0 });
});

test("the first ask carries the retry flag when the visitor pressed Try again", async () => {
  const asked = [];
  await fetchWhileBusy(async (retry) => { asked.push(retry); return answer(200, { rows: [] }); }, () => true, () => {}, async () => {}, { ...jobsFor([], []), retry: true });
  assert.deepEqual(asked, [true]);
});

test("the bar's value is a percent of rows done, 99 at most until it has succeeded, and unknown while queued", () => {
  assert.equal(percentOf({ state: "running", done: 300, total: 1000 }), 30);
  assert.equal(percentOf({ state: "running", done: 1000, total: 1000 }), 99);
  assert.equal(percentOf({ state: "queued" }), null);
  assert.equal(percentOf({ state: "succeeded" }), 100);
});
