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

test("the page address is the path with the same pairs the no-script links use", () => {
  assert.equal(pageAddress("/phenomena", STATE, "name"), "/phenomena?sort=type&order=desc&type=nebula&type=rogue_planet");
  assert.equal(pageAddress("/phenomena", { sort: "name", descending: false, filters: {} }, "name"), "/phenomena");
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
