// tests/js/savedlist.test.mjs -- static/savedlist.js, the saved-list store (UX.129): the same behaviour over the
// localStorage backend and over a slow, asynchronous fake backend, a second list kind with no new storage code,
// the instance-stamp hook, and bookmarks.js on top of it.

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/galaxy");
const S = await import(new URL("savedlist.js", STATIC_URL).href);

const slow = () => new Promise((resolve) => setTimeout(resolve, 5));

// A backend that keeps its list in memory and answers a few milliseconds late, like an account over the network.
function fakeAsyncBackend(initial) {
  const state = { stored: initial || [], listeners: [], saves: 0 };
  return {
    state,
    available: () => true,
    load: async () => { await slow(); return JSON.parse(JSON.stringify(state.stored)); },
    save: async (entries) => { await slow(); state.saves += 1; state.stored = JSON.parse(JSON.stringify(entries)); return true; },
    subscribe: (fn) => { state.listeners.push(fn); },
  };
}

function localStore(options) {
  const db = { name: "galaxy_a" };
  const store = S.createStore(Object.assign({
    backend: S.localBackend({ prefix: "planetgen.test.", suffix: () => db.name }),
    keyOf: (e) => (e && typeof e.id === "string" && typeof e.name === "string" ? e.id : null),
    max: 3,
  }, options || {}));
  return { store, db };
}

function asyncStore(backend, options) {
  return S.createStore(Object.assign({
    backend,
    keyOf: (e) => (e && typeof e.id === "string" && typeof e.name === "string" ? e.id : null),
    max: 3,
  }, options || {}));
}

// The same checks over either backend.
async function behaves(store) {
  await store.ready;
  assert.deepEqual(store.list(), []);
  const changes = [];
  store.onChange(() => changes.push(store.list().length));

  assert.deepEqual(await store.add({ id: "a", name: "Alpha" }), { ok: true, entry: { id: "a", name: "Alpha" } });
  assert.equal(store.find("a").name, "Alpha");
  assert.equal(store.find("zzz"), null);
  assert.deepEqual((await store.add({ id: "a", name: "Other" })).entry, { id: "a", name: "Alpha" }, "already saved is ok, as it was");
  assert.equal((await store.add({ id: "b", name: "Beta" })).ok, true);
  assert.deepEqual(store.list().map((e) => e.id), ["a", "b"], "oldest first");

  assert.equal((await store.add({ nope: 1 })).ok, false);
  assert.equal(await store.rename("b", "  Bravo  "), true);
  assert.equal(store.find("b").name, "Bravo");
  assert.equal(await store.rename("b", "   "), false);
  assert.equal(await store.rename("missing", "x"), false);

  assert.equal((await store.add({ id: "c", name: "C" })).ok, true);
  const full = await store.add({ id: "d", name: "D" });
  assert.equal(full.ok, false);
  assert.match(full.reason, /3/);

  assert.equal(await store.remove("a"), true);
  assert.equal(await store.remove("a"), false);
  assert.deepEqual(store.list().map((e) => e.id), ["b", "c"]);
  assert.ok(changes.length >= 5, "onChange after every change: " + changes);
}

test("the local backend behaves", async () => {
  window.localStorage.clear();
  const { store } = localStore();
  await behaves(store);
  assert.deepEqual(JSON.parse(window.localStorage.getItem("planetgen.test.galaxy_a")).map((e) => e.id), ["b", "c"]);
});

test("a fake asynchronous backend passes the same checks", async () => {
  const backend = fakeAsyncBackend();
  await behaves(asyncStore(backend));
  assert.deepEqual(backend.state.stored.map((e) => e.id), ["b", "c"]);
});

test("a store with no database named keeps nothing and says so", async () => {
  window.localStorage.clear();
  const { store, db } = localStore({ words: { blocked: "blocked" } });
  db.name = "";
  await store.ready;
  assert.equal(store.available(), false);
  const result = await store.add({ id: "a", name: "A" });
  assert.deepEqual(result, { ok: false, reason: "blocked" });
  assert.deepEqual(store.list(), []);
});

test("blocked storage gives no entries and a refusal", async () => {
  const backend = S.localBackend({ prefix: "p.", suffix: () => "db", storage: () => { throw new Error("blocked"); } });
  const store = S.createStore({ backend, keyOf: (e) => (e && e.id) || null, max: 3 });
  await store.ready;
  assert.deepEqual(store.list(), []);
  assert.equal((await store.add({ id: "a", name: "A" })).ok, false);
});

test("a second kind of list needs no new storage code", async () => {
  window.localStorage.clear();
  const courses = S.createStore({
    backend: S.localBackend({ prefix: "planetgen.courses.", suffix: () => "galaxy_a" }),
    keyOf: (e) => (e && typeof e.from === "string" && typeof e.to === "string" ? e.from + ">" + e.to : null),
    max: 50,
  });
  const marks = localStore().store;
  await Promise.all([courses.ready, marks.ready]);
  await courses.add({ from: "system:1", to: "system:2", name: "Home run" });
  await marks.add({ id: "a", name: "A" });
  assert.equal(courses.list().length, 1);
  assert.equal(marks.list().length, 1);
  assert.ok(window.localStorage.getItem("planetgen.courses.galaxy_a"));
  assert.equal(courses.find("system:1>system:2").name, "Home run");
});

test("the stamp hook records what it gives with each new entry", async () => {
  const backend = fakeAsyncBackend();
  const store = asyncStore(backend, { stamp: () => ({ galaxy: "g1", epoch: 7 }) });
  await store.ready;
  const result = await store.add({ id: "a", name: "A" });
  assert.deepEqual(result.entry.stamp, { galaxy: "g1", epoch: 7 });
  assert.deepEqual(backend.state.stored[0].stamp, { galaxy: "g1", epoch: 7 });
});

test("writes made at once are saved in order and none is lost", async () => {
  const backend = fakeAsyncBackend();
  const store = asyncStore(backend, { max: 10 });
  await store.ready;
  await Promise.all([store.add({ id: "a", name: "A" }), store.add({ id: "b", name: "B" }), store.add({ id: "c", name: "C" })]);
  assert.deepEqual(store.list().map((e) => e.id), ["a", "b", "c"]);
  assert.deepEqual(backend.state.stored.map((e) => e.id), ["a", "b", "c"]);
});

test("a change made elsewhere is picked up, and not overwritten by the next write", async () => {
  const backend = fakeAsyncBackend([{ id: "a", name: "A" }]);
  const store = asyncStore(backend, { max: 10 });
  await store.ready;
  let changed = 0;
  store.onChange(() => { changed += 1; });
  backend.state.stored = [{ id: "a", name: "A" }, { id: "x", name: "Elsewhere" }];
  backend.state.listeners.forEach((fn) => fn());
  await slow(); await slow();
  assert.deepEqual(store.list().map((e) => e.id), ["a", "x"]);
  assert.equal(changed, 1);
  // A change that is not yet announced is still read before the write.
  backend.state.stored = backend.state.stored.concat([{ id: "y", name: "Quiet" }]);
  await store.add({ id: "b", name: "B" });
  assert.deepEqual(backend.state.stored.map((e) => e.id), ["a", "x", "y", "b"]);
});

test("bookmarks.js keeps the format and key it always had", async () => {
  installDom("http://localhost/galaxy");
  const body = document.body;
  const holder = document.createElement("div");
  holder.setAttribute("data-bookmark-db", "galaxy_b");
  body.appendChild(holder);
  window.localStorage.setItem("planetgen.bookmarks.galaxy_b", JSON.stringify([
    { name: "Old", kind: "system", value: "system:00001000000000-4000000-000", created: "2026-01-01T00:00:00.000Z" },
    { nonsense: true },
  ]));
  const B = await import(new URL("bookmarks.js", STATIC_URL).href);
  await B.ready;
  assert.deepEqual(B.list().map((e) => e.name), ["Old"], "a stored entry is read, a malformed one dropped");
  assert.equal(B.find("system", "system:00001000000000-4000000-000").name, "Old");
  const added = await B.add({ name: "  New  ", kind: "sector", value: "FE81000A2B", url: "/sector/7" });
  assert.equal(added.ok, true);
  assert.equal(added.entry.name, "New");
  assert.equal(added.entry.url, "/sector/7");
  assert.deepEqual(Object.keys(added.entry).sort(), ["created", "kind", "name", "url", "value"]);
  assert.equal(await B.rename("sector", "FE81000A2B", "Newer"), true);
  const stored = JSON.parse(window.localStorage.getItem("planetgen.bookmarks.galaxy_b"));
  assert.deepEqual(stored.map((e) => e.name), ["Old", "Newer"]);
  assert.equal(await B.remove("system", "system:00001000000000-4000000-000"), true);
  assert.deepEqual(B.list().map((e) => e.name), ["Newer"]);
  assert.equal((await B.add({ name: "", kind: "", value: "" })).reason, "That can't be bookmarked.");
});
