// tests/js/levelprogress.test.mjs -- static/levelprogress.js and the progressive reader in static/tilepack.js (MAP.364): a bar for each level of
// nodes while they load, none for a fast or cached load, plain words when a load stalls.

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, importPage } from "./fakedom.mjs";

installDom();
const P = await importPage("levelprogress.js");
const T = await importPage("tilepack.js");

function clock() {
  let now = 0;
  let id = 0;
  const pending = new Map();
  return {
    timers: {
      set(fn, ms) { pending.set(++id, { fn, at: now + ms }); return id; },
      clear(handle) { pending.delete(handle); },
    },
    advance(ms) {
      now += ms;
      for (const [handle, item] of Array.from(pending.entries())) {
        if (item.at <= now) { pending.delete(handle); item.fn(); }
      }
    },
  };
}

function host() {
  const el = document.createElement("div");
  document.body.appendChild(el);
  return el;
}

const KEYS = ["243.1.0.0", "243.1.0.1", "27.4.0.0", "1.9.0.0"];

test("a node's level is its first number, and a level is named in plain words", () => {
  assert.equal(P.levelOf("243.1.0.0"), 243);
  assert.equal(P.levelOf("1.9.0.0"), 1);
  assert.equal(P.levelOf("x"), null);
  assert.equal(P.levelName(1), "sectors");
  assert.equal(P.levelName(27), "blocks 27 shells deep");
});

test("no bar shows for a load that ends inside the delay, which is every load served from a cache", () => {
  const c = clock();
  const el = host();
  const p = P.createLevelProgress(el, { timers: c.timers });
  p.start(KEYS);
  c.advance(P.SHOW_DELAY_MS - 1);
  assert.equal(el.querySelector(".map-level-progress"), null);
  p.arrived(KEYS);
  p.finish();
  c.advance(P.SHOW_DELAY_MS * 3);
  assert.equal(el.querySelector(".map-level-progress"), null);
});

test("a slow load shows one progressbar for each level, named and counting nodes received over nodes requested", () => {
  const c = clock();
  const el = host();
  const p = P.createLevelProgress(el, { timers: c.timers });
  p.start(KEYS);
  c.advance(P.SHOW_DELAY_MS);
  const bars = el.querySelectorAll('[role="progressbar"]');
  assert.equal(bars.length, 3);
  const bar = bars.find((b) => b.getAttribute("aria-label") === "Loading blocks 243 shells deep");
  assert.ok(bar);
  assert.equal(bar.getAttribute("aria-valuemax"), "2");
  assert.equal(bar.getAttribute("aria-valuenow"), "0");
  p.arrived(["243.1.0.0"]);
  assert.equal(bar.getAttribute("aria-valuenow"), "1");
  assert.match(bar.getAttribute("aria-valuetext"), /1 of 2/);
  p.arrived(["243.1.0.1", "243.1.0.1"]);
  assert.equal(bar.getAttribute("aria-valuenow"), "2");
  p.arrived(["27.4.0.0", "1.9.0.0"]);
  p.finish();
  c.advance(P.LINGER_MS);
  assert.equal(el.querySelector(".map-level-progress"), null, "the bars are gone when the levels are complete");
});

test("a failed load says so in plain words and when it tries again, instead of a bar stuck part way", () => {
  const c = clock();
  const el = host();
  const p = P.createLevelProgress(el, { timers: c.timers });
  p.start(KEYS);
  c.advance(P.SHOW_DELAY_MS);
  p.arrived(["243.1.0.0", "243.1.0.1"]);
  p.stalled(4);
  assert.equal(el.querySelectorAll('[role="progressbar"]').length, 0);
  const lines = el.querySelectorAll('[role="status"]');
  assert.equal(lines.length, 2, "the finished level has no line");
  assert.match(lines[0].textContent, /^Could not load all of the blocks 27 shells deep\. Trying again in 4 seconds\.$/);
  p.start(["27.4.0.0"]);
  c.advance(P.SHOW_DELAY_MS);
  assert.equal(el.querySelectorAll('[role="status"]').length, 0, "a new request clears the lines");
});

test("an incomplete node keeps its level's line after the load ends", () => {
  const c = clock();
  const el = host();
  const p = P.createLevelProgress(el, { timers: c.timers });
  p.start(KEYS);
  c.advance(P.SHOW_DELAY_MS);
  p.arrived(KEYS);
  p.incomplete(["1.9.0.0"]);
  p.finish();
  c.advance(P.LINGER_MS * 2);
  const lines = el.querySelectorAll('[role="status"]');
  assert.equal(lines.length, 1);
  assert.match(lines[0].textContent, /still being worked out/);
});

function packed(listed, sizes) {
  const header = new TextEncoder().encode(JSON.stringify({ nodes: listed }));
  const body = new Uint8Array(8 + header.length + sizes.reduce((a, b) => a + b, 0));
  body.set(new TextEncoder().encode("PGTB"), 0);
  new DataView(body.buffer).setUint32(4, header.length, true);
  body.set(header, 8);
  return body;
}

test("the reader reports each node as its bytes arrive, and not before the header is whole", () => {
  const body = packed([["27.4.0.0", 10], ["1.9.0.0", 6]], [10, 6]);
  const headerEnd = body.length - 16;
  assert.deepEqual(T.nodesSoFar(body.subarray(0, 5), null), { listed: null, keys: [] });
  assert.deepEqual(T.nodesSoFar(body.subarray(0, headerEnd - 1), null).keys, []);
  const first = T.nodesSoFar(body.subarray(0, headerEnd + 9), null);
  assert.deepEqual(first.keys, []);
  const second = T.nodesSoFar(body.subarray(0, headerEnd + 10), first.listed);
  assert.deepEqual(second.keys, ["27.4.0.0"]);
  assert.deepEqual(T.nodesSoFar(body, second.listed).keys, ["27.4.0.0", "1.9.0.0"]);
  assert.throws(() => T.nodesSoFar(new Uint8Array(20), null), /not a packed node response/);
});

test("readNodeBody reads a streamed response whole and tells of each node once", async () => {
  const body = packed([["27.4.0.0", 10], ["1.9.0.0", 6]], [10, 6]);
  const chunks = [body.slice(0, 30), body.slice(30, body.length - 7), body.slice(body.length - 7)];
  let i = 0;
  const response = { body: { getReader: () => ({ read: async () => (i < chunks.length ? { done: false, value: chunks[i++] } : { done: true }) }) } };
  const told = [];
  const whole = await T.readNodeBody(response, (keys) => told.push(keys));
  assert.equal(whole.byteLength, body.length);
  assert.deepEqual(told.flat(), ["27.4.0.0", "1.9.0.0"]);
  const plain = await T.readNodeBody({ arrayBuffer: async () => body.buffer }, (keys) => told.push(keys));
  assert.equal(plain.byteLength, body.length);
});
