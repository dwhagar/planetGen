// tests/js/starbuffer.test.mjs -- static/starbuffer.js (MAP.267): the relative-to-eye pairs keep a star under a quarter
// pixel down to 1e-6 pc at any distance from the origin (float32 emulated with Math.fround, as the GPU rounds), and the
// slot store keeps the star buffer dense as nodes come and go.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const B = await importPage("starbuffer.js");

// A small deterministic generator (mulberry32), so the placements are the same every run.
function random(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const FOCAL_PX = 800 / (2 * Math.tan((25 * Math.PI) / 180));

// The worst error in pixels of a vertex `distance` from the camera, over `trials` placements near `origin` (one axis
// each, the camera's position anywhere in a 2 pc cube there), for the pair arithmetic and for a single float32.
function worst(distance, origin, trials, rand) {
  let pair = 0;
  let single = 0;
  for (let i = 0; i < trials; i += 1) {
    const camera = origin + (rand() - 0.5) * 2;
    const star = camera + distance * (rand() > 0.5 ? 1 : -1) * (0.5 + rand() / 2);
    const truth = star - camera;
    const [hi, lo] = B.splitPair(star);
    const [camHi, camLo] = B.splitPair(camera);
    const got = B.relativeOffset(hi, lo, camHi, camLo);
    pair = Math.max(pair, (Math.abs(got - truth) / Math.abs(truth)) * Math.abs(truth) / distance * FOCAL_PX);
    single = Math.max(single, (Math.abs(B.singleOffset(star, camera) - truth) / distance) * FOCAL_PX);
  }
  return { pair, single };
}

for (const origin of [8200, 262000]) {
  test(`the pair holds a star under a quarter pixel from 1e5 down to 1e-6 pc, ${origin} pc from the origin`, () => {
    const rand = random(origin);
    for (const exponent of [5, 3, 1, 0, -1, -3, -5, -6]) {
      const found = worst(Math.pow(10, exponent), origin, 400, rand);
      assert.ok(found.pair < 0.25, `d=1e${exponent}: ${found.pair} px`);
    }
  });
}

test("one float32 would not: a star a parsec away is a dozen pixels off 262 kpc out", () => {
  const found = worst(1, 262000, 400, random(7));
  assert.ok(found.single > 4, `${found.single} px`);
  assert.ok(found.pair < 0.25);
});

test("a pair adds back to the value it was made from, to a float32 of what is left", () => {
  for (const value of [0, 1, -3.5, 8200.123456789, 262144.987654321, -1e-9]) {
    const [hi, lo] = B.splitPair(value);
    assert.equal(hi, Math.fround(hi));
    assert.equal(lo, Math.fround(lo));
    assert.ok(Math.abs(hi + lo - value) <= Math.abs(value) * 1e-14 + 1e-30, `${value}`);
  }
});

test("runsOf merges adjacent slots", () => {
  assert.deepEqual(B.runsOf([]), []);
  assert.deepEqual(B.runsOf([2, 3, 4, 7, 9, 10]), [[2, 3], [7, 1], [9, 2]]);
});

// The slot store against a plain model: every operation leaves the buffer dense with each chunk's stars where it says.
function check(store, model) {
  assert.equal(store.n, model.length);
  assert.equal(store.owners.length, store.n);
  store.keys().forEach((key) => {
    const slots = store.slotsOf(key);
    slots.forEach((slot) => assert.equal(store.owners[slot], key));
    assert.ok(slots.every((slot) => slot < store.n));
  });
  const counted = new Map();
  store.owners.forEach((key) => counted.set(key, (counted.get(key) || 0) + 1));
  counted.forEach((count, key) => assert.equal(store.slotsOf(key).length, count));
}

test("the slot store stays dense through adds and removes, and the moves carry the right stars", () => {
  const store = new B.SlotStore();
  const rand = random(3);
  // `model` is the buffer: each slot holds a label "key:index"; the moves are applied to it as the page would.
  let model = [];
  const held = new Map();
  let next = 0;
  for (let step = 0; step < 400; step += 1) {
    if (held.size === 0 || rand() < 0.55) {
      const key = `n${next++}`;
      const count = Math.floor(rand() * 9);
      const start = store.add(key, count, step);
      assert.equal(start, model.length);
      for (let i = 0; i < count; i += 1) model.push(`${key}:${i}`);
      held.set(key, count);
    } else {
      const keys = Array.from(held.keys());
      const key = keys[Math.floor(rand() * keys.length)];
      const moves = store.remove(key);
      moves.forEach((move) => { model[move.to] = model[move.from]; });
      model.length = store.n;
      held.delete(key);
      assert.ok(!model.some((label) => label.startsWith(key + ":")));
    }
    check(store, model);
    // Every held chunk's stars are all there, once each.
    held.forEach((count, key) => {
      const labels = store.slotsOf(key).map((slot) => model[slot]).sort();
      assert.deepEqual(labels, Array.from({ length: count }, (_, i) => `${key}:${i}`).sort());
    });
  }
});

test("a chunk is added once, and removing an unknown one does nothing", () => {
  const store = new B.SlotStore();
  store.add("a", 2, 1);
  assert.throws(() => store.add("a", 1, 2));
  assert.deepEqual(store.remove("zzz"), []);
  assert.equal(store.tagOf("a"), 1);
  assert.equal(store.tagOf("b"), undefined);
});
