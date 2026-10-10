// tests/js/starfade.test.mjs -- static/starfade.js, the rank birth-radius
// fade that makes the Galaxy Map's stars come in with the zoom (MAP.153):
// the arithmetic on its own, then a synthetic star field laid out in tile
// levels the way the server lists them, zoomed in 9% at a time, to bound
// how much the picture changes in one step.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const S = await importPage("starfade.js");

const ROOT = 65536;
const FACTOR = 1.6;
const edgeOf = (level) => ROOT / 2 ** level;

test("a star starts to appear at its birth radius and is whole one octave of zoom closer", () => {
  const birth = S.birthRadius(edgeOf(10), 100);
  assert.equal(S.zoomShare(birth, birth), 0);
  assert.equal(S.zoomShare(birth, birth * 1.5), 0);
  assert.equal(S.zoomShare(birth, birth / 2), 1);
  assert.equal(S.zoomShare(birth, birth / 8), 1);
  assert.ok(Math.abs(S.zoomShare(birth, birth / Math.SQRT2) - 0.5) < 1e-12);
  assert.equal(S.zoomShare(0, 1), 0);
});

test("the birth radius follows the design's formula and falls with the rank", () => {
  // Level 12 (edge 16 pc) serves camera radii from 5 pc: R* = 5, times 2^W.
  assert.ok(Math.abs(S.birthRadius(16, 400) - 10) < 1e-9);
  assert.ok(Math.abs(S.birthRadius(16, 1) - 10 * 400 ** (1 / 3)) < 1e-9);
  assert.ok(S.birthRadius(16, 10) > S.birthRadius(16, 11));
  assert.ok(S.birthRadius(32, 400) > S.birthRadius(16, 400));
});

test("the level glide runs 0 to 1 across the octave a level serves", () => {
  const edge = edgeOf(8);
  assert.equal(S.levelGlide(edge, edge / FACTOR), 0);
  assert.equal(S.levelGlide(edge, edge / FACTOR * 2), 0);
  assert.ok(Math.abs(S.levelGlide(edge, edge / (2 * FACTOR)) - 1) < 1e-12);
  assert.equal(S.levelGlide(edge, edge / (4 * FACTOR)), 1);
});

test("at a level boundary a star shows what its parent's list gave it, from either side", () => {
  const edge = edgeOf(8);
  const parent = S.birthRadius(edge * 2, 120);
  const own = S.birthRadius(edge, 300);
  const fade = [own, parent, edge, 0];
  const boundary = edge / FACTOR;
  assert.ok(Math.abs(S.starZoomOpacity(fade, boundary) - S.zoomShare(parent, boundary)) < 1e-12);
  assert.ok(Math.abs(S.starZoomOpacity(fade, boundary / 2) - S.zoomShare(own, boundary / 2)) < 1e-12);
  // A star the parent did not list is not there at the zoomed-out end.
  assert.equal(S.starZoomOpacity([own, 0, edge, 0], boundary), 0);
});

test("a detail star uses its own rank alone, and a point is always whole", () => {
  assert.equal(S.starZoomOpacity([S.birthRadius(16, 50), 0, 16, 1], 1e6), 0);
  assert.equal(S.starZoomOpacity([S.birthRadius(16, 50), 0, 16, 1], 4), 1);
  assert.equal(S.starZoomOpacity(S.POINT_FADE, 6000), 1);
  assert.equal(S.starZoomOpacity(S.POINT_FADE, 4), 1);
});

test("tileRanks counts each list in the order the server sent it", () => {
  const ranks = S.tileRanks({ stars: [{ id: 7 }, { id: 3 }], generated: [{ id: 3 }, { id: 9 }] });
  assert.equal(ranks.get("b7"), 1);
  assert.equal(ranks.get("b3"), 2);
  assert.equal(ranks.get("s3"), 1);
  assert.equal(ranks.get("s9"), 2);
});

// --- A synthetic galaxy corner, listed the server's way -------------------

function mulberry32(seed) {
  return function () {
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Stars strewn around a target with a density that falls as r^-1.5 from
// it (a dense core), luminosity a power law; every tile level lists the
// most luminous of the stars in the cube of its edge around the target
// (400 bright at each, the finest listing every star within 8 pc).
function field(count, seed) {
  const random = mulberry32(seed);
  const stars = [];
  for (let i = 0; i < count; i++) {
    const r = 3 * Math.pow(1 + random() * 3000, 1) ** 0.9 / 10;
    const cos = 2 * random() - 1;
    const phi = 2 * Math.PI * random();
    const s = Math.sqrt(1 - cos * cos);
    stars.push({
      id: i, x: r * s * Math.cos(phi), y: r * s * Math.sin(phi), z: r * cos,
      luminosity: Math.pow(random(), 5) * 1000,
    });
  }
  return stars;
}

function levelLists(stars, cap) {
  const lists = [];
  for (let level = 0; level <= 12; level++) {
    const half = edgeOf(level) / 2;
    lists[level] = stars
      .filter((s) => Math.abs(s.x) < half && Math.abs(s.y) < half && Math.abs(s.z) < half)
      .sort((p, q) => q.luminosity - p.luminosity)
      .slice(0, level === 12 ? 4000 : cap);
  }
  return lists;
}

// What the page draws at camera radius `radius`: [fade] for every star.
function drawn(lists, radius) {
  const view = FACTOR * radius;
  const level = Math.max(0, Math.min(12, Math.floor(Math.log2(ROOT / view))));
  const rank = lists.map((list) => new Map(list.map((s, i) => [s.id, i + 1])));
  const out = new Map();
  for (const star of lists[level]) {
    const own = S.birthRadius(edgeOf(level), rank[level].get(star.id), FACTOR);
    const parentRank = level > 0 ? rank[level - 1].get(star.id) : 0;
    const parent = parentRank ? S.birthRadius(edgeOf(level - 1), parentRank, FACTOR) : 0;
    out.set(star.id, [own, parent, edgeOf(level), 0]);
  }
  if (level > 0) {
    for (const star of lists[level - 1]) {
      if (!out.has(star.id)) {
        out.set(star.id, [0, S.birthRadius(edgeOf(level - 1), rank[level - 1].get(star.id), FACTOR), edgeOf(level), 0]);
      }
    }
  }
  if (level < 12 && view <= 160) {
    for (const star of lists[12]) {
      if (Math.hypot(star.x, star.y, star.z) <= Math.min(view, 8) && !out.has(star.id)) {
        out.set(star.id, [S.birthRadius(edgeOf(12), rank[12].get(star.id), FACTOR), 0, edgeOf(12), 1]);
      }
    }
  }
  return out;
}

test("zooming in 9% at a time never adds or drops a burst of stars", () => {
  const stars = field(20000, 7);
  const lists = levelLists(stars, 400);
  let radius = 6000;
  let before = drawn(lists, radius);
  let worstRise = 0;
  let worstFall = 0;
  let steps = 0;
  while (radius / 2 ** (1 / 8) > 4) {
    radius /= 2 ** (1 / 8);
    const after = drawn(lists, radius);
    const ids = new Set([...before.keys(), ...after.keys()]);
    let rise = 0;
    let fall = 0;
    let total = 0;
    for (const id of ids) {
      const a = before.has(id) ? S.starZoomOpacity(before.get(id), radius * 2 ** (1 / 8), FACTOR) : 0;
      const b = after.has(id) ? S.starZoomOpacity(after.get(id), radius, FACTOR) : 0;
      rise += Math.max(0, b - a);
      fall += Math.max(0, a - b);
      total += a;
    }
    if (total >= 20) {
      worstRise = Math.max(worstRise, rise / total);
      worstFall = Math.max(worstFall, fall / total);
      steps++;
    }
    before = after;
  }
  assert.ok(steps > 20, `only ${steps} steps had 20 stars on screen`);
  assert.ok(worstRise < 1.5, `one 9% zoom step raised the opacity on screen by ${(100 * worstRise).toFixed(0)}%`);
  assert.ok(worstFall < 0.25, `one 9% zoom step dropped the opacity on screen by ${(100 * worstFall).toFixed(0)}%`);
});

test("the picture depends on the camera radius alone, in or out", () => {
  const lists = levelLists(field(8000, 11), 400);
  const there = (radius) => S.opacitySum(Array.from(drawn(lists, radius).values()), radius, FACTOR);
  const inward = [];
  for (let r = 3000; r > 8; r /= 1.1) inward.push([r, there(r)]);
  for (const [r, sum] of inward.reverse()) {
    assert.ok(Math.abs(there(r) - sum) < 1e-9);
  }
});

test("the same zoom without the fade (every listed star whole) is far burstier", () => {
  const lists = levelLists(field(20000, 7), 400);
  let radius = 6000;
  let before = drawn(lists, radius);
  let worstFade = 0;
  let worstHard = 0;
  while (radius / 2 ** (1 / 8) > 4) {
    const farther = radius;
    radius /= 2 ** (1 / 8);
    const after = drawn(lists, radius);
    const hardBefore = new Set(before.keys());
    const hardAfter = new Set(after.keys());
    let added = 0;
    for (const id of hardAfter) if (!hardBefore.has(id)) added++;
    let risen = 0;
    let total = 0;
    for (const [id, fade] of after) {
      const a = before.has(id) ? S.starZoomOpacity(before.get(id), farther, FACTOR) : 0;
      risen += Math.max(0, S.starZoomOpacity(fade, radius, FACTOR) - a);
    }
    for (const fade of before.values()) total += S.starZoomOpacity(fade, farther, FACTOR);
    if (hardBefore.size >= 20) worstHard = Math.max(worstHard, added / hardBefore.size);
    if (total >= 20) worstFade = Math.max(worstFade, risen / total);
    before = after;
  }
  assert.ok(worstHard > 2 * worstFade, `without the fade ${(100 * worstHard).toFixed(0)}% vs with it ${(100 * worstFade).toFixed(0)}%`);
});
