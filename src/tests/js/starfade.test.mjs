// tests/js/starfade.test.mjs -- static/starfade.js: the tile fade that makes
// the Galaxy Map's stars come in and go out as the camera nears or leaves the
// tile that lists them (MAP.153, redone for tiles chosen by distance in
// MAP.152): the arithmetic on its own, then a synthetic star field listed
// tile by tile the server's way and approached from afar, to bound how much
// the picture changes in one small step.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const S = await importPage("starfade.js");
const C = await importPage("galaxytilecut.js");

const ROOT = 65536;
const edgeOf = (level) => ROOT / 2 ** level;

const parts = (key) => key.split("/").map(Number);
const FINEST_LEVEL = 11;

test("a tile gives way to its halves as the camera comes in, and is whole beyond", () => {
  const edge = edgeOf(8);
  const leaving = (u) => S.leavingShare(u * edge / S.CUT_EDGES, edge);
  assert.equal(leaving(0), 0);
  assert.equal(leaving(0.5), 0);
  assert.equal(leaving(1.3), 1);
  assert.equal(leaving(5), 1);
  assert.ok(leaving(1.0) > 0 && leaving(1.0) < 1);
});

test("a tile and its halves trade shares: they add up to one across the hand-over", () => {
  // A star at the centre of a tile, the camera approaching along an axis: the
  // tile and the half that holds the star, across the whole hand-over.
  const star = [8000.5, 0.5, 0.5];
  const family = C.fetchCutTiles([9000, 0, 0], [8000, 0, 0], {
    rootEdge: ROOT, finestLevel: FINEST_LEVEL, bounds: [[-12600, -12600, -960], [12600, 12600, 960]], halfAngleDeg: 180,
  }).filter((key) => {
    const box = C.tileBoxOf(key, ROOT);
    return [0, 1, 2].every((a) => star[a] >= box.lo[a] && star[a] < box.lo[a] + box.edge);
  });
  assert.ok(family.length >= 1);
  for (let distance = 5000; distance > 30; distance /= 1.04) {
    const camera = [star[0] + distance, star[1], star[2] + distance * 0.2];
    // The star's whole lineage, coarsest to finest.
    const lineage = [];
    for (let level = 0; level <= FINEST_LEVEL; level++) {
      const edge = edgeOf(level);
      lineage.push([level, ...[0, 1, 2].map((a) => Math.floor((star[a] + ROOT / 2) / edge))]);
    }
    const total = lineage.reduce((sum, p) => sum + S.tileShare(p, camera, ROOT, FINEST_LEVEL), 0);
    assert.ok(Math.abs(total - 1) < 1e-9, `at ${distance.toFixed(0)} pc the shares sum to ${total}`);
  }
});

test("the coarsest tile keeps its share far out and the finest of the cut close in", () => {
  assert.equal(S.tileShare([0, 0, 0, 0], [1e7, 0, 0], ROOT, FINEST_LEVEL), 1);
  // Inside a tile of the finest level (x 8192 to 8224, y and z 0 to 32).
  assert.equal(S.tileShare([FINEST_LEVEL, 1281, 1024, 1024], [8210, 5, 5], ROOT, FINEST_LEVEL), 1);
});

test("a detail tile is whole within 60 pc and gone at 120", () => {
  assert.equal(S.detailShare(10), 1);
  assert.equal(S.detailShare(60), 1);
  assert.equal(S.detailShare(120), 0);
  assert.equal(S.detailShare(500), 0);
  assert.ok(S.detailShare(90) > 0 && S.detailShare(90) < 1);
});

test("a star's share comes from its tile: box distance, a detail tile, no tile", () => {
  const key = parts("12/2048/2048/2048"); // the 16 pc cube at the origin
  assert.equal(S.boxDistance(key, [8, 8, 8], ROOT), 0);
  assert.ok(Math.abs(S.boxDistance(key, [-3, -4, 8], ROOT) - 5) < 1e-12);
  assert.ok(Math.abs(S.boxDistance(key, [19, 8, 20], ROOT) - 5) < 1e-12);
  assert.equal(S.starShare(null, [1e9, 0, 0], ROOT, FINEST_LEVEL), 1);
  assert.equal(S.starShare(S.POINT_FADE, [5, 5, 5], ROOT, FINEST_LEVEL), 1);
  assert.equal(S.starShare({ parts: key, detail: true }, [30, 0, 0], ROOT, FINEST_LEVEL), 1);
  assert.equal(S.starShare({ parts: key, detail: true }, [1000, 0, 0], ROOT, FINEST_LEVEL), 0);
});

test("a cloud sprite fades in between 1 and 4 px of radius", () => {
  assert.equal(S.cloudSizeShare(0.5), 0);
  assert.equal(S.cloudSizeShare(1), 0);
  assert.ok(Math.abs(S.cloudSizeShare(2.5) - 0.5) < 1e-9);
  assert.equal(S.cloudSizeShare(4), 1);
  assert.equal(S.cloudSizeShare(40), 1);
  let last = 0;
  for (let px = 1; px <= 4; px += 0.25) {
    const share = S.cloudSizeShare(px);
    assert.ok(share >= last);
    last = share;
  }
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

// Stars strewn through a 4 kpc cube round (8000, 0, 0), luminosity a power law.
function field(count, seed) {
  const random = mulberry32(seed);
  const stars = [];
  for (let i = 0; i < count; i++) {
    stars.push({
      id: i, x: 8000 + (random() - 0.5) * 4000, y: (random() - 0.5) * 4000, z: (random() - 0.5) * 600,
      luminosity: Math.pow(random(), 5) * 1000,
    });
  }
  return stars;
}

const BOUNDS = [[-12600, -12600, -960], [12600, 12600, 960]];
const CUT = { rootEdge: ROOT, finestLevel: 11, bounds: BOUNDS };

// What the page draws with the camera at `camera` looking at `look`: one
// entry [star id, fade] for every star every fetched tile lists (its `cap`
// most luminous).
function drawn(stars, cap, camera, look) {
  const out = [];
  for (const key of C.fetchCutTiles(camera, look, CUT)) {
    const box = C.tileBoxOf(key, ROOT);
    const inside = stars.filter((s) => s.x >= box.lo[0] && s.x < box.lo[0] + box.edge
      && s.y >= box.lo[1] && s.y < box.lo[1] + box.edge && s.z >= box.lo[2] && s.z < box.lo[2] + box.edge);
    inside.sort((p, q) => q.luminosity - p.luminosity).slice(0, cap).forEach((s) => {
      out.push([s.id, { parts: key.split("/").map(Number), detail: false }]);
    });
  }
  return out;
}

const weigh = (entries, camera) => entries.reduce((total, [, fade]) => total + S.starShare(fade, camera, ROOT, FINEST_LEVEL), 0);

test("approaching a star field 9% at a time never adds or drops a burst of stars", () => {
  const stars = field(30000, 7);
  const look = [8000, 0, 0];
  const direction = [0.6, -0.5, 0.62].map((v, i, all) => v / Math.hypot(...all));
  let distance = 20000;
  let worstRise = 0;
  let worstFall = 0;
  let steps = 0;
  let previous = null;
  while (distance > 30) {
    const camera = look.map((v, i) => v + direction[i] * distance);
    const entries = drawn(stars, 60, camera, look);
    const total = weigh(entries, camera);
    if (previous !== null && previous >= 50) {
      worstRise = Math.max(worstRise, (total - previous) / previous);
      worstFall = Math.max(worstFall, (previous - total) / previous);
      steps++;
    }
    previous = total;
    distance /= 2 ** (1 / 8);
  }
  assert.ok(steps > 20, `only ${steps} steps had 50 stars' worth on screen`);
  // Every halving of the distance brings eight times the tiles, so the stars on screen
  // rise by up to about 1.9 times over a 9% step in the busiest octave; the magnitude
  // limit (static/starmagnitude.js) holds the count on screen. They never drop in a burst.
  assert.ok(worstRise < 1.2, `one 9% zoom step raised the opacity on screen by ${(100 * worstRise).toFixed(0)}%`);
  assert.ok(worstFall < 0.25, `one 9% zoom step dropped the opacity on screen by ${(100 * worstFall).toFixed(0)}%`);
});

test("a star listed in a tile and its parent keeps its whole share through the hand-over", () => {
  // The brightest star of the field is in every tile above it.
  const stars = [{ id: 0, x: 8000.5, y: 0.5, z: 0.5, luminosity: 1e6 }].concat(field(2000, 3));
  const look = [8000, 0, 0];
  let lowest = 1;
  for (let distance = 8000; distance > 40; distance /= 1.05) {
    const camera = [look[0] + distance, look[1], look[2] + distance * 0.3];
    const entries = drawn(stars, 60, camera, look).filter(([id]) => id === 0);
    assert.ok(entries.length >= 1);
    lowest = Math.min(lowest, weigh(entries, camera));
  }
  assert.ok(lowest > 0.7, `the star dimmed to ${lowest.toFixed(2)} of its brightness`);
});

