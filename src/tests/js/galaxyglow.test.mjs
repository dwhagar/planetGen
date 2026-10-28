// tests/js/galaxyglow.test.mjs -- static/galaxyglow.js (MAP.266): the unresolved glow's blobs and how strongly they show.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const G = await importPage("galaxyglow.js");
const N = await importPage("galaxynodecut.js");

const EDGE = 4;
const block = { m: 243, shell: 3, band: 1, wedge: 2 };
const box = N.geometry(block, EDGE);
const middle = [0, 1, 2].map((a) => (box.lo[a] + box.hi[a]) / 2);
const star = (lum, x, y, z) => ({ luminosity_sol: lum, x, y, z });

test("a node with no aggregate light has no blob", () => {
  assert.equal(G.blobOf("k", block, { stars: [star(10, 0, 0, 0)], aggregate: {} }, EDGE), null);
  assert.equal(G.blobOf("k", block, null, EDGE), null);
});

test("the blob is the stored light less the listed stars, centred on what is left", () => {
  const a = star(100, middle[0] - 10, middle[1], middle[2]);
  const b = star(300, middle[0] + 10, middle[1], middle[2]);
  // 1,000 in all, centre at the middle: the listed 400 weigh the centre to x of (100 * -10 + 300 * 10) / 400 = +5
  const node = {
    stars: [a, b],
    aggregate: { lights: { star_luminosity_sol: 1000, star_centre: middle, star_temperature_k: 4100 } },
  };
  const blob = G.blobOf("k", block, node, EDGE);
  assert.equal(blob.luminosity, 600);
  assert.equal(blob.temperatureK, 4100);
  // 1000 * middle - (100 a + 300 b) over 600: the unlisted part sits where the listed part is not.
  const expectedX = (1000 * middle[0] - (100 * a.x + 300 * b.x)) / 600;
  assert.ok(Math.abs(blob.x - expectedX) < 1e-9);
  assert.ok(Math.abs(blob.y - middle[1]) < 1e-9 && Math.abs(blob.z - middle[2]) < 1e-9);
});

test("generated sectors' light joins at the middle of the block, less the generated stars listed", () => {
  const node = {
    stars: [], generated: [star(50, 0, 0, 0)],
    aggregate: { stats: { luminosity_sol: 250 } },
  };
  const blob = G.blobOf("k", block, node, EDGE);
  assert.equal(blob.luminosity, 200);
  assert.deepEqual([blob.x, blob.y, blob.z].map((v, a) => Math.abs(v - middle[a]) < 1e-9), [true, true, true]);
  assert.equal(blob.temperatureK, G.GLOW_DEFAULT_TEMPERATURE_K);
});

test("listed stars brighter than the stored total leave nothing, never a negative light", () => {
  const node = { stars: [star(500, middle[0], middle[1], middle[2])],
    aggregate: { lights: { star_luminosity_sol: 400, star_centre: middle } } };
  assert.equal(G.blobOf("k", block, node, EDGE), null);
});

test("a blob is brighter nearer, never past the cap, and faint far away", () => {
  const blob = { x: 0, y: 0, z: 0, luminosity: 1e8, radiusPc: 500, temperatureK: 5000 };
  const focal = 900;
  const near = G.glowLook(blob, [0, 0, 5000], focal);
  const far = G.glowLook(blob, [0, 0, 5e6], focal);
  assert.ok(near.alpha > far.alpha);
  assert.ok(near.alpha <= G.GLOW_ALPHA_MAX + 1e-12);
  assert.ok(far.alpha < 0.01);
  assert.ok(near.radiusPx <= G.GLOW_RADIUS_PX[1] && far.radiusPx >= G.GLOW_RADIUS_PX[0]);
  assert.ok(Math.abs(near.radiusPx - Math.min(G.GLOW_RADIUS_PX[1], 500 * focal / 5000)) < 1e-9);
});

test("a glow fades as the stars it stands for resolve and is gone when its block fills the screen (MAP.324)", () => {
  const blob = { x: 0, y: 0, z: 0, luminosity: 1e8, radiusPc: 100, temperatureK: 5000 };
  const focal = 900;
  // The camera at a distance that puts the block's radius at the given pixels.
  const at = (px) => G.glowLook(blob, [0, 0, (blob.radiusPc * focal) / px], focal);
  const [start, end] = G.GLOW_RESOLVE_PX;
  const small = at(start * 0.5);
  assert.ok(small.alpha > 0.9 * G.GLOW_ALPHA_MAX, "an unresolved far glow keeps its strength");
  const middle = at((start + end) / 2);
  assert.ok(middle.alpha < small.alpha && middle.alpha > 0, "it fades in between");
  assert.equal(at(end).alpha, 0);
  assert.equal(at(end * 4).alpha, 0, "and stays gone at the zooms where the haze flooded the map");
  let last = Infinity;
  for (let px = start * 0.5; px <= end * 1.2; px += 20) {
    const alpha = at(px).alpha;
    assert.ok(alpha <= last + 1e-12, "never stronger as the block grows");
    last = alpha;
  }
});

test("a glow is never stronger than the cap, and the cap leaves a star's core above it (MAP.324)", () => {
  assert.ok(G.GLOW_ALPHA_MAX <= 0.4);
  const blob = { x: 0, y: 0, z: 0, luminosity: 1e12, radiusPc: 10, temperatureK: 5000 };
  assert.ok(G.glowLook(blob, [0, 0, 100], 900).alpha <= G.GLOW_ALPHA_MAX + 1e-12);
});

test("a hand-over keeps the light: a parent's share and its children's add up to one", () => {
  const blob = { luminosity: 1000 };
  for (const share of [0, 0.25, 0.5, 1]) {
    assert.equal(G.sharedLight(blob, share) + G.sharedLight(blob, 1 - share), 1000);
  }
});

test("blobsOf keeps the order of the nodes and leaves out the ones with no glow", () => {
  const lit = { stars: [], aggregate: { stats: { luminosity_sol: 10 } } };
  const dark = { stars: [], aggregate: {} };
  const found = G.blobsOf([{ key: "a", block, node: lit }, { key: "b", block, node: dark }, { key: "c", block, node: lit }], EDGE);
  assert.deepEqual(found.map((b) => b.key), ["a", "c"]);
});
