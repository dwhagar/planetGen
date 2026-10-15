// tests/js/galaxytilecut.test.mjs -- static/galaxytilecut.js (MAP.152): the
// distance cut of tiles the Galaxy Map fetches, key for key what
// planetgen.galaxy.viewport computes (the cases come from Python in
// PLANETGEN_JS_FIXTURES), and the properties that make it a cut.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const C = await importPage("galaxytilecut.js");
const F = await importPage("starfade.js");
const fixtures = JSON.parse(process.env.PLANETGEN_JS_FIXTURES).tileCut;

const config = (extra = {}) => Object.assign({ rootEdge: fixtures.root, finestLevel: fixtures.finestLevel }, extra);

test("the constants equal Python's", () => {
  assert.equal(C.DISTANCE_CUT_EDGES, fixtures.edges);
  assert.equal(F.CUT_EDGES, fixtures.edges);
  assert.equal(F.FADE_BAND, fixtures.band);
});

test("the cut and the fetch set are Python's, key for key", () => {
  for (const one of fixtures.cases) {
    const cfg = config({ bounds: one.bounds });
    assert.deepEqual(C.distanceCutTiles(one.camera, one.look, cfg), one.cut, JSON.stringify(one.camera));
    assert.deepEqual(C.fetchCutTiles(one.camera, one.look, cfg), one.fetch, JSON.stringify(one.camera));
  }
});

test("a cut with no cone covers the bounds once: no tile holds another", () => {
  const camera = [8000, 0, 0];
  const keys = C.distanceCutTiles(camera, camera, config({
    bounds: [[-12600, -12600, -960], [12600, 12600, 960]], halfAngleDeg: 180,
  }));
  assert.ok(keys.length > 50);
  const boxes = keys.map((key) => C.tileBoxOf(key, fixtures.root));
  boxes.forEach((a, i) => boxes.forEach((b, j) => {
    if (i >= j) return;
    const overlap = [0, 1, 2].every((axis) => a.lo[axis] < b.lo[axis] + b.edge && b.lo[axis] < a.lo[axis] + a.edge);
    assert.ok(!overlap, `${keys[i]} overlaps ${keys[j]}`);
  }));
});

test("tiles get smaller toward the camera", () => {
  const camera = [8000, 0, 0];
  const keys = C.distanceCutTiles(camera, camera, config({
    bounds: [[-12600, -12600, -960], [12600, 12600, 960]], halfAngleDeg: 180,
  }));
  keys.forEach((key) => {
    const box = C.tileBoxOf(key, fixtures.root);
    const gap = [0, 1, 2].map((a) => Math.max(box.lo[a] - camera[a], camera[a] - (box.lo[a] + box.edge), 0));
    const distance = Math.hypot(...gap);
    const level = Number(key.split("/")[0]);
    if (level < fixtures.finestLevel) assert.ok(box.edge <= C.DISTANCE_CUT_EDGES * distance + 1e-6, key);
  });
});

test("the cone leaves out tiles behind the camera", () => {
  const camera = [8000, 0, 0];
  const bounds = [[-12600, -12600, -960], [12600, 12600, 960]];
  const forward = C.distanceCutTiles(camera, [camera[0] + 1, 0, 0], config({ bounds }));
  const all = C.distanceCutTiles(camera, [camera[0] + 1, 0, 0], config({ bounds, halfAngleDeg: 180 }));
  assert.ok(forward.length < all.length * 0.6);
  assert.ok(forward.every((key) => all.includes(key)));
});
