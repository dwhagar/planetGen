// tests/js/lumpymesh.test.mjs -- static/lumpymesh.js (MAP.279): which bodies are drawn as lumps, that the shape is a pure function of the
// look record, and that the mesh keeps the mean radius of the unit sphere.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";
import * as THREE from "../../html/static/vendor/three.module.min.js";

const L = await importPage("lumpymesh.js");

const ROCK = { v: 1, kind: "moon", body: "solid", seed: 424242, color: "#8a8478", radiusKm: 12 };

test("a body under roundness 1 is a lump, a round one is not, and a missing ratio falls back on the radius", () => {
  assert.equal(L.isLumpy({ ...ROCK, roundness: 0.3 }), true);
  assert.equal(L.isLumpy({ ...ROCK, roundness: 1 }), false);
  assert.equal(L.isLumpy({ ...ROCK, roundness: 2.4, radiusKm: 150 }), false);
  assert.equal(L.isLumpy({ ...ROCK, radiusKm: 120 }), true);
  assert.equal(L.isLumpy({ ...ROCK, radiusKm: 2000 }), false);
  assert.equal(L.isLumpy({ ...ROCK, roundness: 0.2, body: "gas" }), false);
  assert.equal(L.isLumpy(null), false);
});

test("the shape is the same for the same record and differs between seeds", () => {
  const a = L.shapeOf({ ...ROCK, roundness: 0.2 });
  const b = L.shapeOf({ ...ROCK, roundness: 0.2 });
  assert.deepEqual(a, b);
  const other = L.shapeOf({ ...ROCK, roundness: 0.2, seed: 7 });
  assert.notDeepEqual(a.axes, other.axes);
});

test("the axes keep the volume of the unit sphere and stretch more the farther the body is from round", () => {
  for (const roundness of [0.9, 0.5, 0.05]) {
    const { axes } = L.shapeOf({ ...ROCK, roundness });
    assert.ok(Math.abs(axes[0] * axes[1] * axes[2] - 1) < 1e-9);
    assert.ok(axes[0] >= axes[1] && axes[1] >= axes[2]);
  }
  const near = L.shapeOf({ ...ROCK, roundness: 0.95 });
  const far = L.shapeOf({ ...ROCK, roundness: 0.02 });
  assert.ok(far.axes[0] / far.axes[2] > near.axes[0] / near.axes[2]);
  assert.ok(far.axes[0] / far.axes[2] <= 1 + L.MAX_ELONGATION + 1e-9);
  assert.ok(far.amplitude > near.amplitude);
});

test("the surface stays within the extent and never folds through the centre", () => {
  for (const seed of [1, 99, 123456789, 4000000000]) {
    const shape = L.shapeOf({ ...ROCK, roundness: 0.01, seed });
    for (let i = 0; i < 400; i++) {
      const z = 2 * ((i * 0.6180339887) % 1) - 1;
      const turn = i * 2.399963;
      const flat = Math.sqrt(1 - z * z);
      const p = L.surfacePoint(shape, flat * Math.cos(turn), flat * Math.sin(turn), z);
      const r = Math.hypot(...p);
      assert.ok(r > 0.3 && r <= shape.extent + 1e-9, `radius ${r} outside 0.3 to ${shape.extent}`);
    }
  }
});

test("the geometry is a closed sphere of vertices pushed to the surface, with normals", () => {
  const look = { ...ROCK, roundness: 0.1 };
  const { geometry, extent } = L.lumpyGeometry(THREE, look, 24, 16);
  const position = geometry.attributes.position;
  const shape = L.shapeOf(look);
  let longest = 0;
  for (let i = 0; i < position.count; i++) {
    longest = Math.max(longest, Math.hypot(position.getX(i), position.getY(i), position.getZ(i)));
  }
  assert.ok(longest <= extent + 1e-9 && longest > 1 - 1e-9);
  assert.equal(extent, shape.extent);
  assert.ok(geometry.attributes.normal.count === position.count);
  assert.ok(geometry.attributes.uv.count === position.count);
  const again = L.lumpyGeometry(THREE, look, 24, 16).geometry.attributes.position;
  assert.deepEqual(Array.from(again.array), Array.from(position.array));
});
