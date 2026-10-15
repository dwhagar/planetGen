// tests/js/gridopacity.test.mjs -- static/gridopacity.js (MAP.188): the share that keeps
// the blocks' faces and lines from adding up to a near-opaque cube grid.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const G = await importPage("gridopacity.js");

// A grid of `n` x `n` x `layers` cells, 10 across and 10 high, all at `opacity`, about the equator
// and far from the center (so the cells are nearly cubes): n radii, n longitudes, `layers`
// colatitudes.
function grid(n, layers, opacity) {
  const cells = [];
  const r = 100000;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      for (let k = 0; k < layers; k++) {
        const r0 = r + 10 * i;
        const c0 = Math.PI / 2 - 0.0001 * (k + 1);
        cells.push({ r0, r1: r0 + 10, c0, c1: c0 + 0.0001, t0: (j * 10) / r, t1: ((j + 1) * 10) / r, opacity });
      }
    }
  }
  return cells;
}

test("no cells give no stats and one crossing", () => {
  assert.equal(G.cellStats([]), null);
  assert.equal(G.crossings(null, [0, 0, -1]), 1);
});

test("a ray along a stack of layers crosses each of them, a flat slab seen from above one", () => {
  const stack = G.cellStats(grid(10, 10, 0.3));
  const slab = G.cellStats(grid(10, 1, 0.3));
  assert.ok(Math.abs(G.crossings(stack, [0, 0, -1]) - 10) < 0.5, String(G.crossings(stack, [0, 0, -1])));
  assert.ok(Math.abs(G.crossings(slab, [0, 0, -1]) - 1) < 0.2);
  // Along the slab it crosses its cells one after another.
  assert.ok(G.crossings(slab, [1, 0, 0]) > 8);
  // A slant sits between.
  const slanted = G.crossings(stack, [Math.SQRT1_2, 0, -Math.SQRT1_2]);
  assert.ok(slanted > 10 && slanted < 21, String(slanted));
});

test("the share leaves one face alone and keeps many under the cap", () => {
  assert.equal(G.gridScale(1, 0.3), 1);
  assert.equal(G.gridScale(5, 0), 1);
  for (const n of [2, 5, 10, 25, 60]) {
    for (const a of [0.1, 0.3, 0.5]) {
      const s = G.gridScale(n, a);
      assert.ok(s > 0 && s <= 1);
      assert.ok(G.accumulated(n, a * s) <= G.GRID_ACCUMULATION_CAP + 1e-9, `n=${n} a=${a}`);
    }
  }
  // Barely visible when deep: 60 cells of 30% faces are scaled well below half.
  assert.ok(G.gridScale(60, 0.3) < 0.1);
  // Unscaled, they would be near opaque.
  assert.ok(G.accumulated(60, 0.3) > 0.99);
});

test("the deepest free view of a filled stage stays under the cap", () => {
  const stats = G.cellStats(grid(12, 12, 0.5));
  for (const direction of [[0, 0, -1], [1, 0, 0], [0.6, 0.6, -0.52], [0.3, 0.2, -0.93]]) {
    const n = G.crossings(stats, direction);
    assert.ok(G.accumulated(n, stats.maxOpacity * G.gridScale(n, stats.maxOpacity)) <= G.GRID_ACCUMULATION_CAP + 1e-9);
  }
});
