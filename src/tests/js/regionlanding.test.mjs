// tests/js/regionlanding.test.mjs -- static/regionlanding.js (MAP.336): the cluster a "densest" region landing frames.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const L = await import(new URL("regionlanding.js", STATIC_URL).href);

const EDGE = 4;
const sector = (x, y, z, massKg) => ({ x, y, z, sizePc: EDGE, massKg });
const region = { center: [0, 0, 0], radius: 40 };

test("a region with no ranked sector has no cluster, so the landing falls back to the frame", () => {
  assert.equal(L.densestCluster([], region, EDGE), null);
  assert.equal(L.densestCluster([sector(0, 0, 0, 0)], region, EDGE), null);
  assert.equal(L.densestCluster([sector(500, 0, 0, 9)], region, EDGE), null);
  assert.equal(L.densestCluster([sector(0, 0, 0, 1)], null, EDGE), null);
  assert.equal(L.densestCluster([sector(0, 0, 0, 1)], { center: [0, 0, 0], radius: 0 }, EDGE), null);
});

test("one ranked sector is framed alone, its cube inside the sphere", () => {
  const found = L.densestCluster([sector(8, 0, 0, 5)], region, EDGE);
  assert.deepEqual(found.center, [8, 0, 0]);
  assert.equal(found.count, 1);
  assert.ok(found.radius >= 0.5 * Math.sqrt(3) * EDGE - 1e-9);
});

test("the cluster is the heavier group, not the single heaviest sector far from the rest", () => {
  const lone = sector(-30, 0, 0, 10);
  const group = [sector(10, 0, 0, 6), sector(14, 0, 0, 6), sector(10, 4, 0, 6)];
  const found = L.densestCluster([lone].concat(group), region, EDGE);
  assert.equal(found.count, 3);
  assert.ok(found.center[0] > 9 && found.center[0] < 15);
  assert.equal(found.massKg, 18);
  // Every member's cube is inside the sphere that frames it.
  group.forEach((s) => assert.ok(Math.hypot(s.x - found.center[0], s.y - found.center[1], s.z - found.center[2]) < found.radius));
});

test("the middle is weighted by mass and sectors outside the region are left out", () => {
  const found = L.densestCluster([sector(0, 0, 0, 3), sector(4, 0, 0, 1), sector(900, 0, 0, 100)], region, EDGE);
  assert.equal(found.count, 2);
  assert.ok(Math.abs(found.center[0] - 1) < 1e-9);
});
