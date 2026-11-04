// tests/js/maplogruler.test.mjs -- static/maplogruler.js (MAP.336): the named scales of the map's log ruler.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const R = await import(new URL("maplogruler.js", STATIC_URL).href);

const FAR = 36000;
const FLOOR = 1e-10;

test("the stops run from the whole galaxy at the top to a moon at the bottom, each in its place", () => {
  const stops = R.stopsBetween(FAR, FLOOR);
  assert.deepEqual(stops.map((s) => s.key), ["galaxy", "arm", "region", "sector", "neighbourhood", "system", "planet", "moon"]);
  assert.equal(stops[0].share, 0);
  assert.ok(Math.abs(stops[7].share - 1) < 1e-12);
  for (let n = 1; n < stops.length; n++) assert.ok(stops[n].share > stops[n - 1].share);
  // A decade is the same length of ruler anywhere: sector (4 pc) to neighbourhood (1 pc) is log10(4) decades.
  const decade = (stops[4].share - stops[3].share) / Math.log10(4);
  assert.ok(Math.abs(decade - (stops[5].share - stops[4].share) / 3) < 1e-9);
});

test("a stop nearer than the floor or farther than the galaxy is left off", () => {
  const keys = R.stopsBetween(500, 1e-6).map((s) => s.key);
  assert.deepEqual(keys, ["galaxy", "region", "sector", "neighbourhood", "system"]);
});

test("a share and a distance are each other's inverse", () => {
  [0, 0.25, 0.5, 0.9, 1].forEach((share) => {
    assert.ok(Math.abs(R.shareOf(R.distanceAt(share, FAR, FLOOR), FAR, FLOOR) - share) < 1e-12);
  });
  assert.ok(Math.abs(R.distanceAt(2, FAR, FLOOR) / FLOOR - 1) < 1e-9);
  assert.equal(R.shareOf(1e9, FAR, FLOOR), 0);
});

test("the readout names the distance in the reader's units and the nearest thing along the line of sight", () => {
  const pc = (v) => v + " pc (" + (v * 3.2616).toFixed(1) + " ly)";
  assert.equal(R.readoutText(2.4, null, pc), "Distance from focus: 2.4 pc (7.8 ly)");
  assert.equal(R.readoutText(2.4, { name: "Kanbee 3", distance: 0.9 }, pc), "Distance from focus: 2.4 pc (7.8 ly); nearest: Kanbee 3, 0.9 pc (2.9 ly)");
  assert.equal(R.readoutText(1e-5, { name: "", distance: 1e-5 }, pc), "Distance from focus: 0.00001 pc (0.0 ly)");
  assert.equal(R.scaleNamed(3, R.stopsBetween(FAR, FLOOR)).key, "sector");
  assert.equal(R.scaleNamed(1e7, R.stopsBetween(FAR, FLOOR)).key, "galaxy");
  assert.equal(R.scaleNamed(1e-10, R.stopsBetween(FAR, FLOOR)).key, "moon");
});

test("a surface stop stands at the new floor, 30 km, and is left off a ruler that stops short of it", () => {
  const deep = R.stopsBetween(FAR, 1e-12);
  assert.equal(deep[deep.length - 1].key, "surface");
  assert.ok(Math.abs(deep[deep.length - 1].share - 1) < 1e-12);
  assert.equal(R.scaleNamed(1e-12, deep).key, "surface");
  assert.ok(!R.stopsBetween(FAR, FLOOR).some((s) => s.key === "surface"));
});
