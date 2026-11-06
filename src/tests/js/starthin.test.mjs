// tests/js/starthin.test.mjs -- static/starthin.js (MAP.366): a tile that reached its star cap keeps the density on screen
// near that of the bigger tiles in the view, fading its dimmer stars out by rank.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const T = await importPage("starthin.js");

const star = (luminosity, extra = {}) => Object.assign({ luminosity_sol: luminosity }, extra);

test("a node's stars are ranked by luminosity and a phenomenon gets no rank", () => {
  const stars = [star(1), star(100), star(10), star(5, { phenomenon: true })];
  const found = T.rankStars(stars, { m: 27 });
  assert.deepEqual(stars.map((s) => s.nodeRank), [2, 0, 1, undefined]);
  assert.deepEqual(found, { listed: 3, capped: false });
});

test("a block node is capped once either of its lists reaches the cap, a sector node never", () => {
  const bright = Array.from({ length: T.NODE_LIST_CAP }, (_, i) => star(i + 1));
  assert.equal(T.rankStars(bright.map((s) => ({ ...s })), { m: 3 }).capped, true);
  assert.equal(T.rankStars(bright.map((s) => ({ ...s })), { m: 1 }).capped, false);
  const generated = Array.from({ length: T.NODE_LIST_CAP }, (_, i) => star(i + 1, { generated: true }));
  assert.equal(T.rankStars(generated, { m: 9 }).capped, true);
  const mixed = bright.slice(0, 50).concat(generated.slice(0, 50).map((s) => ({ ...s })));
  assert.equal(T.rankStars(mixed, { m: 9 }).capped, false, "50 and 50 reach neither cap");
});

test("the reference size is the quantile of the capped tiles' sizes, at least the minimum, none without capped tiles", () => {
  assert.equal(T.referencePx([]), null);
  assert.equal(T.referencePx([10, 20, 30]), T.REFERENCE_MIN_PX);
  const sizes = [200, 400, 600, 800];
  assert.equal(T.referencePx(sizes), 600, "the 0.6 quantile of four is the third");
  assert.equal(T.referencePx([5000, 6000]), T.PX_CAP, "a tile over the screen counts as the cap");
});

test("a tile keeps the square of its size against the reference, between the floor and one", () => {
  assert.equal(T.keptShare(800, null), 1);
  assert.equal(T.keptShare(800, 800), 1);
  assert.equal(T.keptShare(1200, 800), 1);
  assert.ok(Math.abs(T.keptShare(560, 800) - 0.49) < 1e-12);
  assert.equal(T.keptShare(100, 800), T.DENSITY_FLOOR);
});

test("the stars past the kept count fade out by rank and the first few are never lost", () => {
  assert.equal(T.rankShare(95, 96, 1), 1, "a node that keeps everything keeps every rank");
  const kept = 0.5;
  const shares = Array.from({ length: 96 }, (_, rank) => T.rankShare(rank, 96, kept));
  assert.equal(shares[0], 1);
  assert.ok(shares[47] > 0 && shares[47] < 1, "the edge of the kept count fades");
  assert.equal(shares[95], 0);
  for (let i = 1; i < shares.length; i++) assert.ok(shares[i] <= shares[i - 1], "never rises with rank");
  const smallest = Array.from({ length: 96 }, (_, rank) => T.rankShare(rank, 96, T.DENSITY_FLOOR));
  assert.equal(smallest.slice(0, T.KEEP_MIN).every((s) => s > 0.5), true, "at least KEEP_MIN stars stay");
});

test("two capped tiles a fifth apart in size differ in density by the floor's share of the plain ratio", () => {
  const reference = T.referencePx([200, 400, 1000, 1000]);
  assert.equal(reference, 1000);
  const density = (px) => 96 * T.keptShare(px, reference) / (px * px);
  const plain = (px) => 96 / (px * px);
  assert.ok(Math.abs(plain(200) / plain(1000) - 25) < 1e-9, "unthinned, 25 times as thick");
  assert.ok(Math.abs(density(200) / density(1000) - 7.5) < 1e-9, "thinned to the floor (0.3), 7.5 times");
});
