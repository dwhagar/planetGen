// tests/js/galaxyregions.test.mjs -- static/galaxyregions.js (MAP.316): the charted regions' blobs, look, order and names.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const R = await import(new URL("galaxyregions.js", STATIC_URL).href);

const region = (key, share, extra) => Object.assign({
  key, m: 9, shell: 1, band: 0, wedge: 1, centre: [100, 0, 0], radius: 12, sector_count: 3, capacity: 780, share,
  name: "Fardoo Beez Region", parent: "Fardoo Reach",
}, extra);

test("a region's blob is its extent sphere, with the share's color", () => {
  const [blob] = R.regionBlobs([region("9.1.0.1", 0.14)]);
  assert.deepEqual([blob.x, blob.y, blob.z, blob.radiusPc], [100, 0, 0, 12]);
  assert.equal(blob.key, "9.1.0.1");
  assert.deepEqual(blob.block, { m: 9, shell: 1, band: 0, wedge: 1 });
  assert.equal(blob.name, "Fardoo Beez Region");
  assert.deepEqual(blob.color, R.shareColor(0.14));
});

test("a region without an extent is left out, and every share above zero is kept", () => {
  const list = R.regionBlobs([region("a", 0.0001), region("b", 0.5, { centre: null }), region("c", 0.2, { radius: 0 }), null]);
  assert.deepEqual(list.map((b) => b.key), ["a"]);
});

test("the color goes from cool for a lone sector to warm for a filled region, in range", () => {
  const cool = R.shareColor(0.0013);
  const warm = R.shareColor(1);
  assert.ok(cool[2] > cool[0], "a lone sector is blue-ish");
  assert.ok(warm[0] > warm[2], "a filled region is warm");
  [0, 0.01, 0.14, 0.5, 1, 2, -1].forEach((share) => R.shareColor(share).forEach((v) => assert.ok(v >= 0 && v <= 1)));
  assert.notDeepEqual(R.shareColor(0.14), R.shareColor(0.6));
});

test("a far region is a visible mark and a charted one is stronger than a lone sector", () => {
  const blob = R.regionBlobs([region("a", 0.01, { centre: [0, 0, 0] })])[0];
  const stronger = R.regionBlobs([region("b", 0.9, { centre: [0, 0, 0] })])[0];
  const far = R.regionLook(blob, [0, 0, 60000], 800);
  assert.equal(far.radiusPx, R.REGION_RADIUS_PX[0], "never under the smallest mark");
  assert.ok(far.alpha >= R.REGION_ALPHA[0]);
  assert.ok(R.regionLook(stronger, [0, 0, 60000], 800).alpha > far.alpha);
});

test("the sphere fades as the camera comes in and is gone once it fills the view", () => {
  const blob = R.regionBlobs([region("a", 0.5, { centre: [0, 0, 0], radius: 12 })])[0];
  const alphaAt = (distance) => R.regionLook(blob, [0, 0, distance], 800).alpha;
  const near = alphaAt(12 * 800 / 100);      // 100 px in radius: not yet fading
  const mid = alphaAt(12 * 800 / 350);       // 350 px: half way
  const inside = alphaAt(12 * 800 / 700);    // past 600 px
  assert.ok(near > mid && mid > 0, [near, mid]);
  assert.equal(inside, 0);
  assert.equal(R.regionLook(blob, [0, 0, 12 * 800 / 100], 800).alpha, near);
  assert.ok(R.regionLook(blob, [0, 0, 1], 800).radiusPx <= R.REGION_RADIUS_PX[1], "capped so it never washes the canvas");
});

test("the list sorts by share charted, or by name, the cell breaking ties", () => {
  const list = R.regionBlobs([
    region("9.1.0.1", 0.02, { name: "Zeta Beez Region" }),
    region("9.2.0.1", 0.5, { name: "Alfa Beez Region" }),
    region("9.3.0.1", 0.5, { name: "Mid Beez Region" }),
    region("9.4.0.1", 0.0004, { name: null }),
  ]);
  assert.deepEqual(R.sortRegions(list, "share").map((b) => b.key), ["9.2.0.1", "9.3.0.1", "9.1.0.1", "9.4.0.1"]);
  assert.deepEqual(R.sortRegions(list, "name").map((b) => b.key), ["9.2.0.1", "9.3.0.1", "9.1.0.1", "9.4.0.1"]);
  assert.deepEqual(R.sortRegions(list.slice().reverse(), "share").map((b) => b.key), ["9.2.0.1", "9.3.0.1", "9.1.0.1", "9.4.0.1"]);
  assert.equal(list[0].key, "9.1.0.1", "a copy");
});

test("the share reads in percent, never rounding a nonzero one down to nothing", () => {
  assert.equal(R.chartedText(0.14), "14 percent charted");
  assert.equal(R.chartedText(0.0013), "under 1 percent charted");
  assert.equal(R.chartedText(1), "100 percent charted");
});

const inhabited = (extra) => Object.assign({ systems: 2, worlds: 3, score: 1.2, tier: "Blue" }, extra);
const KIND = { id: "marker-inhabited", interest: "inhabited", icon: "inhabited-world" };
const alwaysLook = () => ({ alpha: 0.4 });
const yes = () => true;
const no = () => false;

test("a region's blob carries its inhabited worlds with the tier as an index, and none when it has none", () => {
  const [lived, empty, odd] = R.regionBlobs([
    region("a", 0.2, { inhabited: inhabited() }), region("b", 0.2), region("c", 0.2, { inhabited: { systems: 1, worlds: 1 } })]);
  assert.deepEqual(lived.inhabited, { systems: 2, worlds: 3, tier: 0 });
  assert.equal(empty.inhabited, null);
  assert.deepEqual(odd.inhabited, { systems: 1, worlds: 1, tier: null });
});

test("a region with inhabited worlds gets one token at its middle, in its best tier, named for the region", () => {
  const blobs = R.regionBlobs([region("a", 0.2, { inhabited: inhabited({ tier: "Green" }) }), region("b", 0.2)]);
  const found = R.regionMarkers(KIND, blobs, alwaysLook, yes, no);
  assert.equal(found.length, 1);
  const { kind, item } = found[0];
  assert.equal(kind.id, "marker-inhabited");
  assert.deepEqual([item.x, item.y, item.z, item.n, item.systems, item.tier, item.regionKey], [100, 0, 0, 3, 2, 1, "a"]);
  assert.equal(kind.describe(item), "Fardoo Beez Region, 3 inhabited worlds in 2 systems, best Green");
});

test("a region's token waits for its sphere to show, for the view, and for no system's own mark to stand for it", () => {
  const blobs = R.regionBlobs([region("a", 0.2, { inhabited: inhabited() })]);
  assert.equal(R.regionMarkers(KIND, blobs, () => ({ alpha: R.REGION_TOKEN_MIN_ALPHA / 2 }), yes, no).length, 0);
  assert.equal(R.regionMarkers(KIND, blobs, alwaysLook, no, no).length, 0);
  assert.equal(R.regionMarkers(KIND, blobs, alwaysLook, yes, yes).length, 0);
  assert.equal(R.regionMarkers(KIND, blobs, alwaysLook, yes, no).length, 1);
});

test("the inhabited text is singular for one and says nothing without a tier or a region", () => {
  assert.equal(R.inhabitedText({ worlds: 1, systems: 1, tier: null }), "1 inhabited world in 1 system");
  assert.equal(R.inhabitedText({ worlds: 2, systems: 1, tier: 2 }), "2 inhabited worlds in 1 system, best Yellow");
  assert.equal(R.inhabitedText(null), "");
});
