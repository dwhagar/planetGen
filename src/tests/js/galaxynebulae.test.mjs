// tests/js/galaxynebulae.test.mjs -- static/galaxynebulae.js (MAP.176): the nebula regions' blobs, look and budget.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const N = await import(new URL("galaxynebulae.js", STATIC_URL).href);
const L = await import(new URL("nebulalook.js", STATIC_URL).href);

const looks = L.nebulaLooks(false);
const block = { m: 81, shell: 4, band: 1, wedge: 2 };

const region = (key, extra) => Object.assign({
  key, count: 12, volume_pc3: 5000, centre: [100, 20, 5], radius_pc: 40, shares: { emission: 0.7, reflection: 0.3 }, class: "H",
}, extra);
const node = (regions) => ({ key: "81.4.1.2", block, node: { nebula_regions: regions } });

test("a region's blob is the group's sphere in the colour of its dominant family, with a ring when it holds several", () => {
  const [blob] = N.regionBlobs([node([region("27.12.3.1")])], looks, false, null);
  assert.deepEqual([blob.x, blob.y, blob.z, blob.radiusPc, blob.count], [100, 20, 5, 40, 12]);
  assert.equal(blob.key, "27.12.3.1");
  assert.equal(blob.family, "emission");
  assert.deepEqual(blob.block, block);
  assert.deepEqual(blob.color, N.colourOf(looks.emission[0]));
  assert.equal(blob.ring, true);
  const [one] = N.regionBlobs([node([region("27.1.1.1", { count: 1 })])], looks, false, null);
  assert.equal(one.ring, false);
});

test("a malformed group, a repeated key and a hidden family are left out", () => {
  const list = N.regionBlobs([node([
    region("a"), region("a"), region("b", { centre: null }), region("c", { radius_pc: 0 }), region("d", { count: 0 }),
    region("e", { shares: { dark: 1 } }), region("f", { shares: {} }), null,
  ])], looks, false, (family) => family !== "dark");
  assert.deepEqual(list.map((b) => b.key), ["a"]);
});

test("the dominant family is the largest share, a tie going to the earlier family", () => {
  assert.equal(N.dominantFamily({ emission: 0.2, dark: 0.8 }), "dark");
  assert.equal(N.dominantFamily({ diffuse: 0.5, emission: 0.5 }), "diffuse");
  assert.equal(N.dominantFamily({}), null);
  assert.equal(N.dominantFamily(null), null);
});

test("every family has a colour on both themes and a label", () => {
  [false, true].forEach((light) => N.FAMILIES.forEach((family) => {
    const colour = N.familyColour(family, L.nebulaLooks(light), light);
    assert.equal(colour.length, 3);
    colour.forEach((v) => assert.ok(v >= 0 && v <= 1));
    assert.ok(N.FAMILY_LABELS[family]);
  }));
  assert.notDeepEqual(N.familyColour("emission", looks, false), N.familyColour("reflection", looks, false));
});

test("a region shows from 2 px across, is stronger with more nebulae and is gone when it fills the screen", () => {
  const at = (distance, extra) => N.regionLook(Object.assign({ x: distance, y: 0, z: 0, radiusPc: 40, count: 12 }, extra), [0, 0, 0], 1000);
  assert.equal(at(1e6).alpha, 0, "a speck under a pixel is not drawn");
  assert.ok(at(40000 / 1.5).alpha > 0, "from a pixel and a half of radius");
  assert.ok(at(500).alpha > 0.1);
  assert.ok(at(500, { count: 200 }).alpha > at(500, { count: 2 }).alpha);
  assert.ok(at(500, { count: 200 }).alpha <= N.ALPHA[1] + 1e-9);
  assert.equal(at(40000 / N.RESOLVE_PX[1] / 2).alpha, 0, "inside its nebulae the haze is gone");
  const near = at(40000 / 100);
  assert.ok(near.radiusPx >= N.RADIUS_PX[0] && near.radiusPx <= N.RADIUS_PX[1]);
  assert.equal(at(0).radiusPx, N.RADIUS_PX[1], "a camera in the middle gets the largest allowed radius, not infinity");
});

test("the budget keeps the strongest regions, fades the last tenth and never changes what fits", () => {
  const small = [0.3, 0, 0.2, 0.1];
  assert.deepEqual(N.cutToBudget(small.slice(), 4), small);
  const alphas = [];
  for (let i = 0; i < 20; i++) alphas.push(0.5 - i * 0.01);
  const cut = N.cutToBudget(alphas.slice(), 10);
  assert.equal(cut.filter((a) => a > 0).length, 10);
  assert.deepEqual(cut.slice(0, 9), alphas.slice(0, 9), "the strong ones are untouched");
  assert.ok(cut[9] > 0 && cut[9] < alphas[9], "the last rank fades towards the cut");
  assert.deepEqual(cut.slice(10), new Array(10).fill(0));
});

test("3,000 candidates: exactly the strongest 600 are drawn, the same ones every time, and quickly", () => {
  let seed = 12345;
  const random = () => {
    seed = (seed * 1103515245 + 12345) % 2147483648;
    return seed / 2147483648;
  };
  const regions = [];
  for (let i = 0; i < 3000; i++) {
    regions.push(region("27." + i + ".0.0", {
      centre: [(random() - 0.5) * 20000, (random() - 0.5) * 20000, (random() - 0.5) * 600],
      radius_pc: 8 + random() * 120, count: 2 + Math.floor(random() * 300),
      shares: { [N.FAMILIES[i % N.FAMILIES.length]]: 1 },
    }));
  }
  const blobs = N.regionBlobs([node(regions)], looks, false, null);
  assert.equal(blobs.length, 3000);
  const camera = [0, 0, 9000];
  const alphasOf = () => blobs.map((blob) => N.regionLook(blob, camera, 800).alpha);
  const drawable = alphasOf().filter((a) => a > 0).length;
  assert.ok(drawable > 600, "the fixture really has more than the budget to draw: " + drawable);
  const began = process.hrtime.bigint();
  const first = N.cutToBudget(alphasOf(), N.REGION_BUDGET);
  const millis = Number(process.hrtime.bigint() - began) / 1e6;
  assert.ok(millis < 500, "the cut of 3,000 took " + millis + " ms");
  const shown = first.filter((a) => a > 0);
  assert.equal(shown.length, N.REGION_BUDGET, "exactly the budget is drawn");
  const raw = alphasOf();
  const strongest = raw.map((a, i) => [a, i]).sort((p, q) => (q[0] - p[0]) || (p[1] - q[1])).slice(0, N.REGION_BUDGET).map((p) => p[1]).sort((p, q) => p - q);
  assert.deepEqual(first.map((a, i) => (a > 0 ? i : -1)).filter((i) => i >= 0), strongest, "they are the strongest ones");
  assert.deepEqual(N.cutToBudget(alphasOf(), N.REGION_BUDGET), first, "the same input draws the same regions");
});

test("a region's words say how many nebulae and of what family", () => {
  assert.equal(N.regionText({ count: 12, family: "emission" }), "12 nebulae, mostly emission nebulae");
  assert.equal(N.regionText({ count: 1, family: "supernova-remnant" }), "1 nebula, mostly supernova remnants");
});
