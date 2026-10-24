// tests/js/starfade.test.mjs -- static/starfade.js: the node fade that makes the Galaxy Map's stars come in and go out as
// the camera nears or leaves the node that lists them (MAP.153, redone for sector-block nodes in MAP.261): the
// arithmetic on its own, then a synthetic star field listed node by node the server's way and approached from afar,
// to bound how much the picture changes in one small step.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const S = await importPage("starfade.js");
const C = await importPage("galaxynodecut.js");
const G = await importPage("galaxyprisms.js");

const EDGE = 4;
const FOCAL = C.DEFAULT_FOCAL_PX;
const BOUNDS = [[-12600, -12600, -960], [12600, 12600, 960]];

test("a star in no node is whole, and a node's star takes the node's weight", () => {
  assert.equal(S.starShare(null, () => 0.25), 1);
  assert.equal(S.starShare(S.POINT_FADE, () => 0.25), 1);
  assert.equal(S.starShare({ block: { m: 3, shell: 0, band: 0, wedge: 0 } }, () => 0.25), 0.25);
  assert.equal(S.opacitySum([null, { block: { m: 3, shell: 0, band: 0, wedge: 0 } }], () => 0.5), 1.5);
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

// Stars strewn through a 4 kpc cube round (8000, 0, 0), luminosity a power law, each with the keys of the blocks over it.
function field(count, seed) {
  const random = mulberry32(seed);
  const stars = [];
  for (let i = 0; i < count; i++) {
    const star = {
      id: i, x: 8000 + (random() - 0.5) * 4000, y: (random() - 0.5) * 4000, z: (random() - 0.5) * 600,
      luminosity: Math.pow(random(), 5) * 1000,
    };
    const address = G.sectorAddressAt(star.x, star.y, star.z, EDGE);
    star.chain = G.drillChainOf(address.shell, address.band, address.wedge).map((block) => C.formatNodeKey(block));
    stars.push(star);
  }
  return stars;
}

// Which stars each node lists: the `cap` most luminous of those in it, found once.
function lister(stars, cap) {
  const byNode = new Map();
  for (const star of stars) {
    for (const key of star.chain) {
      if (!byNode.has(key)) byNode.set(key, []);
      byNode.get(key).push(star);
    }
  }
  for (const [key, list] of byNode) byNode.set(key, list.sort((p, q) => q.luminosity - p.luminosity).slice(0, cap));
  return byNode;
}

// What the page draws with the camera at `camera` looking at `look`: one entry [star id, fade] for every star every
// node of the cut lists, and the weight calculator for the camera.
function drawn(listing, camera, look, halfAngle) {
  const cut = C.nodeCut(camera, look, { edgePc: EDGE, focalPx: FOCAL, bounds: BOUNDS, halfAngleDeg: halfAngle });
  const entries = [];
  for (const key of cut.keys) {
    const fade = { block: C.parseNodeKey(key) };
    for (const star of listing.get(key) || []) entries.push([star.id, fade]);
  }
  return { entries, weigh: C.weightCalculator(camera, EDGE, FOCAL, cut.expanded) };
}

const total = ({ entries, weigh }) => entries.reduce((sum, [, fade]) => sum + S.starShare(fade, weigh), 0);

test("approaching a star field 9% at a time never adds or drops a burst of stars", () => {
  const listing = lister(field(30000, 7), 60);
  const look = [8000, 0, 0];
  const direction = [0.6, -0.5, 0.62].map((v, i, all) => v / Math.hypot(...all));
  let distance = 20000;
  let worstRise = 0;
  let worstFall = 0;
  let steps = 0;
  let previous = null;
  while (distance > 30) {
    const camera = look.map((v, i) => v + direction[i] * distance);
    const sum = total(drawn(listing, camera, look, C.CONE_HALF_ANGLE_DEG));
    if (previous !== null && previous >= 50) {
      worstRise = Math.max(worstRise, (sum - previous) / previous);
      worstFall = Math.max(worstFall, (previous - sum) / previous);
      steps++;
    }
    previous = sum;
    distance /= 2 ** (1 / 8);
  }
  assert.ok(steps > 20, `only ${steps} steps had 50 stars' worth on screen`);
  // Blocks hand over to their children over a band of size on screen, so a small step changes the stars on screen by
  // a few tens of percent at most (the first steps from afar are the largest, up to 70%: the field's blocks open one
  // level at a time and each child lists as many stars as its parent did); the magnitude limit
  // (static/starmagnitude.js) holds the count on screen.
  assert.ok(worstRise < 0.8, `one 9% zoom step raised the opacity on screen by ${(100 * worstRise).toFixed(0)}%`);
  assert.ok(worstFall < 0.4, `one 9% zoom step dropped the opacity on screen by ${(100 * worstFall).toFixed(0)}%`);
});

test("a star listed in a node and its children keeps its whole share through the hand-over", () => {
  // The brightest star of the field is in every node above it.
  const bright = { id: 0, x: 8000.5, y: 0.5, z: 0.5, luminosity: 1e6 };
  const address = G.sectorAddressAt(bright.x, bright.y, bright.z, EDGE);
  bright.chain = G.drillChainOf(address.shell, address.band, address.wedge).map((block) => C.formatNodeKey(block));
  const listing = lister([bright].concat(field(2000, 3)), 60);
  const look = [8000, 0, 0];
  let lowest = 1;
  let highest = 1;
  for (let distance = 8000; distance > 40; distance /= 1.05) {
    const camera = [look[0] + distance, look[1], look[2] + distance * 0.3];
    const view = drawn(listing, camera, look, 180);
    const mine = { entries: view.entries.filter(([id]) => id === 0), weigh: view.weigh };
    assert.ok(mine.entries.length >= 1);
    const sum = total(mine);
    lowest = Math.min(lowest, sum);
    highest = Math.max(highest, sum);
  }
  assert.ok(lowest > 0.999 && highest < 1.001, `the star's shares added up to between ${lowest} and ${highest}`);
});
