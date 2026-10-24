// tests/js/sectorreach.test.mjs -- static/sectorreach.js (MAP.152): which sector's scene the camera
// wants up, and how much of it shows. The spec is docs/design/interface/fly-through-view-distance.md 4.7: a
// sector's own scene fades in within 1.5 edges of its centre, is whole by 1 edge and is dropped past 3.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const R = await importPage("sectorreach.js");
const EDGE = 4;

test("the distances are the design's", () => {
  assert.equal(R.FADE_IN_EDGES, 1.5);
  assert.equal(R.FULL_EDGES, 1);
  assert.equal(R.UNLOAD_EDGES, 3);
});

test("the share is 0 from 1.5 edges out, 1 from 1 edge in, and rises smoothly between", () => {
  assert.equal(R.sectorShare(1.5 * EDGE, EDGE), 0);
  assert.equal(R.sectorShare(5 * EDGE, EDGE), 0);
  assert.equal(R.sectorShare(EDGE, EDGE), 1);
  assert.equal(R.sectorShare(0, EDGE), 1);
  let last = 0;
  for (let d = 1.5; d >= 1; d -= 0.05) {
    const share = R.sectorShare(d * EDGE, EDGE);
    assert.ok(share >= last - 1e-12, "never falls as the camera comes nearer");
    last = share;
  }
  assert.ok(Math.abs(R.sectorShare(1.25 * EDGE, EDGE) - 0.5) < 1e-9, "halfway at the middle of the band");
  assert.equal(R.sectorShare(1, 0), 0, "no edge, no scene");
});

const near = (x) => ({ key: "near", center: [x, 0, 0] });
const far = (x) => ({ key: "far", center: [x, 0, 0] });

test("nothing is wanted until the camera is within 1.5 edges of a sector", () => {
  assert.equal(R.chooseSector([0, 0, 0], [near(1.6 * EDGE)], EDGE, null), null);
  const found = R.chooseSector([0, 0, 0], [near(1.4 * EDGE)], EDGE, null);
  assert.equal(found.key, "near");
  assert.ok(found.share > 0 && found.share < 0.2);
});

test("the nearest sector in range is the one chosen", () => {
  const found = R.chooseSector([0, 0, 0], [far(1.4 * EDGE), near(0.9 * EDGE)], EDGE, null);
  assert.equal(found.key, "near");
  assert.equal(found.share, 1);
});

test("a scene up stays between 1.5 and 3 edges (hysteresis) with no share, and goes past 3", () => {
  const kept = R.chooseSector([0, 0, 0], [near(2.5 * EDGE)], EDGE, "near");
  assert.equal(kept.key, "near");
  assert.equal(kept.share, 0);
  assert.equal(R.chooseSector([0, 0, 0], [near(3.01 * EDGE)], EDGE, "near"), null);
  // ... and the same camera does not bring it back in from outside.
  assert.equal(R.chooseSector([0, 0, 0], [near(2.5 * EDGE)], EDGE, null), null);
});

test("a camera hovering at the edge of the range does not flip the scene on and off", () => {
  let current = null;
  const states = [];
  for (const d of [1.51, 1.49, 1.51, 1.7, 1.49, 2.9, 1.52, 3.1, 1.52]) {
    const found = R.chooseSector([0, 0, 0], [near(d * EDGE)], EDGE, current);
    current = found ? found.key : null;
    states.push(current);
  }
  assert.deepEqual(states, [null, "near", "near", "near", "near", "near", "near", null, null]);
});

test("one scene is up at a time: a sector takes over only when clearly nearer", () => {
  const a = { key: "a", center: [-EDGE, 0, 0] };
  const b = { key: "b", center: [EDGE, 0, 0] };
  // Equal distances: the one up stays.
  assert.equal(R.chooseSector([0, 0, 0], [a, b], EDGE, "a").key, "a");
  assert.equal(R.chooseSector([0, 0, 0], [a, b], EDGE, "b").key, "b");
  // The camera has moved toward b: b is under 85% of a's distance.
  assert.equal(R.chooseSector([0.4 * EDGE, 0, 0], [a, b], EDGE, "a").key, "b");
  // Slightly nearer is not enough.
  assert.equal(R.chooseSector([0.05 * EDGE, 0, 0], [a, b], EDGE, "a").key, "a");
});

test("a scene that is up for a sector no longer in view is let go", () => {
  assert.equal(R.chooseSector([0, 0, 0], [], EDGE, "gone"), null);
  const other = R.chooseSector([0, 0, 0], [{ key: "other", center: [EDGE, 0, 0] }], EDGE, "gone");
  assert.equal(other.key, "other");
});
