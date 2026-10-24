// tests/js/galaxynodecut.test.mjs -- static/galaxynodecut.js (MAP.261) against planetgen.galaxy.nodecut: the same cut
// of sector blocks and the same weights for the cameras tests/test_js_unit.py's fixtures hold, and the hand-over that
// waits for the children of a block to arrive.

import { createHash } from "node:crypto";
import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

const C = await importPage("galaxynodecut.js");
const cut = F.nodeCut;

test("the constants are the server's", () => {
  assert.deepEqual(C.HANDOVER_PX, cut.handover);
  assert.equal(C.FETCH_FROM_PX, cut.fetchFrom);
  assert.equal(C.MAX_NODES, cut.maxNodes);
  assert.equal(C.CONE_HALF_ANGLE_DEG, cut.halfAngle);
  assert.ok(Math.abs(C.DEFAULT_FOCAL_PX - cut.focal) < 1e-9);
});

for (const [index, one] of cut.cases.entries()) {
  test(`case ${index}: the cut holds the server's nodes in the server's order`, () => {
    const ours = C.nodeCut(one.camera, one.look, {
      edgePc: cut.edge, focalPx: cut.focal, bounds: one.bounds, maxNodes: cut.smallBudget,
    });
    assert.deepEqual(ours.keys, one.keys);
    assert.deepEqual([...ours.expanded].sort(), one.expanded);
  });

  test(`case ${index}: the full budget gives the server's cut too`, () => {
    const ours = C.nodeCut(one.camera, one.look, { edgePc: cut.edge, focalPx: cut.focal, bounds: one.bounds });
    assert.equal(ours.keys.length, one.fullCount);
    assert.equal(createHash("sha256").update(ours.keys.join(",")).digest("hex"), one.fullDigest);
  });

  test(`case ${index}: the weights are the server's`, () => {
    const weigh = C.weightCalculator(one.camera, cut.edge, cut.focal, new Set(one.expanded));
    for (const key of one.keys) {
      const got = weigh(C.parseNodeKey(key));
      const want = one.weights[key];
      assert.ok(Math.abs(got - want) <= 1e-9 * Math.max(1, Math.abs(want)), `${key}: ${got} vs ${want}`);
    }
  });
}

test("the weights of the blocks over a point add up to one", () => {
  const one = cut.cases[1];
  const keys = new Set(one.keys);
  const expanded = new Set(one.expanded);
  const weigh = C.weightCalculator(one.camera, cut.edge, cut.focal, expanded);
  // The chain from a root down along expanded blocks to the first shut one, through the first child kept at each step.
  let block = C.roots().find((root) => keys.has(C.formatNodeKey(root)));
  let sum = 0;
  while (block) {
    sum += weigh(block);
    const key = C.formatNodeKey(block);
    block = expanded.has(key) ? C.childrenOf(block).find((child) => keys.has(C.formatNodeKey(child))) : null;
  }
  assert.ok(Math.abs(sum - 1) < 1e-9, `the chain's weights add up to ${sum}`);
});

test("a block whose children are not all in hand does not hand over", () => {
  const one = cut.cases[1];
  const expanded = new Set(one.expanded);
  const shown = { keys: one.keys, expanded };
  const everything = C.settledExpanded(shown, () => true);
  assert.deepEqual([...everything].sort(), one.expanded);
  const nothing = C.settledExpanded(shown, () => false);
  const held = new Set(one.keys);
  for (const key of nothing) {
    const kept = C.childrenOf(C.parseNodeKey(key)).filter((child) => held.has(C.formatNodeKey(child)));
    assert.equal(kept.length, 0, `${key} hands over with children still to come`);
  }
  assert.ok(nothing.size < everything.size);
  const weigh = C.weightCalculator(one.camera, cut.edge, cut.focal, nothing);
  for (const root of C.roots()) {
    if (one.keys.includes(C.formatNodeKey(root))) assert.equal(weigh(root), 1, "a root with nothing under it stays whole");
  }
});

test("a sector gives nothing to children and a block over 140 pixels gives everything", () => {
  assert.equal(C.handover(1000, { m: 1, shell: 0, band: 0, wedge: 0 }), 0);
  assert.equal(C.handover(60, { m: 729, shell: 0, band: 0, wedge: 0 }), 0);
  assert.equal(C.handover(140, { m: 729, shell: 0, band: 0, wedge: 0 }), 1);
  assert.ok(Math.abs(C.handover(105, { m: 729, shell: 0, band: 0, wedge: 0 }) - 0.5) < 1e-12);
});
