// tests/js/nodeview.test.mjs -- static/nodeview.js (MAP.249): a view is known from its nodes alone. The guard here is
// that the pass never reads a node's stars, so its cost follows the number of nodes (a view holds about a thousand)
// and not the number of stars, and that it still tells a changed view from an unchanged one.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const V = await importPage("nodeview.js");

// A node whose star arrays count every read, so a pass over the stars shows.
function node(starCount, reads) {
  const guard = (list) => new Proxy(list, { get(target, name) { reads.count += 1; return target[name]; } });
  return { stars: guard(new Array(starCount).fill(null).map((_, i) => ({ id: i }))), generated: guard([]), points: guard([]) };
}

function view(nodeCount, starsPerNode) {
  const reads = { count: 0 };
  const nodes = new Map();
  const keys = [];
  for (let i = 0; i < nodeCount; i += 1) {
    const key = `3.${i}.0.0`;
    keys.push(key);
    nodes.set(key, node(starsPerNode, reads));
  }
  return { keys, nodes, reads, getNode: (key) => nodes.get(key) };
}

test("equal views have equal signatures, in any order, and a replaced or missing node changes it", () => {
  const v = view(6, 3);
  const first = V.drawnSignature(v.keys, v.getNode);
  assert.deepEqual(first.missing, []);
  assert.equal(V.drawnSignature([...v.keys].reverse(), v.getNode).signature, first.signature);
  const replaced = new Map(v.nodes);
  replaced.set(v.keys[2], { stars: [], generated: [], points: [] });
  assert.notEqual(V.drawnSignature(v.keys, (key) => replaced.get(key)).signature, first.signature);
  const gone = V.drawnSignature(v.keys, (key) => (key === v.keys[1] ? undefined : v.nodes.get(key)));
  assert.deepEqual(gone.missing, [v.keys[1]]);
  assert.notEqual(gone.signature, first.signature);
});

test("a view of 1,200 nodes and 100,000 stars is signed without reading a star, in a few milliseconds", () => {
  const v = view(1200, 80);
  V.drawnSignature(v.keys, v.getNode);   // the revisions are numbered the first time
  const started = performance.now();
  const rounds = 20;
  for (let i = 0; i < rounds; i += 1) V.drawnSignature(v.keys, v.getNode);
  const each = (performance.now() - started) / rounds;
  assert.equal(v.reads.count, 0, "the signature read a node's stars");
  assert.ok(each < 8, `a pass over 1,200 nodes took ${each.toFixed(2)} ms`);
});
