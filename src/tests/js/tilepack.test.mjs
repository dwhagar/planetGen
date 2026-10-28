// tests/js/tilepack.test.mjs -- static/tilepack.js, the reader of MAP.159's
// packed Galaxy Map nodes, against what planetgen/web/lib/tilepack.py
// writes (the packed answer and the nodes Python itself reads from it come
// from tests/test_js_unit.py's fixtures).

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

const P = await importPage("tilepack.js");

const bytes = () => new Uint8Array(Buffer.from(F.nodePack.body, "base64"));

function close(actual, expected, label) {
  const scale = Math.max(Math.abs(expected), 1e-12);
  assert.ok(Math.abs(actual - expected) <= 1e-9 * scale, `${label}: ${actual} vs ${expected}`);
}

test("a packed answer reads to the nodes the server's own reader finds", () => {
  const { header, nodes, packed } = P.unpackResponse(bytes());
  assert.deepEqual(header, F.nodePack.header);
  assert.deepEqual(Object.keys(nodes), Object.keys(F.nodePack.nodes));
  for (const key of Object.keys(nodes)) {
    const want = F.nodePack.nodes[key];
    const got = nodes[key];
    assert.deepEqual(got.clouds, want.clouds);
    assert.deepEqual(got.points, want.points);
    assert.deepEqual(got.incomplete, want.incomplete);
    for (const section of ["stars", "generated"]) {
      assert.equal(got[section].length, want[section].length, `${key} ${section}`);
      got[section].forEach((star, i) => {
        const expected = want[section][i];
        assert.deepEqual(Object.keys(star).sort(), Object.keys(expected).sort(), `${key} ${section} ${i} fields`);
        for (const field of Object.keys(expected)) {
          if (typeof expected[field] === "number") close(star[field], expected[field], `${key} ${section} ${i} ${field}`);
          else assert.equal(star[field], expected[field], `${key} ${section} ${i} ${field}`);
        }
      });
    }
    assert.equal(packed.get(key).length, F.nodePack.sizes[key]);
  }
});

test("a star keeps its place in the list and an id past 32 bits stays exact", () => {
  const { nodes } = P.unpackResponse(bytes());
  assert.deepEqual(nodes["3.1.2.3"].stars.map((star) => star.id), [7, 2 ** 40 + 3, 9]);
  assert.equal(nodes["3.1.2.3"].stars[1].star_type, "M");
  assert.equal(nodes["3.1.2.3"].stars[1].charted, false);
  assert.equal(nodes["3.1.2.3"].stars[0].charted, true);
  // MAP.252: a satellite galaxy's star carries its tier, the others none.
  assert.equal(nodes["3.1.2.3"].stars[2].tier, 2);
  assert.equal(nodes["3.1.2.3"].stars[2].charted, false);
  assert.equal("tier" in nodes["3.1.2.3"].stars[0], false);
});

test("the names of named stars come back on their stars, by id (MAP.298)", () => {
  const { nodes } = P.unpackResponse(bytes());
  assert.equal(nodes["3.1.2.3"].stars[0].name, "Named Giant");
  assert.equal(nodes["3.1.2.3"].generated[0].name, "Printed Name");
  assert.equal("name" in nodes["3.1.2.3"].stars[1], false);
  assert.equal("names" in nodes["3.1.2.4"], false);
});

test("printed object IDs and mixed ids come back as they went in", () => {
  const { nodes } = P.unpackResponse(bytes());
  assert.deepEqual(nodes["3.1.2.3"].generated.map((star) => star.id), ["FE81000A2B-0000005-000", "FE81000A2B-0000006-000"]);
  assert.deepEqual(nodes["3.1.2.6"].stars.map((star) => star.id), [3, "odd"]);
});

test("a star without a temperature or radius comes out without", () => {
  const { nodes } = P.unpackResponse(bytes());
  const star = nodes["3.1.2.3"].stars[1];
  assert.equal("temperature_k" in star, false);
  assert.equal("radius_sol" in star, false);
});

test("a stored node reads on its own", () => {
  const { packed } = P.unpackResponse(bytes());
  const node = P.unpackNode(packed.get("3.1.2.3"));
  assert.equal(node.stars.length, 3);
  assert.equal(node.generated.length, 2);
  assert.equal(P.unpackNode(packed.get("3.1.2.5")).stars.length, 0);
});

test("bytes that are not a whole node are refused", () => {
  const { packed } = P.unpackResponse(bytes());
  const whole = packed.get("3.1.2.3");
  assert.throws(() => P.unpackNode(whole.slice(0, whole.length - 3)));
  assert.throws(() => P.unpackNode(new Uint8Array([123, 34, 115, 116, 97, 114, 115, 34, 58, 91])));
  assert.throws(() => P.unpackResponse(bytes().slice(0, 40)));
});
