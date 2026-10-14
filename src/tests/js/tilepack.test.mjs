// tests/js/tilepack.test.mjs -- static/tilepack.js, the reader of MAP.159's
// packed Galaxy Map tiles, against what planetgen/web/lib/tilepack.py
// writes (the packed answer and the tiles Python itself reads from it come
// from tests/test_js_unit.py's fixtures).

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

const P = await importPage("tilepack.js");

const bytes = () => new Uint8Array(Buffer.from(F.tilePack.body, "base64"));

function close(actual, expected, label) {
  const scale = Math.max(Math.abs(expected), 1e-12);
  assert.ok(Math.abs(actual - expected) <= 1e-9 * scale, `${label}: ${actual} vs ${expected}`);
}

test("a packed answer reads to the tiles the server's own reader finds", () => {
  const { header, tiles, packed } = P.unpackResponse(bytes());
  assert.deepEqual(header, F.tilePack.header);
  assert.deepEqual(Object.keys(tiles), Object.keys(F.tilePack.tiles));
  for (const key of Object.keys(tiles)) {
    const want = F.tilePack.tiles[key];
    const got = tiles[key];
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
    assert.equal(packed.get(key).length, F.tilePack.sizes[key]);
  }
});

test("a star keeps its place in the list and an id past 32 bits stays exact", () => {
  const { tiles } = P.unpackResponse(bytes());
  assert.deepEqual(tiles["3/1/2/3"].stars.map((star) => star.id), [7, 2 ** 40 + 3, 9]);
  assert.equal(tiles["3/1/2/3"].stars[1].star_type, "M");
  assert.equal(tiles["3/1/2/3"].stars[1].charted, false);
  assert.equal(tiles["3/1/2/3"].stars[0].charted, true);
});

test("printed object IDs and mixed ids come back as they went in", () => {
  const { tiles } = P.unpackResponse(bytes());
  assert.deepEqual(tiles["3/1/2/3"].generated.map((star) => star.id), ["FE81000A2B-0000005-000", "FE81000A2B-0000006-000"]);
  assert.deepEqual(tiles["3/1/2/6"].stars.map((star) => star.id), [3, "odd"]);
});

test("a star without a temperature or radius comes out without", () => {
  const { tiles } = P.unpackResponse(bytes());
  const star = tiles["3/1/2/3"].stars[1];
  assert.equal("temperature_k" in star, false);
  assert.equal("radius_sol" in star, false);
});

test("a stored tile reads on its own", () => {
  const { packed } = P.unpackResponse(bytes());
  const tile = P.unpackTile(packed.get("3/1/2/3"));
  assert.equal(tile.stars.length, 3);
  assert.equal(tile.generated.length, 2);
  assert.equal(P.unpackTile(packed.get("3/1/2/5")).stars.length, 0);
});

test("bytes that are not a whole tile are refused", () => {
  const { packed } = P.unpackResponse(bytes());
  const whole = packed.get("3/1/2/3");
  assert.throws(() => P.unpackTile(whole.slice(0, whole.length - 3)));
  assert.throws(() => P.unpackTile(new Uint8Array([123, 34, 115, 116, 97, 114, 115, 34, 58, 91])));
  assert.throws(() => P.unpackResponse(bytes().slice(0, 40)));
});
