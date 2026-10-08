// tests/js/nebulamesh.test.mjs -- static/nebulamesh.js, the nebula shape
// mesh the maps draw (MAP.103): the shape URL, the one-request-per-nebula
// loader and the scaling of a served mesh.

import { test } from "node:test";
import assert from "node:assert/strict";

import { createShapeLoader, flatFaces, scaledPositions, shapeUrl } from "../../html/static/nebulamesh.js";

const SHAPE = { vertices: [[1, 0, 0], [0, 0.5, 0], [0, 0, -1]], faces: [[0, 1, 2]] };

test("the shape URL fills in the id and the level of detail", () => {
  assert.equal(shapeUrl("/galaxy/nebula/{id}/shape", 7), "/galaxy/nebula/7/shape?lod=low");
  assert.equal(shapeUrl("/galaxy/nebula/{id}/shape", 7, "full"), "/galaxy/nebula/7/shape?lod=full");
});

test("a nebula's shape is asked for once per level of detail", async () => {
  const asked = [];
  const load = createShapeLoader("/s/{id}", (url) => { asked.push(url); return Promise.resolve(SHAPE); });
  await Promise.all([load(3), load(3), load(3, "full"), load(4)]);
  assert.deepEqual(asked, ["/s/3?lod=low", "/s/3?lod=full", "/s/4?lod=low"]);
});

test("a failed request is forgotten, so the next call asks again", async () => {
  let calls = 0;
  const load = createShapeLoader("/s/{id}", () => (++calls === 1 ? Promise.reject(new Error("down")) : Promise.resolve(SHAPE)));
  await assert.rejects(load(1));
  assert.equal(await load(1), SHAPE);
  assert.equal(calls, 2);
});

test("vertices scale from nebula-radius units to the map's units", () => {
  assert.deepEqual(Array.from(scaledPositions(SHAPE, 10)), [10, 0, 0, 0, 5, 0, 0, 0, -10]);
  assert.deepEqual(Array.from(flatFaces(SHAPE)), [0, 1, 2]);
});
