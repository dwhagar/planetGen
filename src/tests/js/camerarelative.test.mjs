// tests/js/camerarelative.test.mjs -- static/camerarelative.js (MAP.240): buffers built as offsets from an origin of their own, so a
// shape 262 kpc from the galaxy's centre is as exact as one at the centre.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const C = await importPage("camerarelative.js");

const FAR = 262143.456789;   // parsecs, at the edge of the reach

test("a point far from the origin keeps its place as an offset, where one float32 would be 0.016 pc off", () => {
  const point = [FAR, -FAR / 3, FAR / 7];
  const built = C.relativeBuffer(point, point.map((v) => v - 12.345678901));
  // The shape's own float32 holds the absolute value to no better than the spacing there.
  const absolute = Math.abs(Math.fround(point[0]) - point[0]);
  assert.ok(absolute > 1e-3, `absolute float32 error ${absolute}`);
  for (let axis = 0; axis < 3; axis++) {
    const back = built.origin[axis] + built.positions[axis];
    assert.ok(Math.abs(back - point[axis]) < 1e-6, `axis ${axis}: ${back} vs ${point[axis]}`);
  }
});

test("the origin is the middle of the points' box and the offsets are about the origin", () => {
  const flat = [100, 0, 0, 300, 20, -40];
  assert.deepEqual(C.centreOf(flat), [200, 10, -20]);
  assert.deepEqual(Array.from(C.relativeBuffer(flat).positions), [-100, -10, 20, 100, 10, -20]);
  assert.deepEqual(C.centreOf([]), [0, 0, 0]);
});

test("scattered points are split into cubes, each with an origin inside it and every offset within the cube", () => {
  const flat = [];
  for (let i = 0; i < 200; i++) flat.push(FAR - i * 37.5, 1000 * Math.sin(i), -FAR + (i % 7) * 900);
  const chunks = C.groupByCell(flat, 1000);
  assert.ok(chunks.length > 3);
  const seen = chunks.flatMap((chunk) => chunk.indices).sort((a, b) => a - b);
  assert.deepEqual(seen, Array.from({ length: 200 }, (_, n) => n), "every point is in exactly one chunk");
  chunks.forEach((chunk) => {
    assert.ok(C.extentOf(chunk.positions) <= 500 + 1e-3, `offsets stay within half a cube: ${C.extentOf(chunk.positions)}`);
    chunk.indices.forEach((at, n) => {
      for (let axis = 0; axis < 3; axis++) {
        assert.ok(Math.abs(chunk.origin[axis] + chunk.positions[3 * n + axis] - flat[3 * at + axis]) < 1e-3);
      }
    });
  });
});

test("a long course becomes short pieces that join end to end, each exact about its own middle", () => {
  const stops = [[FAR, 0, 0], [FAR - 30000, 5000, 100], [FAR - 30000, 5000, 20100]];
  const pieces = C.linePieces(stops, 4000);
  assert.ok(pieces.length >= 8 + 5);
  pieces.forEach((piece) => assert.ok(C.extentOf(piece.positions) <= 2000 + 1e-3));
  const ends = (piece) => [0, 1].map((n) => [0, 1, 2].map((axis) => piece.origin[axis] + piece.positions[3 * n + axis]));
  assert.ok(Math.abs(ends(pieces[0])[0][0] - FAR) < 1e-3, "starts at the first stop");
  for (let n = 0; n + 1 < pieces.length; n++) {
    ends(pieces[n])[1].forEach((v, axis) => assert.ok(Math.abs(v - ends(pieces[n + 1])[0][axis]) < 1e-3), `piece ${n} meets ${n + 1}`);
  }
  assert.equal(C.linePieces([[0, 0, 0]], 10).length, 0);
});
