// tests/js/lumpyview.test.mjs -- static/lumpyview.js (MAP.363): the lump of a small body on the object page's 2D canvas.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const V = await importPage("lumpyview.js");

const LONG = { v: 1, kind: "moon", body: "solid", seed: 31337, color: "#8a8478", radiusKm: 30, roundness: 0.05 };
const SIZE = 160;

function draw(look, angle, tiltDeg = 0) {
  const mesh = V.lumpMesh(look);
  const covered = new Uint8Array(SIZE * SIZE);
  const depth = new Float32Array(SIZE * SIZE);
  const calls = [];
  V.rasterise(mesh, { angle: angle, tiltDeg: tiltDeg }, SIZE, SIZE / 2, SIZE / 2, (SIZE / 2) / mesh.extent * 0.95, covered, depth, (pixel, lon, lat, nx, ny, nz) => calls.push([pixel, lon, lat, nx, ny, nz]));
  return { mesh, covered, calls };
}

function box(covered) {
  let minX = SIZE, maxX = -1, minY = SIZE, maxY = -1, area = 0;
  for (let y = 0; y < SIZE; y++) for (let x = 0; x < SIZE; x++) {
    if (!covered[y * SIZE + x]) continue;
    area += 1;
    minX = Math.min(minX, x); maxX = Math.max(maxX, x); minY = Math.min(minY, y); maxY = Math.max(maxY, y);
  }
  return { width: maxX - minX + 1, height: maxY - minY + 1, area: area };
}

test("the grid has a point for every corner and outward unit normals", () => {
  const mesh = V.lumpMesh(LONG);
  assert.equal(mesh.count, (V.GRID_WIDTH + 1) * (V.GRID_HEIGHT + 1));
  assert.equal(mesh.indices.length, V.GRID_WIDTH * V.GRID_HEIGHT * 6);
  for (let n = 0; n < mesh.count; n++) {
    const a = n * 3;
    assert.ok(Math.abs(Math.hypot(mesh.normals[a], mesh.normals[a + 1], mesh.normals[a + 2]) - 1) < 1e-4);
    assert.ok(mesh.normals[a] * mesh.directions[a] + mesh.normals[a + 1] * mesh.directions[a + 1] + mesh.normals[a + 2] * mesh.directions[a + 2] > 0);
  }
  assert.ok(mesh.extent > 1);
});

test("the lump fills a solid shape inside the picture and every pixel is painted once with a unit normal and a longitude in range", () => {
  const { covered, calls } = draw(LONG, 0);
  const seen = covered.reduce((a, b) => a + b, 0);
  assert.ok(seen > 2000 && seen < SIZE * SIZE);
  for (const [pixel, lon, lat, nx, ny, nz] of calls) {
    assert.ok(covered[pixel] === 1 && Math.abs(lon) <= Math.PI + 1e-6 && Math.abs(lat) <= Math.PI / 2 + 1e-6);
    assert.ok(Math.abs(Math.hypot(nx, ny, nz) - 1) < 1e-4);
  }
});

test("an elongated lump is wide when its long axis crosses the view and narrow when it points at the viewer", () => {
  const side = box(draw(LONG, 0).covered);
  const end = box(draw(LONG, Math.PI / 2).covered);
  assert.ok(side.width > end.width * 1.1, `${side.width} against ${end.width}`);
  assert.ok(Math.abs(draw(LONG, 0).mesh.extent - draw(LONG, 1).mesh.extent) < 1e-9);
});

test("a turn of a whole circle gives the same picture, and the same record the same lump", () => {
  const a = draw(LONG, 0.4).covered;
  const b = draw(LONG, 0.4 + 2 * Math.PI).covered;
  assert.deepEqual(Array.from(a), Array.from(b));
  assert.deepEqual(Array.from(V.lumpMesh(LONG).positions), Array.from(V.lumpMesh({ ...LONG }).positions));
  assert.notDeepEqual(Array.from(V.lumpMesh(LONG).positions), Array.from(V.lumpMesh({ ...LONG, seed: 5 }).positions));
});

test("the side of the lump facing the viewer faces the light that comes from there; a tilted body still fills the picture", () => {
  const { calls } = draw(LONG, 0);
  assert.ok(calls.filter(([, , , , , nz]) => nz > 0.3).length > calls.length * 0.3);
  const tilted = box(draw(LONG, 0, 40).covered);
  assert.ok(tilted.area > 1000);
});

test("the texture longitude at the middle of the picture is the one facing the viewer, and moves east as the body turns", () => {
  const near = (angle) => {
    const { calls } = draw({ ...LONG, roundness: 0.9 }, angle);
    const middle = (SIZE / 2) * SIZE + SIZE / 2;
    return calls.filter(([pixel]) => pixel === middle).pop()[1];
  };
  assert.ok(Math.abs(near(0)) < 0.25);
  assert.ok(near(0.5) < near(0) - 0.2 || near(0.5) > near(0) + 6);       // the longitude shown is the body's less how far it has turned
});
