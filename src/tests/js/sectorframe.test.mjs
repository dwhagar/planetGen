// tests/js/sectorframe.test.mjs -- static/sectorframe.js (MAP.327): the walls of a sector's cell and where its neighbours are.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const F = await importPage("sectorframe.js");
const P = await importPage("galaxyprisms.js");

const EDGE = 4;
const addressKey = (a) => [a.shell, a.band, a.wedge].join(".");
const norm = (p) => Math.hypot(p[0], p[1], p[2]);

// A sector well inside the disc, away from the galaxy's poles and centre.
const SECTOR = { shell: 200, band: 0, wedge: 5 };
const bounds = P.sectorBounds(SECTOR.shell, SECTOR.band, SECTOR.wedge, EDGE);

test("a cell inside the galaxy has all six walls, each with a neighbour just beyond it", () => {
  const walls = F.wallCentres(bounds, 0.25 * EDGE);
  assert.deepEqual(walls.map((w) => w.dir).sort(), ["east", "in", "north", "out", "south", "west"]);
  const here = addressKey(SECTOR);
  const seen = new Set([here]);
  walls.forEach((wall) => {
    const there = P.sectorAddressAt(wall.beyond[0], wall.beyond[1], wall.beyond[2], EDGE);
    assert.notEqual(addressKey(there), here, wall.dir + " wall's neighbour is another sector");
    seen.add(addressKey(there));
  });
  // The wall centres themselves are on the cell's own boundary: the point just inside is this sector.
  assert.equal(seen.size, 7, "six different neighbours");
});

test("the outer and inner walls are at the cell's radii and the others between them", () => {
  const walls = Object.fromEntries(F.wallCentres(bounds, 1).map((w) => [w.dir, w]));
  assert.ok(Math.abs(norm(walls.out.point) - bounds.r1) < 1e-6);
  assert.ok(Math.abs(norm(walls.in.point) - bounds.r0) < 1e-6);
  ["north", "south", "east", "west"].forEach((dir) => {
    const r = norm(walls[dir].point);
    assert.ok(r > bounds.r0 && r < bounds.r1, dir);
  });
  assert.ok(norm(walls.out.beyond) > norm(walls.out.point) && norm(walls.in.beyond) < norm(walls.in.point));
});

test("the wall on the galaxy's own boundary has no icon: no inner wall at the centre, no wall past a pole, no longitude wall all the way round", () => {
  const core = { r0: 0, r1: 4, c0: 0, c1: Math.PI, t0: 0, t1: 2 * Math.PI };
  assert.deepEqual(F.wallCentres(core, 1).map((w) => w.dir), ["out"]);
  const cap = { r0: 100, r1: 104, c0: 0, c1: 0.1, t0: 0, t1: 2 * Math.PI };
  assert.deepEqual(F.wallCentres(cap, 1).map((w) => w.dir).sort(), ["in", "out", "south"]);
  const polar = P.sectorBounds(200, P.sectorAddressAt(0.01, 0, 800, EDGE).band, 0, EDGE);
  assert.ok(F.wallCentres(polar, 1).every((w) => w.dir !== "north"), "no north wall at the north pole");
});

test("the walls move with the cell: the same bounds always give the same points", () => {
  const a = F.wallCentres(bounds, 1);
  const b = F.wallCentres(bounds, 1);
  assert.deepEqual(a, b);
});

test("every wall has an arrow and words", () => {
  ["out", "in", "north", "south", "west", "east"].forEach((dir) => {
    assert.ok(F.WALL_LABELS[dir].arrow && F.WALL_LABELS[dir].words, dir);
  });
});
