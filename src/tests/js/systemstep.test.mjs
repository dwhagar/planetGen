// tests/js/systemstep.test.mjs -- static/systemstep.js (MAP.121): the nearest system in a direction.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const S = await importPage("systemstep.js");

// A point 1,000 pc from the galaxy's centre on the x axis, on the plane: out is +x, east +y, north +z (colatitude shrinks
// toward the pole), south -z.
const HERE = { x: 1000, y: 0, z: 0 };
const near = (name, dx, dy, dz) => ({ name, x: HERE.x + dx, y: HERE.y + dy, z: HERE.z + dz });
const NEIGHBOURS = [near("out", 2, 0, 0), near("in", -3, 0, 0), near("east", 0, 4, 0), near("west", 0, -5, 0), near("north", 0, 0, 6), near("south", 0, 0, -7)];

test("each direction is a unit vector in the galaxy's frame at the point", () => {
  const expected = { out: [1, 0, 0], in: [-1, 0, 0], east: [0, 1, 0], west: [0, -1, 0], north: [0, 0, 1], south: [0, 0, -1] };
  Object.entries(expected).forEach(([direction, vector]) => {
    const got = S.directionAt([HERE.x, HERE.y, HERE.z], direction);
    vector.forEach((value, i) => assert.ok(Math.abs(got[i] - value) < 1e-12, direction + " " + got));
  });
  const off = S.directionAt([300, 400, 500], "east");
  assert.ok(Math.abs(Math.hypot(off[0], off[1], off[2]) - 1) < 1e-12);
  assert.ok(Math.abs(off[0] * 300 + off[1] * 400 + off[2] * 500) < 1e-9, "east is level with the point's own radius");
});

test("a step goes to the system lying that way, not to the nearest of all", () => {
  NEIGHBOURS.forEach((n) => assert.equal(S.nearestInDirection(HERE, NEIGHBOURS, n.name), n, n.name));
  // Inside the cone (55 degrees off east) and nearer: it wins over the one straight east.
  assert.equal(S.nearestInDirection(HERE, [near("close", 0.5, 0.5, 0.5), near("east", 0, 4, 0)], "east").name, "close");
});

test("of several in the cone the nearest wins, and what lies outside the cone is not offered", () => {
  const field = [near("far east", 0, 9, 0), near("near east", 1, 3, 0), near("behind", 0, -1, 0)];
  assert.equal(S.nearestInDirection(HERE, field, "east").name, "near east");
  assert.equal(S.nearestInDirection(HERE, [near("behind", 0, -1, 0)], "east"), null);
});

test("the system stood at is skipped, and the poles and the centre have no sideways", () => {
  const me = { x: HERE.x, y: HERE.y, z: HERE.z, name: "me" };
  assert.equal(S.nearestInDirection(HERE, [me], "out"), null, "zero distance is the system itself");
  assert.equal(S.nearestInDirection(HERE, [me, near("out", 2, 0, 0)], "out", (c) => c === me).name, "out");
  assert.equal(S.directionAt([0, 0, 500], "east"), null);
  assert.equal(S.directionAt([0, 0, 0], "out"), null);
  assert.equal(S.directionAt([0, 0, 500], "out")[2], 1);
  assert.equal(S.directionAt([1, 1, 1], "sideways"), null);
});
