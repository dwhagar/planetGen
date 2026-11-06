// tests/js/systemcourse.test.mjs -- static/systemcourse.js, a NAV course inside a system as the 3D view draws it (NAV.22).

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const { drawnLine, drawnCourse, PIECES_PER_SEGMENT } = await import(new URL("systemcourse.js", STATIC_URL).href);
const Scale = await import(new URL("systemscale.js", STATIC_URL).href);

const AU = 1.495978707e8;
const scene = {
  stars: [{ ref: "star:1", kind: "star", radius_km: 6.96e5 }],
  planets: [
    { ref: "planet:1", kind: "planet", radius_km: 6.4e3, orbit: { around: "barycenter", distance_km: 1 * AU }, moons: [] },
    { ref: "planet:2", kind: "planet", radius_km: 6.0e4, orbit: { around: "barycenter", distance_km: 30 * AU }, moons: [] },
  ],
  comets: [], belts: [], tree: null,
};

test("at true scale a straight leg stays one straight piece and is in million kilometers", () => {
  const layout = Scale.createLayout(scene, Scale.MODE_TRUE);
  const line = drawnLine([[0, 0, 0], [3e6, 0, 0]], layout);
  assert.deepEqual(line, [[0, 0, 0], [3, 0, 0]]);
});

test("in the compressed scale a leg is cut into pieces that follow the log scale out from the middle", () => {
  const layout = Scale.createLayout(scene, Scale.MODE_COMPRESSED);
  const line = drawnLine([[0, 0, 0], [30 * AU, 0, 0]], layout);
  assert.equal(line.length, PIECES_PER_SEGMENT + 1);
  const radii = line.map((p) => Math.hypot(p[0], p[1], p[2]));
  assert.ok(radii.every((r, n) => n === 0 || r > radii[n - 1]), "it runs outward");
  assert.ok(Math.abs(radii[radii.length - 1] - layout.distance("barycenter", 30 * AU)) < 1e-9, "its end is where a body 30 AU out is drawn");
  assert.ok(radii[PIECES_PER_SEGMENT / 2] > 0.5 * radii[radii.length - 1], "the middle of the leg is past half way out: the scale is a log");
});

test("a course draws each leg with its ends, and where a moving destination is aimed", () => {
  const layout = Scale.createLayout(scene, Scale.MODE_COMPRESSED);
  const heliopause = 120 * AU;
  const course = { navUrl: "/nav", legs: [
    { kind: "into", label: "Into the system", from: "the heliopause", to: "Beta", points: [[-heliopause, 0, 0], [0, 0, 0], [0.3 * AU, 0, 0]], aimed: [0.31 * AU, 0, 0] },
    { kind: "out", label: "Out of the system", from: "Beta", to: "the heliopause", points: [[0.3 * AU, 0, 0]], aimed: null },
  ] };
  const drawn = drawnCourse(course, layout);
  assert.equal(drawn.length, 1, "a leg with fewer than two points is not drawn");
  assert.equal(drawn[0].kind, "into");
  assert.deepEqual(drawn[0].start, drawnLine([[-heliopause, 0, 0], [0, 0, 0]], layout)[0]);
  assert.ok(Math.abs(Math.hypot(...drawn[0].start) - layout.distance("barycenter", heliopause)) < 1e-9, "it starts on the heliopause sphere");
  assert.deepEqual(drawn[0].aimed, layout.place("barycenter", [0.31 * AU, 0, 0]));
  assert.deepEqual(drawnCourse(null, layout), []);
});
