// tests/js/systemscale.test.mjs -- static/systemscale.js (MAP.71): true
// scale and the compressed scale that keeps every body visible.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

const S = await import(new URL("systemscale.js", STATIC_URL).href);
const P = await import(new URL("orbitpositions.js", STATIC_URL).href);

const AU_KM = 149597870.7;
const scene = JSON.parse(JSON.stringify(F.orbitPositions.scene));
// Give the hand-built scene the fields the layout reads.
for (const star of scene.stars) Object.assign(star, { kind: "star", radius_km: 695700 });
for (const planet of scene.planets) {
  Object.assign(planet, { kind: "planet", radius_km: planet.ref === "planet:1" ? 6371 : 69911 });
  for (const moon of planet.moons) Object.assign(moon, { kind: "moon", radius_km: 1737 });
}
for (const comet of scene.comets) Object.assign(comet, { kind: "comet", radius_km: 5 });
scene.belts = [{ ref: "belt:1", around: "barycenter", inner_km: 2.5 * AU_KM, outer_km: 3.2 * AU_KM }];

test("true scale draws real kilometres over a million", () => {
  const layout = S.createLayout(scene, S.MODE_TRUE);
  assert.deepEqual(layout.place("barycenter", [1e9, 0, -2e9]), [1000, 0, -2000]);
  assert.equal(layout.radiusOf("star:1"), 0.6957);
  assert.match(layout.note, /True scale/);
});

test("compressed scale keeps order, separates the planets and enlarges bodies", () => {
  const layout = S.createLayout(scene, S.MODE_COMPRESSED);
  const near = Math.hypot(...layout.place("barycenter", [1.5 * AU_KM, 0, 0]));
  const far = Math.hypot(...layout.place("barycenter", [5 * AU_KM, 0, 0]));
  assert.ok(near > 0 && far > near, "outer orbits stay outside inner ones");
  assert.ok(far < 400, "even 5 AU fits in the view");
  assert.ok(layout.radiusOf("planet:1") > 6371 / S.UNIT_KM * 10, "bodies are larger than life");
  assert.deepEqual(layout.place("barycenter", [0, 0, 0]), [0, 0, 0]);
  assert.match(layout.note, /Compressed/);
});

test("a moon stays outside its planet's drawn size in both modes", () => {
  const moonRel = [0.002 * AU_KM, 0, 0];
  const compressed = S.createLayout(scene, S.MODE_COMPRESSED);
  assert.ok(Math.hypot(...compressed.place("planet:1", moonRel)) > compressed.radiusOf("planet:1") + compressed.radiusOf("moon:1"));
});

test("layoutPositions sums offsets up the parents", () => {
  const layout = S.createLayout(scene, S.MODE_COMPRESSED);
  const relative = P.relativeAt(scene, 0.4);
  const placed = S.layoutPositions(layout, relative);
  const planet = placed["planet:1"];
  const moon = placed["moon:1"];
  const offset = layout.place("planet:1", relative["moon:1"].rel);
  assert.ok(Math.hypot(moon[0] - planet[0] - offset[0], moon[1] - planet[1] - offset[1], moon[2] - planet[2] - offset[2]) < 1e-9);
});

test("orbit paths close up for planets and open out for a parabolic comet", () => {
  const circle = P.orbitPath(scene.planets[0].orbit, 64);
  assert.equal(circle.length, 65);
  assert.ok(Math.hypot(...circle[0].map((v, i) => v - circle[64][i])) < 1);
  assert.ok(circle.every((p) => Math.abs(Math.hypot(...p) - 1.5 * AU_KM) < 1e-3 * AU_KM));
  const comet = P.orbitPath(scene.comets[1].orbit, 64, 500);
  const reach = Math.max(...comet.map((p) => Math.hypot(...p))) / AU_KM;
  assert.ok(reach <= 500.01 && reach > 100, `a parabola out to ${reach} AU`);
  const ellipse = P.orbitPath(scene.comets[0].orbit, 64);
  assert.ok(Math.hypot(...ellipse[0]) / AU_KM > 29, "an ellipse starts at aphelion (true anomaly -180 degrees)");
});
