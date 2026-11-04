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

test("a planet's path is the ellipse of its eccentricity, periapsis first (GEN.727)", () => {
  const orbit = { ...scene.planets[0].orbit, eccentricity: 0.4, arg_periapsis_deg: 90 };
  const path = P.orbitPath(orbit, 64);
  assert.equal(path.length, 65);
  const radii = path.map((p) => Math.hypot(...p) / AU_KM);
  assert.ok(Math.abs(Math.min(...radii) - 1.5 * 0.6) < 1e-9, "periapsis is a (1 - e)");
  assert.ok(Math.abs(Math.max(...radii) - 1.5 * 1.4) < 1e-9, "apoapsis is a (1 + e)");
  assert.ok(Math.abs(radii[0] - 1.5 * 0.6) < 1e-9, "the path starts at periapsis");
  assert.ok(Math.hypot(...path[0].map((v, i) => v - path[64][i])) < 1, "and closes");
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

// --- MAP.213: a tree of stars is drawn node by node ---------------------------------------------------------------------

function treeScene(name, withPlanets) {
  const tree = JSON.parse(JSON.stringify(F.orbitPositions.trees[name].scene));
  for (const star of tree.stars) Object.assign(star, { kind: "star", radius_km: 695700 });
  for (const planet of tree.planets) {
    Object.assign(planet, { kind: "planet", radius_km: 6371 });
    for (const moon of planet.moons) Object.assign(moon, { kind: "moon", radius_km: 1737 });
  }
  tree.belts = [];
  if (!withPlanets) tree.planets = [];
  return tree;
}

const apart = (placed, a, b) => Math.hypot(...placed[a].map((c, axis) => c - placed[b][axis]));

test("compressed, the stars of a triple clear each other at every time and the outer pair is wider than the inner (MAP.213)", () => {
  const triple = treeScene("triple", false);
  const layout = S.createLayout(triple, S.MODE_COMPRESSED);
  const radius = layout.radiusOf("star:1");
  for (const years of [0, 0.3, 7.5, 80, 400]) {
    const placed = S.layoutPositions(layout, P.relativeAt(triple, years));
    for (const [a, b] of [["star:1", "star:2"], ["star:1", "star:3"], ["star:2", "star:3"]]) {
      assert.ok(apart(placed, a, b) > 2 * radius, `${a} and ${b} overlap at ${years} y`);
    }
  }
  const placed = S.layoutPositions(layout, P.relativeAt(triple, 0));
  const centre = (a, b) => placed[a].map((c, axis) => (c + placed[b][axis]) / 2);
  const innerGap = apart(placed, "star:1", "star:2");
  const outerGap = Math.hypot(...centre("star:1", "star:2").map((c, axis) => c - placed["star:3"][axis]));
  assert.ok(outerGap > innerGap * 1.5, "the outer orbit is drawn wider than the inner one");
});

test("true scale draws a tree at its real distances (MAP.213)", () => {
  const triple = treeScene("triple", false);
  const layout = S.createLayout(triple, S.MODE_TRUE);
  const placed = S.layoutPositions(layout, P.relativeAt(triple, 0));
  const real = P.positionsAt(triple, 0);
  for (const ref of ["star:1", "star:2", "star:3"]) {
    real[ref].forEach((c, axis) => assert.ok(Math.abs(placed[ref][axis] - c / S.UNIT_KM) <= 1e-9 * Math.abs(c / S.UNIT_KM) + 1e-9));
  }
});

test("compressed, a wide pair with planets round its stars keeps the planets' orbits clear of the other star (MAP.213)", () => {
  const wide = treeScene("wide", true);
  wide.planets[0].orbit.around = "star:1";
  const layout = S.createLayout(wide, S.MODE_COMPRESSED);
  const placed = S.layoutPositions(layout, P.relativeAt(wide, 0));
  const reachOfPlanets = Math.max(...[0.001, 1, 10].map((au) => Math.hypot(...layout.place("star:1", [au * AU_KM, 0, 0]))));
  assert.ok(apart(placed, "star:1", "star:2") > reachOfPlanets, "the second star lies outside the first star's planets");
});

test("a node's ellipse keeps its shape: drawn distances are in proportion to the real ones (MAP.213)", () => {
  const wide = treeScene("wide", false);
  const layout = S.createLayout(wide, S.MODE_COMPRESSED);
  const near = Math.hypot(...layout.place("node:0", [0.4 * 30 * AU_KM, 0, 0]));
  const far = Math.hypot(...layout.place("node:0", [1.6 * 30 * AU_KM, 0, 0]));
  assert.ok(Math.abs(far / near - 4) < 1e-9);
});

test("compressed, a planet round the inner pair of a triple is drawn outside the pair's own orbit (a close pair below a further star)", () => {
  const triple = treeScene("triple", false);
  triple.planets = [{
    ref: "planet:1", kind: "planet", radius_km: 6371, moons: [],
    orbit: { around: "node:0", distance_km: 3 * AU_KM, period_years: 3, inclination_deg: 0, ascending_node_deg: 0,
      position_km: [3 * AU_KM, 0, 0], velocity_kms: [0, 17, 0], epoch_unix: triple.epoch_unix },
  }];
  const layout = S.createLayout(triple, S.MODE_COMPRESSED);
  const inner = triple.tree.nodes.find((node) => node.node_no === 0);
  const reach = Math.hypot(...layout.place("node:0", [inner.separation_km * (1 + inner.eccentricity), 0, 0]));
  assert.ok(Math.hypot(...layout.place("node:0", [3 * AU_KM, 0, 0])) > reach, "the planet lies beyond the stars of the pair");
  const placed = S.layoutPositions(layout, P.relativeAt(triple, 0));
  const away = Math.hypot(...placed["planet:1"].map((c, axis) => c - placed["node:0"][axis]));
  assert.ok(Math.abs(away - Math.hypot(...layout.place("node:0", [3 * AU_KM, 0, 0]))) < 1e-9, "placed about the pair's centre");
});
