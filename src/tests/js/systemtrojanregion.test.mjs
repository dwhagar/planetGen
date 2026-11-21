// tests/js/systemtrojanregion.test.mjs -- a trojan region (GEN.1150, MAP.433) in the 3D system view: static/systemview3d.js `buildSystemScene` draws a belt row
// whose scene entry carries `region` as an arc about its giant's L4 or L5 (turned to the giant's longitude plus or minus 60 degrees plus the libration), not as a ring.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage, installDom } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

installDom("http://localhost/system");
const View = await importPage("systemview3d.js");
const S = await importPage("systemscale.js");

const AU_KM = 149597870.7;
const SECONDS_PER_YEAR = 365.25 * 86400;

function sceneWithRegions(regions) {
  const scene = JSON.parse(JSON.stringify(F.orbitPositions.scene));
  for (const star of scene.stars) Object.assign(star, { kind: "star", radius_km: 695700 });
  for (const planet of scene.planets) {
    Object.assign(planet, { kind: "planet", radius_km: planet.ref === "planet:1" ? 6371 : 69911 });
    for (const moon of planet.moons) Object.assign(moon, { kind: "moon", radius_km: 1737 });
  }
  for (const comet of scene.comets) Object.assign(comet, { kind: "comet", radius_km: 5 });
  scene.system = scene.system || { id: 1, ref: "system:1", name: "Test", heliopause_au: null };
  scene.belts = regions.map((region, index) => ({
    ref: "belt:" + (index + 1), id: index + 1, around: "barycenter", inner_km: 4.9 * AU_KM, outer_km: 5.1 * AU_KM, members: [],
    region: Object.assign({ host: "planet:2", point: "L4", amplitude_deg: 0, phase_deg: 0, libration_period_years: null, half_arc_deg: 8, mean_inclination_deg: 3 }, region),
  }));
  return scene;
}

function longitudes(points) {
  const attribute = points.geometry.attributes.position;
  const out = [];
  for (let n = 0; n < attribute.count; n += 1) out.push(Math.atan2(attribute.getY(n), attribute.getX(n)) * 180 / Math.PI);
  return out;
}

test("a region is an arc: every particle lies within its half width of the arc's center, none all the way round", () => {
  const scene = sceneWithRegions([{}]);
  const built = View.buildSystemScene(scene, S.createLayout(scene, S.MODE_COMPRESSED), { glow: false });
  const points = built.group.children.find((child) => child.userData.ref === "belt:1");
  const spread = longitudes(points);
  assert.ok(spread.every((deg) => Math.abs(deg) <= 8 + 1e-6), "inside the half width about longitude zero");
  assert.ok(Math.max(...spread) - Math.min(...spread) > 8, "spread across the arc");
  built.dispose();
});

test("the arc is turned to the giant's longitude plus 60 degrees for L4 and minus 60 for L5", () => {
  const scene = sceneWithRegions([{ point: "L4" }, { point: "L5" }]);
  const built = View.buildSystemScene(scene, S.createLayout(scene, S.MODE_COMPRESSED), { glow: false });
  const placed = built.update(0);
  const host = built.byRef["planet:2"];
  const giantDeg = Math.atan2(host.y, host.x) * 180 / Math.PI;
  const turned = (ref) => built.group.children.find((child) => child.userData.ref === ref).rotation.z * 180 / Math.PI;
  const wrap = (deg) => ((deg % 360) + 540) % 360 - 180;
  assert.ok(Math.abs(wrap(turned("belt:1") - giantDeg - 60)) < 1e-6);
  assert.ok(Math.abs(wrap(turned("belt:2") - giantDeg + 60)) < 1e-6);
  void placed;
  built.dispose();
});

test("the libration moves the arc about its point at the pace of the stored period and phase, counted from Unix time zero", () => {
  const region = { point: "L4", amplitude_deg: 20, phase_deg: 30, libration_period_years: 150 };
  const scene = sceneWithRegions([region]);
  const built = View.buildSystemScene(scene, S.createLayout(scene, S.MODE_COMPRESSED), { glow: false });
  const at = (years) => {
    built.update(years);
    const host = built.byRef["planet:2"];
    const giantDeg = Math.atan2(host.y, host.x) * 180 / Math.PI;
    const points = built.group.children.find((child) => child.userData.ref === "belt:1");
    return ((points.rotation.z * 180 / Math.PI - giantDeg + 540) % 360) - 180;
  };
  for (const years of [0, 12.5, 40]) {
    const t = scene.epoch_unix / SECONDS_PER_YEAR + years;
    const want = 60 + 20 * Math.sin(2 * Math.PI * t / 150 + 30 * Math.PI / 180);
    assert.ok(Math.abs(at(years) - want) < 1e-6, "at " + years + " years");
  }
  built.dispose();
});

test("a plain belt has no region and stays a ring", () => {
  const scene = sceneWithRegions([{}]);
  delete scene.belts[0].region;
  const built = View.buildSystemScene(scene, S.createLayout(scene, S.MODE_COMPRESSED), { glow: false });
  const points = built.group.children.find((child) => child.userData.ref === "belt:1");
  const spread = longitudes(points);
  assert.ok(Math.max(...spread) - Math.min(...spread) > 300, "all the way round");
  assert.equal(points.rotation.z, 0);
  built.dispose();
});
