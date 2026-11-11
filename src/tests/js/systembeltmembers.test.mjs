// tests/js/systembeltmembers.test.mjs -- the belt members the 3D system view draws one by one (MAP.414): static/systemview3d.js
// `buildSystemScene` over a scene whose belts carry `members` (the scene's `systemscene._belt_members`), and the radius they are drawn at.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage, installDom } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

installDom("http://localhost/system");
const View = await importPage("systemview3d.js");
const S = await importPage("systemscale.js");

const AU_KM = 149597870.7;

function sceneWithBelt(members) {
  const scene = JSON.parse(JSON.stringify(F.orbitPositions.scene));
  for (const star of scene.stars) Object.assign(star, { kind: "star", radius_km: 695700 });
  for (const planet of scene.planets) {
    Object.assign(planet, { kind: "planet", radius_km: planet.ref === "planet:1" ? 6371 : 69911 });
    for (const moon of planet.moons) Object.assign(moon, { kind: "moon", radius_km: 1737 });
  }
  for (const comet of scene.comets) Object.assign(comet, { kind: "comet", radius_km: 5 });
  scene.system = scene.system || { id: 1, ref: "system:1", name: "Test", heliopause_au: null };
  scene.belts = [{ ref: "belt:1", id: 1, around: "barycenter", inner_km: 2.5 * AU_KM, outer_km: 3.2 * AU_KM, members: members }];
  return scene;
}

function member(designation, radiusKm, lumpy, at) {
  return {
    designation: designation, rank: 1, radius_km: radiusKm, position_km: at,
    look: { v: 1, kind: "belt_member", body: "solid", seed: 12345 + radiusKm, color: "#8a8478", quality: "composition", radiusKm: radiusKm, roundness: lumpy ? 0.3 : 1.5, lumpy: lumpy },
  };
}

test("a belt member is drawn larger when it is larger and smaller than any moon, in the compressed scale", () => {
  const layout = S.createLayout(sceneWithBelt([]), S.MODE_COMPRESSED);
  const small = layout.radiusFor("belt_member", 2);
  const large = layout.radiusFor("belt_member", 470);
  assert.ok(small > 0 && large > small);
  assert.ok(large < layout.radiusFor("moon", 0), "smaller than the smallest moon");
});

test("in true scale a member is its radius in scene units", () => {
  const layout = S.createLayout(sceneWithBelt([]), S.MODE_TRUE);
  assert.equal(layout.radiusFor("belt_member", 470), 470 / S.UNIT_KM);
});

test("each member of a belt is a mesh at the place its orbit gives, hung from what the belt goes round, and it is not an entry", () => {
  const members = [member("A.0001", 470, false, [2.8 * AU_KM, 0, 0]), member("A.0002", 40, true, [0, 2.9 * AU_KM, 1e6])];
  const scene = sceneWithBelt(members);
  const layout = S.createLayout(scene, S.MODE_COMPRESSED);
  const built = View.buildSystemScene(scene, layout, { glow: false });
  const bare = View.buildSystemScene(sceneWithBelt([]), S.createLayout(sceneWithBelt([]), S.MODE_COMPRESSED), { glow: false });
  assert.equal(built.entries.length, bare.entries.length, "members are not picked, labelled or dotted");
  for (const one of members) {
    const mesh = built.objects["belt_member:" + one.designation];
    assert.ok(mesh, one.designation);
    assert.equal(mesh.userData.member, one.designation);
    assert.equal(mesh.scale.x, layout.radiusFor("belt_member", one.radius_km));
    const local = layout.place("barycenter", one.position_km);
    assert.deepEqual([mesh.position.x, mesh.position.y, mesh.position.z].map((v) => Math.fround(v)), local.map((v) => Math.fround(v)));
    const place = built.byRef["belt_member:" + one.designation];
    assert.ok(Math.abs(place.x - mesh.position.x) < 1e-6 && Math.abs(place.radius - mesh.scale.x) < 1e-12);
  }
  const lump = built.objects["belt_member:A.0002"].geometry, sphere = built.objects["belt_member:A.0001"].geometry;
  assert.notEqual(lump.type, undefined);
  assert.ok(lump.attributes.position.count > 0 && sphere.attributes.position.count > 0);
  built.dispose();
  bare.dispose();
});

test("a lumpy member is not the sphere a round one is", () => {
  const scene = sceneWithBelt([member("A.0001", 470, false, [2.8 * AU_KM, 0, 0]), member("A.0002", 40, true, [0, 2.9 * AU_KM, 0])]);
  const built = View.buildSystemScene(scene, S.createLayout(scene, S.MODE_COMPRESSED), { glow: false });
  const radii = (name) => {
    const attribute = built.objects[name].geometry.attributes.position;
    const out = new Set();
    for (let n = 0; n < attribute.count; n += 1) out.add(Math.hypot(attribute.getX(n), attribute.getY(n), attribute.getZ(n)).toFixed(4));
    return out.size;
  };
  assert.equal(radii("belt_member:A.0001"), 1, "round: every point at one radius");
  assert.ok(radii("belt_member:A.0002") > 5, "lumpy: points at many radii");
  built.dispose();
});

test("a belt with no members draws as before", () => {
  const scene = sceneWithBelt(undefined);
  delete scene.belts[0].members;
  const built = View.buildSystemScene(scene, S.createLayout(scene, S.MODE_COMPRESSED), { glow: false });
  assert.ok(Object.keys(built.objects).every((ref) => !ref.startsWith("belt_member:")));
  built.dispose();
});

test("the same members draw the same meshes twice: places, sizes and shapes", () => {
  const members = [member("A.0001", 470, false, [2.8 * AU_KM, 0, 0]), member("A.0002", 40, true, [0, 2.9 * AU_KM, 1e6]), member("A.0003", 9, true, [-2.6 * AU_KM, 1e7, 0])];
  const digest = () => {
    const scene = sceneWithBelt(members);
    const built = View.buildSystemScene(scene, S.createLayout(scene, S.MODE_COMPRESSED), { glow: false });
    const out = members.map((one) => {
      const mesh = built.objects["belt_member:" + one.designation];
      return [mesh.position.x, mesh.position.y, mesh.position.z, mesh.scale.x, mesh.geometry.attributes.position.count,
        Array.from(mesh.geometry.attributes.position.array.slice(0, 30))];
    });
    built.dispose();
    return JSON.stringify(out);
  };
  assert.equal(digest(), digest());
});
