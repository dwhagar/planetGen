// tests/js/starpoint.test.mjs -- static/starpoint.js (MAP.208): a star or a
// phenomenon is a point of light until the camera is about a sector away and
// grows smoothly from there, on ONE near-sector constant the Galaxy Map and the
// sector scene share, so the hand-off between the two shows no jump in size.

import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";

import { installDom, importPage, STATIC_URL } from "./fakedom.mjs";

const P = await importPage("starpoint.js");

const source = (name) => readFileSync(new URL(name, STATIC_URL), "utf8");

test("a star is a point at far and mid distances and fully grown near a sector", () => {
  const full = 24;
  // The whole-galaxy view: a point POINT_PX across, whatever the star's size.
  assert.equal(P.pointSizePx(full, 20000), P.POINT_PX);
  assert.equal(P.pointSizePx(full, P.POINT_FAR_PC), P.POINT_PX);
  assert.equal(P.pointSizePx(4, 20000), P.POINT_PX);
  assert.ok(P.POINT_PX >= 1 && P.POINT_PX <= 2, "one to two pixels");
  // Mid distance (a few sectors away): still a point, not yet a disc.
  const mid = P.pointSizePx(full, P.POINT_FAR_PC / 2);
  assert.ok(mid < full / 2, `mid ${mid}`);
  // Near a sector: the full size, and nearer stays the same (never sized by distance).
  assert.equal(P.pointSizePx(full, P.NEAR_SECTOR_PC), full);
  assert.equal(P.pointSizePx(full, P.NEAR_SECTOR_PC / 10), full);
  assert.equal(P.pointSizePx(full, 0), full);
});

test("the size rises smoothly and monotonically over the last zoom steps, never snapping", () => {
  let last = -1;
  let biggestStep = 0;
  for (let d = P.POINT_FAR_PC * 1.2; d >= P.NEAR_SECTOR_PC * 0.8; d *= 0.97) {
    const size = P.pointSizePx(30, d);
    assert.ok(size >= last - 1e-9, `grows as the camera nears (${d})`);
    if (last >= 0) biggestStep = Math.max(biggestStep, size - last);
    last = size;
  }
  assert.ok(biggestStep < 2.5, `largest step ${biggestStep} px for a 3% closing`);
  assert.equal(P.POINT_FAR_PC, P.NEAR_SECTOR_PC * P.ZOOM_STEP ** P.POINT_RAMP_STEPS);
});

test("an arriving star (MAP.155) stays a point until its zoom share rises", () => {
  assert.equal(P.pointSizePx(20, 1, 0), P.POINT_PX);
  assert.equal(P.pointSizePx(20, 1, 1), 20);
  assert.ok(P.pointSizePx(20, 1, 0.5) > P.POINT_PX && P.pointSizePx(20, 1, 0.5) < 20);
});

test("a point is all core and the core share eases back to the full look", () => {
  assert.equal(P.pointCoreShare(4, 20, 20000), 1);
  assert.ok(Math.abs(P.pointCoreShare(4, 20, 1) - 4 / 20) < 1e-12);
  const between = P.pointCoreShare(4, 20, P.POINT_FAR_PC / 3);
  assert.ok(between < 1 && between > 4 / 20);
});

test("brightness alone tells far stars apart: dimmer for a fainter star", () => {
  assert.ok(P.pointBrightness(0) < P.pointBrightness(0.5));
  assert.ok(P.pointBrightness(0.5) < P.pointBrightness(1));
  assert.equal(P.pointBrightness(1), 1);
  assert.equal(P.pointBrightness(-3), P.pointBrightness(0));
});

test("the hand-off: both scenes size a light from the same constant and function", () => {
  const galaxy = source("galaxymap3d.js");
  const sector = source("sectorscene.js");
  [galaxy, sector].forEach((text) => {
    assert.match(text, /starpoint\.js/, "imports the shared module");
    assert.match(text, /STAR_POINT_GLSL/, "uses the shared pointGrow shader function");
    assert.match(text, /pointGrow\(/, "calls it");
    assert.match(text, /POINT_PX/, "uses the shared point width");
    assert.doesNotMatch(text, /NEAR_SECTOR_PC\s*=/, "defines no near-sector constant of its own");
  });
  // The shader twin is built from the module's constants.
  assert.ok(P.STAR_POINT_GLSL.includes(P.NEAR_SECTOR_PC.toFixed(4)));
  assert.ok(P.STAR_POINT_GLSL.includes(Math.log(P.POINT_FAR_PC / P.NEAR_SECTOR_PC).toFixed(6)));
  // At the hand-off distance both draw the full size the sector scene's light carries.
  assert.equal(P.pointSizePx(21, P.NEAR_SECTOR_PC), 21);
});

test("the shader twin matches the JS function", () => {
  // Evaluate the GLSL by hand: same expression, same constants.
  const near = P.NEAR_SECTOR_PC;
  const span = Math.log(P.POINT_FAR_PC / near);
  for (const d of [0.5, 5, 16, 30, 60, 100, 144, 500, 1e5]) {
    const t = Math.min(1, Math.max(0, Math.log(Math.max(d, 1e-6) / near) / span));
    assert.ok(Math.abs((1 - t * t * (3 - 2 * t)) - P.pointGrow(d)) < 1e-12, `d=${d}`);
  }
});

test("a phenomenon is as bright as a star of its mass, never brighter than the brightest star", () => {
  assert.equal(P.phenomenonLuminositySol(5000, 1e9), 5000);
  assert.equal(P.phenomenonLuminositySol(5e9, 4000), 4000, "capped at the brightest star");
  assert.equal(P.phenomenonLuminositySol(5e9, 0), 5e9, "no star: no cap");
  assert.equal(P.phenomenonLuminositySol(0, 0), 1, "no luminosity on record: a one solar mass star");
  assert.equal(P.phenomenonLuminositySol(undefined, 7), 1);
});

test("Charted only dims a star or a scattered phenomenon outside the filled sectors, nothing inside", () => {
  assert.equal(P.isUncharted({ id: 1 }), true, "a pre-placed star in an uncharted sector");
  assert.equal(P.isUncharted({ id: 1, charted: true }), false);
  assert.equal(P.isUncharted({ id: 1, generated: true }), false, "a generated star is in a filled sector");
  assert.equal(P.isUncharted({ phenomenon: true, scattered: true, charted: false }), true, "an unbuilt scattered black hole");
  assert.equal(P.isUncharted({ phenomenon: true, scattered: true }), true, "a point from a tile that predates the field");
  assert.equal(P.isUncharted({ phenomenon: true, scattered: true, charted: true }), false, "its sector is built");
  assert.equal(P.isUncharted({ phenomenon: true, type: "black_hole" }), false, "a built phenomenon lies in a filled sector");
});

test("Charted only dims in the cut order: satellite stars, halo stars, then two ramps, the brightest and most massive last (MAP.252)", () => {
  const sat = P.unchartedKeep({ id: 1, tier: 2 }, 1);
  const halo = P.unchartedKeep({ id: 2, tier: 1 }, 1);
  const dimStar = P.unchartedKeep({ id: 3 }, 0);
  const brightStar = P.unchartedKeep({ id: 4 }, 1);
  assert.ok(sat < halo && halo < dimStar && dimStar < brightStar, [sat, halo, dimStar, brightStar].join(" < "));
  assert.equal(P.unchartedKeep({ id: 5, charted: true, tier: 2 }, 0), 1, "a charted star keeps all of its light");
  assert.equal(P.unchartedKeep({ id: 6, generated: true }, 0), 1);
  // The two ramps run together: a small phenomenon is cut with the small dim stars, a large one with the bright.
  const small = P.unchartedKeep({ phenomenon: true, scattered: true }, P.massRampShare(1));
  const large = P.unchartedKeep({ phenomenon: true, scattered: true }, P.massRampShare(1e8));
  assert.equal(small, dimStar);
  assert.equal(large, brightStar);
  assert.ok(P.massRampShare(1e4) > P.massRampShare(10) && P.massRampShare(10) > P.massRampShare(1));
  assert.equal(P.massRampShare(undefined), 0.5, "no mass sits in the middle");
  assert.equal(P.unchartedKeep({ id: 7 }, 7), brightStar, "a share past the ramp is clamped");
});

test("Charted only dims the farther from the focus more, and the glsl twin agrees", () => {
  assert.equal(P.unchartedFarShare(10, 100), 0, "nearer than the focus");
  assert.equal(P.unchartedFarShare(100, 100), 0);
  assert.ok(P.unchartedFarShare(200, 100) > 0 && P.unchartedFarShare(200, 100) < P.unchartedFarShare(300, 100));
  assert.equal(P.unchartedFarShare(400, 100), 1);
  assert.equal(P.unchartedFarShare(4000, 100), 1);
  assert.equal(P.unchartedFarShare(50, 0), 0, "no focus, no rule");
  assert.match(P.UNCHARTED_FAR_GLSL, /float unchartedFar\(float distancePc, float focusPc\)/);
  assert.match(P.UNCHARTED_FAR_GLSL, new RegExp(Math.log(P.UNCHARTED_FAR_RATIO).toFixed(6).replace(".", "\\.")));
  assert.ok(P.UNCHARTED_FAR_KEEP > 0 && P.UNCHARTED_FAR_KEEP < 1);
});

test("the Galaxy Map gives every object a keep share and the shader applies it with the distance rule", () => {
  const map = source("galaxymap3d.js");
  assert.match(map, /keeps\[i\] = unchartedKeep\(star, star\.phenomenon \? massRampShare\(star\.mass_solar\) : t\)/);
  assert.match(map, /attribute float starKeep;/);
  assert.match(map, /mix\(1\.0, unchartedFarKeep, unchartedFar\(dist, focusPc\)\)/);
  assert.match(map, /uniforms\.focusPc\.value = nearUniforms\.nearD\.value/);
});

test("the sector scene's far brightness follows the star's own and never exceeds it", async () => {
  installDom("http://localhost/sector/1");
  const SS = await import(new URL("sectorscene.js", STATIC_URL).href);
  assert.equal(SS.lightFarBright({ bright: 0.75 }), P.pointBrightness(0));
  assert.equal(SS.lightFarBright({ bright: 1 }), 1);
  assert.ok(SS.lightFarBright({ bright: 0.9 }) > SS.lightFarBright({ bright: 0.8 }));
  assert.equal(SS.lightFarBright({ bright: 0.2 }), 0.2, "a rogue planet stays dim");
});
