// tests/js/skypick.test.mjs -- static/skypick.js (VIEW.37): which star of the sky a pointer is on, and what is said of it.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const P = await importPage("skypick.js");

// A camera looking down -z with a 90 degree vertical field on a 200 x 100 canvas: the galaxy's axes already the camera's.
const VIEW = { matrix: [1, 0, 0, 0, 1, 0, 0, 0, 1], p00: 0.5, p11: 1, width: 200, height: 100 };

function skyOf(stars) {
  const flat = (key) => Float32Array.from(stars.flatMap((s) => s[key]));
  return {
    count: stars.length, observerVelocityKms: [0, 0, 0], offsetPc: flat("offset"), velocityKms: flat("velocity"),
    magnitude: Float32Array.from(stars.map((s) => s.M)), temperatureK: Float32Array.from(stars.map((s) => s.T)),
    kind: Uint8Array.from(stars.map((s) => s.kind || 0)), index: Float64Array.from(stars.map((s, i) => s.index ?? i)),
  };
}

const star = (offset, M, extra) => Object.assign({ offset, velocity: [0, 0, 0], M, T: 5800 }, extra);

test("the star nearest the pointer within reach is the pick, and none beyond it", () => {
  // Straight ahead (-z) at 10 pc: the middle of the canvas, (100, 50). Ten pc to the right at 10 pc: x_ndc = 0.5 * 1 = 0.5, so 150 px.
  const sky = skyOf([star([0, 0, -10], 4.83), star([10, 0, -10], 4.83), star([0, 0, 10], 0)]);
  const middle = P.pickSkyStar(sky, 0, VIEW, { x: 103, y: 52 }, { limit: 7.5 });
  assert.equal(middle.place, 0);
  assert.ok(Math.abs(middle.distancePc - 10) < 1e-6 && middle.tauYears > 0);
  assert.equal(P.pickSkyStar(sky, 0, VIEW, { x: 150, y: 50 }, { limit: 7.5 }).place, 1);
  assert.equal(P.pickSkyStar(sky, 0, VIEW, { x: 20, y: 90 }, { limit: 7.5 }), null, "nothing there");
});

test("a star behind the camera or dimmer than the limit is not picked, and the brighter of two as near wins", () => {
  const sky = skyOf([star([0, 0, 10], 0), star([0, 0, -10], 20), star([0.0001, 0, -10], 4.83), star([0, 0, -10], 0)]);
  const found = P.pickSkyStar(sky, 0, VIEW, { x: 100, y: 50 }, { limit: 7.5 });
  assert.equal(found.place, 3, "the one behind is out, the faint one is out, the brightest of the rest wins");
});

test("the light time moves the picked star: a star moving across the line of sight is picked where it was", () => {
  const moving = skyOf([star([0, 0, -10], 4.83, { velocity: [20000, 0, 0] })]);
  const years = 0;
  const now = P.pickSkyStar(moving, years, VIEW, { x: 100, y: 50 }, { limit: 7.5, reachPx: 30 });
  assert.ok(now, "found near the middle");
  assert.ok(now.tauYears > 30 && now.tauYears < 35, "10 pc is about 32.6 years of light");
});

test("the spectral class and the 'N years ago' words", () => {
  assert.deepEqual([40000, 12000, 8000, 6500, 5800, 4500, 3000, NaN].map(P.spectralClass), ["O", "B", "A", "F", "G", "K", "M", null]);
  assert.equal(P.agoText(0.2), "seen as it was less than a year ago");
  assert.equal(P.agoText(1), "seen as it was 1 year ago");
  assert.equal(P.agoText(4.4), "seen as it was 4 years ago");
  assert.equal(P.agoText(1234), "seen as it was 1,234 years ago");
});

test("a charted star is named and linked through its printed ID; a scatter star is told but not linked", () => {
  const charted = { 7: ["FE81-0000001-002", "Alpha Test"] };
  const pick = { place: 0, index: 7, kind: 0, distancePc: 1.5, tauYears: 4.89, magnitude: -0.3, temperatureK: 5800, screenPx: 1 };
  const one = P.describeSkyStar(pick, charted, "/object/star:{id}");
  assert.equal(one.title, "Alpha Test");
  assert.equal(one.href, "/object/star:FE81-0000001-002");
  assert.ok(one.lines.join(";").includes("seen as it was 5 years ago") && one.lines.join(";").includes("class G"));
  assert.deepEqual(one.fields.map((f) => f[0]), ["Class", "Distance", "Brightness", "Light"]);
  const other = P.describeSkyStar({ ...pick, index: 9, kind: 1 }, charted, "/object/star:{id}");
  assert.equal(other.title, "Star not yet charted");
  assert.equal(other.href, null);
  assert.ok(P.hoverText(other).startsWith("Star not yet charted: "));
  assert.equal(P.describeSkyStar(pick, undefined, "/object/star:{id}").href, null, "an older list has no names");
});

test("the stars on the screen are the brightest in front of the camera and inside the canvas", () => {
  const sky = skyOf([star([0, 0, -10], 4.83), star([100, 0, -10], 0), star([0, 0, 10], -2), star([0, 0, -20], 1)]);
  const found = P.skyStarsOnScreen(sky, 0, VIEW, { limit: 7.5 });
  assert.deepEqual(found.map((f) => f.place), [3, 0]);
  assert.ok(Math.abs(found[0].x - 100) < 1e-6 && Math.abs(found[0].y - 50) < 1e-6);
});
