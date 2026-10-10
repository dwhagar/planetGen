// tests/js/starmagnitude.test.mjs -- static/starmagnitude.js, the Galaxy
// Map's star visibility law (MAP.148): apparent magnitude, the opacity ramp,
// brightness by flux, and the histogram limit that holds about N stars on
// screen.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const M = await importPage("starmagnitude.js");

test("a Sun is magnitude 4.83 and the apparent magnitude follows the inverse square law", () => {
  assert.equal(M.absoluteMagnitude(1), 4.83);
  assert.ok(Math.abs(M.absoluteMagnitude(100) - (4.83 - 5)) < 1e-12, "100 L_sun is 5 magnitudes brighter");
  assert.equal(M.absoluteMagnitude(0), M.UNKNOWN_MAGNITUDE);
  // At 10 pc the two agree; ten times farther is 5 magnitudes dimmer.
  assert.equal(M.apparentMagnitude(4.83, 10), 4.83);
  assert.ok(Math.abs(M.apparentMagnitude(4.83, 100) - 9.83) < 1e-12);
});

test("the table in the design: a Sun is visible out to 22 pc at m_lim 6.5 and 1.1 kpc at 15", () => {
  const seen = (absolute, limit) => 10 ** ((limit - absolute) / 5 + 1);
  assert.ok(Math.abs(seen(4.83, 6.5) - 21.6) < 0.5);
  assert.ok(Math.abs(seen(4.83, 15) - 1075) < 20);
  // The ramp ends at the limit: share 0 there, 1 a ramp inside it.
  const edge = M.apparentMagnitude(4.83, seen(4.83, 12));
  assert.ok(Math.abs(M.magnitudeShare(edge, 12)) < 1e-9);
  assert.equal(M.magnitudeShare(12 - M.RAMP, 12), 1);
  assert.equal(M.magnitudeShare(5, 12), 1);
  assert.equal(M.magnitudeShare(14, 12), 0);
  assert.ok(Math.abs(M.magnitudeShare(12 - M.RAMP / 2, 12) - 0.5) < 1e-12);
});

test("brightness by flux is capped for a bright star and floors at the faint end of the ramp", () => {
  const limit = 12;
  const bright = M.fluxGain(-10, limit);
  assert.ok(bright > 0.999 && bright <= 1, String(bright));
  const faint = M.fluxGain(limit, limit);
  assert.ok(faint >= M.FLUX_GAIN_MIN && faint < 0.62, String(faint));
  let last = 0;
  for (let m = limit + 2; m >= limit - 8; m -= 0.5) {
    const g = M.fluxGain(m, limit);
    assert.ok(g >= last - 1e-12, "brighter stars never get less gain");
    last = g;
  }
});

function field(count, seed) {
  // Absolute magnitudes of a made-up population, dim stars most common.
  let x = seed;
  const next = () => {
    x = (x * 1664525 + 1013904223) % 4294967296;
    return x / 4294967296;
  };
  return Array.from({ length: count }, () => -4 + 16 * Math.pow(next(), 0.5));
}

test("the histogram limit holds about the target number of stars on screen", () => {
  const mags = field(60000, 7).map((absolute, i) => M.apparentMagnitude(absolute, 5 + (i % 997)));
  const limit = M.histogramLimit(mags, null, 20000);
  const shown = mags.reduce((sum, m) => sum + M.magnitudeShare(m, limit), 0);
  assert.ok(shown > 0.8 * 20000 && shown < 1.25 * 20000, `${shown} stars' worth shown at limit ${limit}`);
  // A smaller target gives a brighter limit.
  assert.ok(M.histogramLimit(mags, null, 8000) < limit);
});

test("with no more stars than the target, everything shows whole", () => {
  const mags = [3, 5, 9];
  const limit = M.histogramLimit(mags, null, 20000);
  mags.forEach((m) => assert.equal(M.magnitudeShare(m, limit), 1));
  assert.equal(M.histogramLimit([], null, 10), M.RAMP);
});

test("weights count a half-shown star as half", () => {
  const mags = Array.from({ length: 1000 }, (_, i) => i / 100);
  const full = M.histogramLimit(mags, null, 500);
  const half = M.histogramLimit(mags, mags.map(() => 0.5), 500);
  assert.ok(half > full, "half-weighted stars let more of them in under the same target");
  const none = M.histogramLimit(mags, mags.map(() => 0), 500);
  assert.equal(none, M.RAMP, "stars with no weight are not counted");
});

test("the calibration limit keeps the dimmest listed star whole at the focus distance", () => {
  const dim = M.dimmestShown([1, 2, 3, 10, 11], [1, 1, 1, 1, 0.2], 1);
  assert.equal(dim, 10, "a barely shown star (weight 0.2) doesn't count");
  assert.equal(M.dimmestShown([], null), null);
  assert.equal(M.dimmestShown([M.UNKNOWN_MAGNITUDE], null), null);
  const limit = M.calibrationLimit(10, 100);
  assert.equal(M.magnitudeShare(M.apparentMagnitude(10, 100), limit), 1);
  // A star twice as far dims; one 10 times farther is gone.
  assert.ok(M.magnitudeShare(M.apparentMagnitude(10, 200), limit) < 1);
  assert.equal(M.magnitudeShare(M.apparentMagnitude(10, 1000), limit), 0);
  assert.equal(M.calibrationLimit(null, 100), -Infinity);
  assert.equal(M.limitingMagnitude(8, 12), 12);
  assert.equal(M.limitingMagnitude(14, 12), 14);
});

test("a landmark stays lit at any distance and the limit eases rather than jumps", () => {
  const far = M.apparentMagnitude(M.PHENOMENON_MAGNITUDE, 1e5);
  assert.equal(M.magnitudeShare(far, 12), 1);
  assert.equal(M.approachLimit(10, 14, 0), 10);
  const later = M.approachLimit(10, 14, 0.3);
  assert.ok(later > 12.5 && later < 12.6, String(later));
  assert.ok(Math.abs(M.approachLimit(10, 14, 30) - 14) < 1e-6);
  assert.equal(M.approachLimit(Infinity, 14, 1), 14);
});

test("the shader's twin carries the same constants", () => {
  const glsl = M.STAR_MAGNITUDE_GLSL;
  assert.match(glsl, /float magnitudeShare\(float m, float limit\)/);
  assert.match(glsl, /float fluxGain\(float m, float limit\)/);
  assert.match(glsl, /float apparentMagnitude\(float M, float distancePc\)/);
  [M.RAMP, M.FLUX_HALF_BELOW, M.FLUX_GAIN_MIN, 1 - M.FLUX_GAIN_MIN].forEach((value) => {
    assert.ok(glsl.includes(value.toFixed(4)), `constant ${value} in the GLSL`);
  });
});
