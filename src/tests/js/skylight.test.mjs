// tests/js/skylight.test.mjs -- static/skylight.js (VIEW.32), the arithmetic that draws one star of the sky behind a system:
// apparent magnitude to opacity, size and colour, the star colour against photometry.py's `star_colour`, and the GLSL twins.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

const S = await importPage("skylight.js");
const M = await importPage("starmagnitude.js");

const near = (a, b, eps = 1e-9) => Math.abs(a - b) <= eps;

test("the peak is 1 at magnitude 4.5 and brighter and 0.2 at 6.5", () => {
  assert.equal(S.skyPeak(4.5), 1);
  assert.equal(S.skyPeak(-1.46), 1);
  assert.ok(near(S.skyPeak(6.5), Math.pow(10, -0.8)));
  assert.ok(near(S.skyPeak(6.5), 0.1585, 1e-4));
});

test("the opacity fades over the last ramp to the limit and is nothing past it", () => {
  const limit = S.SKY_LIMIT;
  assert.equal(S.skyOpacity(limit + 0.01, limit), 0);
  assert.equal(S.skyOpacity(limit, limit), 0);
  assert.ok(near(S.skyOpacity(limit - M.RAMP, limit), S.skyPeak(limit - M.RAMP)));
  const inside = S.skyOpacity(limit - M.RAMP / 2, limit);
  assert.ok(inside > 0 && inside < S.skyPeak(limit - M.RAMP / 2), "the middle of the ramp is half way (smoothstep)");
  assert.ok(near(inside, 0.5 * S.skyPeak(limit - M.RAMP / 2)));
  assert.equal(S.skyOpacity(0, limit), 1);
});

test("a dimmer star is never drawn stronger than a brighter one", () => {
  let last = Infinity;
  for (let m = -2; m <= 9; m += 0.1) {
    const o = S.skyOpacity(m, S.SKY_LIMIT);
    assert.ok(o <= last + 1e-12, `m=${m}`);
    last = o;
  }
});

test("the phone's limit is lower, and a given limit wins", () => {
  assert.ok(S.SKY_LIMIT_PHONE < S.SKY_LIMIT);
  assert.equal(S.skyLimit(false), S.SKY_LIMIT);
  assert.equal(S.skyLimit(true), S.SKY_LIMIT_PHONE);
  assert.equal(S.skyLimit(true, 5), 5);
  assert.equal(S.skyLimit(false, undefined), S.SKY_LIMIT);
});

test("Sirius is a wide point and a star at the limit is small", () => {
  assert.ok(near(S.skySigmaPx(-1.46), 0.7 + 0.45 * (4.5 + 1.46)));
  assert.ok(near(S.skySigmaPx(-1.46), 3.382, 1e-3));
  assert.equal(S.skySizePx(-1.46), S.MAX_SIZE_PX);
  assert.ok(near(S.skySizePx(7.5), 2.8), "a faint star: four sigmas of 0.7 px");
  assert.ok(S.skySizePx(7.5) > S.MIN_SIZE_PX, "and never under the Galaxy Map's smallest core");
  let last = Infinity;
  for (let m = -3; m <= 8; m += 0.25) {
    const px = S.skySizePx(m);
    assert.ok(px <= last + 1e-12 && px >= S.MIN_SIZE_PX && px <= S.MAX_SIZE_PX, `m=${m}`);
    last = px;
  }
});

test("a faint star shows less of its colour, a bright one all of it", () => {
  assert.equal(S.skySaturation(1.5), 1);
  assert.equal(S.skySaturation(-1), 1);
  assert.equal(S.skySaturation(9.5), S.SAT_MIN);
  assert.ok(near(S.skySaturation(5.5), 0.35 + 0.65 * 0));
  const orange = [1, 0.6, 0.2];
  S.skyTint(orange, 0).forEach((c, k) => assert.ok(near(c, orange[k])));
  const dim = S.skyTint(orange, 9.5);
  for (let k = 0; k < 3; k++) assert.ok(dim[k] >= orange[k] && dim[k] <= 1);
  assert.ok(near(dim[2], 1 - 0.8 * 0.35));
  assert.deepEqual(S.skyTint([1, 1, 1], 7), [1, 1, 1]);
});

test("the colour of a star matches photometry.py's star_colour at every temperature", () => {
  for (const row of F.skyColour) {
    const got = S.starColourExact(row.temperatureK);
    for (let k = 0; k < 3; k++) assert.ok(near(got[k], row.linear[k], 1e-6), `${row.temperatureK} K channel ${k}: ${got[k]} vs ${row.linear[k]}`);
  }
});

test("the table lookup stays within 0.01 of the exact colour and is red for a cool star, blue for a hot one", () => {
  for (const row of F.skyColour) {
    const got = S.starColour(row.temperatureK);
    for (let k = 0; k < 3; k++) assert.ok(near(got[k], row.linear[k], 0.01), `${row.temperatureK} K channel ${k}: ${got[k]} vs ${row.linear[k]}`);
  }
  const cool = S.starColour(3000);
  const hot = S.starColour(30000);
  assert.equal(cool[0], 1);
  assert.ok(cool[2] < cool[0]);
  assert.equal(hot[2], 1);
  assert.ok(hot[0] < hot[2]);
  assert.deepEqual(S.starColour(0), [1, 1, 1]);
  assert.deepEqual(S.starColour(NaN), [1, 1, 1]);
});

test("a star put together: a Sun at 10 pc, one at 100 pc, and one with no luminosity", () => {
  const near10 = S.skyStar(1, 5772, 10);
  assert.equal(near10.m, 4.83);
  assert.ok(near10.opacity > 0.7 && near10.opacity <= 1);
  const far = S.skyStar(1, 5772, 100);
  assert.ok(near(far.m, 9.83));
  assert.equal(far.opacity, 0);
  const none = S.skyStar(0, 5772, 10);
  assert.equal(none.opacity, 0);
  assert.equal(S.skyStar(1, 5772, 10, 3).opacity, 0, "a lower limit drops it");
});

test("the GLSL carries the same constants and every function the sky shader calls", () => {
  const g = S.SKY_LIGHT_GLSL;
  for (const name of ["apparentMagnitude", "magnitudeShare", "fluxGain", "skyPeak", "skyOpacity", "skySizePx", "skySaturation", "skyTint"]) {
    assert.ok(g.includes(name + "("), name);
  }
  assert.ok(g.startsWith(M.STAR_MAGNITUDE_GLSL), "built on the Galaxy Map's GLSL, not a copy");
  for (const v of [S.PEAK_MAGNITUDE, S.SIGMA_PX, S.SIGMA_PER_MAG_PX, S.SIZE_SIGMAS, S.MIN_SIZE_PX, S.MAX_SIZE_PX, S.SAT_MIN, S.SAT_SLOPE, S.SAT_REF_MAG, S.SAT_SPAN_MAG]) {
    assert.ok(g.includes(v.toFixed(4)), `constant ${v}`);
  }
});

// A line-by-line evaluation of the GLSL's size and saturation text, so the shader's text and the JS cannot drift apart.
test("the GLSL's functions give the JS functions' values", () => {
  const body = (name) => {
    const m = new RegExp("float " + name + "\\(float m\\) \\{([\\s\\S]*?)\\n\\}").exec(S.SKY_LIGHT_GLSL);
    assert.ok(m, name);
    return m[1];
  };
  const run = (name, m) => {
    const js = body(name)
      .replace(/float (\w+) =/g, "const $1 =")
      .replace(/return /, "return ")
      .replace(/clamp\(/g, "__clamp(")
      .replace(/max\(/g, "Math.max(")
      .replace(/(\d)\.0+(?!\d)/g, "$1");
    return new Function("m", "__clamp", js)(m, (x, lo, hi) => Math.min(hi, Math.max(lo, x)));
  };
  for (let m = -3; m <= 10; m += 0.5) {
    assert.ok(near(run("skySizePx", m), S.skySizePx(m), 1e-3), `size m=${m}`);
    assert.ok(near(run("skySaturation", m), S.skySaturation(m), 1e-3), `sat m=${m}`);
  }
});
