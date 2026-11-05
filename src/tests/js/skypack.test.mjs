// tests/js/skypack.test.mjs -- static/skypack.js, the reader of VIEW.29's packed sky list, against what
// planetgen/web/lib/skypack.py writes (the packed list and the columns Python itself reads from it come from
// tests/test_js_unit.py's `sky_pack` fixture).

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

const P = await importPage("skypack.js");

const bytes = (text) => new Uint8Array(Buffer.from(text, "base64"));
const want = F.skyPack;

function close(actual, expected, label, tolerance = 1e-6) {
  if (Number.isNaN(expected)) return assert.ok(Number.isNaN(actual), `${label}: ${actual} is not NaN`);
  const scale = Math.max(Math.abs(expected), 1e-12);
  assert.ok(Math.abs(actual - expected) <= tolerance * scale, `${label}: ${actual} vs ${expected}`);
}

test("a packed sky reads to the columns the server's own reader finds", () => {
  const sky = P.unpackSky(bytes(want.body));
  const expected = want.sky;
  assert.equal(sky.count, expected.count);
  assert.equal(sky.epochUnix, expected.epoch_unix);
  close(sky.limit, expected.limit, "limit");
  close(sky.velocityScale, expected.velocity_scale, "velocity scale");
  assert.deepEqual(sky.observerVelocityKms, expected.observer_velocity_kms);
  assert.deepEqual(sky.meta, expected.meta);
  for (let i = 0; i < sky.count; i++) {
    for (let axis = 0; axis < 3; axis++) {
      close(sky.offsetPc[3 * i + axis], expected.offset_pc[i][axis], `offset ${i} ${axis}`, 1e-6);
      close(sky.velocityKms[3 * i + axis], expected.velocity_kms[i][axis], `velocity ${i} ${axis}`, 1e-6);
    }
    close(sky.magnitude[i], expected.magnitude[i], `magnitude ${i}`);
    close(sky.temperatureK[i], expected.temperature_k[i], `temperature ${i}`, 1e-6);
    assert.equal(sky.kind[i], expected.kind[i], `kind ${i}`);
    assert.equal(sky.index[i], expected.index[i], `index ${i}`);
  }
  assert.equal(sky.offsetPc.length, 3 * sky.count);
});

test("the list is 26 bytes a star and the kinds and indices come through, including an index past 31 bits", () => {
  const sky = P.unpackSky(bytes(want.body));
  assert.equal(P.BYTES_PER_STAR, want.perStar);
  assert.equal(want.size, 52 + JSON.stringify(want.sky.meta).length + want.sky.count * want.perStar);
  assert.deepEqual(Array.from(sky.index).sort((a, b) => a - b), [0, 2, 7, 15, 2 ** 32 - 1].slice(0, sky.count).sort((a, b) => a - b));
  assert.deepEqual(Array.from(sky.kind), want.sky.kind);
  assert.deepEqual(Array.from(sky.kind).map((k) => P.KIND_NAMES[k]).sort(), want.sky.kind.map((k) => P.KIND_NAMES[k]).sort());
  assert.ok(sky.kind.length === 5 && new Set(sky.kind).size === 4, "all four kinds ride in the fixture");
});

test("an empty sky reads to no stars and its header", () => {
  const sky = P.unpackSky(bytes(want.emptyBody));
  assert.equal(sky.count, 0);
  assert.equal(sky.epochUnix, want.empty.epoch_unix);
  assert.equal(sky.offsetPc.length, 0);
});

test("bytes that are not one whole sky are refused", () => {
  const whole = bytes(want.body);
  assert.throws(() => P.unpackSky(whole.subarray(0, whole.length - 1)), /not a packed sky/);
  assert.throws(() => P.unpackSky(Uint8Array.from([...whole, 0])), /not a packed sky/);
  const wrong = whole.slice();
  wrong[0] = 0x58;
  assert.throws(() => P.unpackSky(wrong), /not a packed sky/);
  assert.throws(() => P.unpackSky(new Uint8Array(3)), /not a packed sky/);
});
