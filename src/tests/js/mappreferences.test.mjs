// tests/js/mappreferences.test.mjs -- static/mappreferences.js (MAP.339): the map's per-browser number settings read with a default.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const P = await import(new URL("mappreferences.js", STATIC_URL).href);

const store = (values) => ({ getItem: (key) => (key in values ? values[key] : null) });

test("a missing key reads as the default: one black hole mark per 120 by 120 pixels", () => {
  assert.equal(P.blackHoleDensityPx(store({})), 120);
  assert.equal(P.blackHoleDensityPx(null), 120);
  assert.equal(P.blackHoleDensityPx(undefined), 120);
});

test("a stored number inside the range is read", () => {
  assert.equal(P.blackHoleDensityPx(store({ "planetgen.pref.map.blackHoleDensity": "60" })), 60);
  assert.equal(P.blackHoleDensityPx(store({ "planetgen.pref.map.blackHoleDensity": " 200.5 " })), 200.5);
});

test("a value that is not a number, or is out of range, reads as the default", () => {
  ["", "abc", "NaN", "Infinity", "-5", "0", "29", "1001", "[1]"].forEach((text) => {
    assert.equal(P.blackHoleDensityPx(store({ "planetgen.pref.map.blackHoleDensity": text })), 120, text);
  });
});

test("storage that throws reads as the default and an unknown key is an error", () => {
  assert.equal(P.blackHoleDensityPx({ getItem() { throw new Error("blocked"); } }), 120);
  assert.throws(() => P.readNumber(store({}), "map.nope"));
});
