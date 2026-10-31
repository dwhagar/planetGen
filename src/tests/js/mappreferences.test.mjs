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

test("the label buffer reads as a number of pixels from 0 to 40 and is 4 when missing or out of range", () => {
  assert.equal(P.labelBufferPx(store({})), 4);
  assert.equal(P.labelBufferPx(null), 4);
  assert.equal(P.labelBufferPx(store({ "planetgen.pref.map.labelBuffer": "12" })), 12);
  assert.equal(P.labelBufferPx(store({ "planetgen.pref.map.labelBuffer": "0" })), 0);
  assert.equal(P.labelBufferPx(store({ "planetgen.pref.map.labelBuffer": "41" })), 4);
  assert.equal(P.labelBufferPx(store({ "planetgen.pref.map.labelBuffer": "wide" })), 4);
});

test("the region landing reads frame unless the stored word is densest", () => {
  assert.equal(P.regionLanding(store({})), "frame");
  assert.equal(P.regionLanding(null), "frame");
  assert.equal(P.regionLanding(store({ "planetgen.pref.map.regionLanding": "densest" })), "densest");
  assert.equal(P.regionLanding(store({ "planetgen.pref.map.regionLanding": " frame " })), "frame");
  ["", "Densest", "closest", "0", "[densest]"].forEach((text) => {
    assert.equal(P.regionLanding(store({ "planetgen.pref.map.regionLanding": text })), "frame", text);
  });
  assert.equal(P.regionLanding({ getItem() { throw new Error("blocked"); } }), "frame");
  assert.throws(() => P.readChoice(store({}), "map.nope"));
});

// --- UX.236: writing the settings, and what the preferences page lists ---

const memory = (values = {}) => ({
  getItem: (key) => (key in values ? values[key] : null),
  setItem: (key, value) => { values[key] = String(value); },
  removeItem: (key) => { delete values[key]; },
  values,
});

test("every setting has its words on the preferences page, in the order they are listed", () => {
  const keys = [...Object.keys(P.NUMBER_SETTINGS), ...Object.keys(P.CHOICE_SETTINGS)].sort();
  assert.deepEqual(Object.keys(P.PREFERENCE_TEXT).sort(), keys);
  const list = P.preferenceList();
  assert.deepEqual(list.map((entry) => entry.key), Object.keys(P.PREFERENCE_TEXT));
  list.forEach((entry) => {
    assert.ok(entry.label && entry.description && entry.group, entry.key);
    if (entry.kind === "choice") {
      assert.ok(Object.keys(entry.choices).includes(entry.fallback), entry.key);
      // the words the page offers are the words the reader accepts
      assert.deepEqual(Object.keys(entry.choices).sort(), P.CHOICE_SETTINGS[entry.key].choices.slice().sort(), entry.key);
    } else {
      assert.equal(entry.kind, "number");
      assert.ok(entry.unit && entry.min < entry.max && entry.fallback >= entry.min && entry.fallback <= entry.max, entry.key);
    }
  });
});

test("the first two settings are where a region lands (frame by default) and the black hole mark density (one per 120 pixels)", () => {
  const [landing, density] = P.preferenceList();
  assert.equal(landing.key, "map.regionLanding");
  assert.equal(landing.fallback, "frame");
  assert.deepEqual(Object.keys(landing.choices), ["frame", "densest"]);
  assert.equal(density.key, "map.blackHoleDensity");
  assert.equal(density.fallback, 120);
  assert.equal(density.min, 30);
  assert.equal(density.max, 1000);
});

test("a written word or number is what the map then reads", () => {
  const store = memory();
  assert.equal(P.writeSetting(store, "map.regionLanding", "densest"), true);
  assert.equal(P.regionLanding(store), "densest");
  assert.equal(store.values["planetgen.pref.map.regionLanding"], "densest");
  assert.equal(P.writeSetting(store, "map.blackHoleDensity", " 60 "), true);
  assert.equal(P.blackHoleDensityPx(store), 60);
  assert.equal(P.readSetting(store, "map.blackHoleDensity"), 60);
  assert.equal(P.readSetting(store, "map.regionLanding"), "densest");
});

test("a word or number the reader would not accept is not written", () => {
  const store = memory();
  ["closest", "", "Densest"].forEach((text) => assert.equal(P.writeSetting(store, "map.regionLanding", text), false, text));
  ["", "abc", "29", "1001", "-5", "NaN"].forEach((text) => assert.equal(P.writeSetting(store, "map.blackHoleDensity", text), false, text));
  assert.deepEqual(store.values, {});
  assert.throws(() => P.writeSetting(store, "map.nope", "1"));
});

test("a store that refuses a write says so", () => {
  const refusing = { getItem: () => null, setItem() { throw new Error("full"); }, removeItem() { throw new Error("blocked"); } };
  assert.equal(P.writeSetting(refusing, "map.regionLanding", "densest"), false);
  assert.equal(P.clearSetting(refusing, "map.regionLanding"), false);
  assert.equal(P.resetSettings(refusing), false);
});

test("clearing one setting, or resetting all of them, brings back the defaults and leaves other keys alone", () => {
  const store = memory({ "planetgen.pref.map.regionLanding": "densest", "planetgen.pref.map.blackHoleDensity": "60",
                         "planetgen.pref.map.labelBuffer": "9", "planetgen.units": "metric" });
  assert.equal(P.clearSetting(store, "map.blackHoleDensity"), true);
  assert.equal(P.blackHoleDensityPx(store), 120);
  assert.equal(P.regionLanding(store), "densest");
  assert.equal(P.resetSettings(store), true);
  assert.equal(P.regionLanding(store), "frame");
  assert.equal(P.labelBufferPx(store), 4);
  assert.deepEqual(store.values, { "planetgen.units": "metric" });
});
