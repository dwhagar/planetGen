// tests/js/starlight.test.mjs -- static/starlight.js, the brightness
// boost for faint stars on the Galaxy Map (MAP.87), against its Python
// twin in lib/starmap.py (star_light_boost, _boost_light), which works out
// the Sector Map's points: one curve for both maps.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

const L = await importPage("starlight.js");

test("the boost matches starmap.py's at every luminosity", () => {
  for (const row of F.starLight) {
    const boost = L.starLightBoost(row.luminositySol);
    assert.ok(Math.abs(boost - row.boost) < 1e-9, `${row.luminositySol} L_sun: ${boost} vs ${row.boost}`);
    const light = L.boostLight({ sizePx: 10, glow: 0.5 }, boost);
    assert.ok(Math.abs(light.sizePx - row.sizePx) < 1e-9 && Math.abs(light.glow - row.glow) < 1e-9, JSON.stringify(row));
  }
});

test("four times at the dim end, none from 1000 suns up, 1 without a luminosity", () => {
  assert.ok(Math.abs(L.starLightBoost(1e-4) - 4) < 1e-12);
  assert.equal(L.starLightBoost(1000), 1);
  assert.equal(L.starLightBoost(5e5), 1);
  assert.equal(L.starLightBoost(null), 1);
  assert.equal(L.starLightBoost(0), 1);
});

test("a star's look matches starmap.py's at every luminosity and radius (MAP.241)", () => {
  assert.ok(F.starLook.length > 10);
  for (const row of F.starLook) {
    const light = L.starLight(row.luminositySol, row.radiusSol);
    const where = `${row.luminositySol} L_sun, ${row.radiusSol} R_sun`;
    for (const key of ["corePx", "sizePx", "bright"]) assert.ok(Math.abs(light[key] - row[key]) < 0.011, `${where}: ${key} ${light[key]} vs ${row[key]}`);
    assert.ok(Math.abs(light.glow - row.glow) < 0.0011, `${where}: glow ${light.glow} vs ${row.glow}`);
  }
});

test("a brighter star glows harder but is not bigger: sizes stay in a narrow band, the glow carries the luminosity (MAP.241)", () => {
  const dim = L.starLight(1e-3, 0.3);
  const bright = L.starLight(1e5, 0.3);
  const sun = L.starLight(1, 1);
  for (const light of [dim, sun, bright]) assert.ok(light.sizePx >= 8 && light.sizePx <= 13, `${light.sizePx}`);
  assert.ok(Math.abs(bright.sizePx - dim.sizePx) < 3, `${bright.sizePx} vs ${dim.sizePx}`);
  assert.ok(bright.glow > 2 * dim.glow, `${bright.glow} vs ${dim.glow}`);
  assert.ok(bright.bright >= dim.bright);
  // A supergiant's core is a little larger than a dwarf's, still small.
  assert.ok(L.starLight(1, 500).corePx <= 2.4 + 1e-9 && L.starLight(1, 500).corePx > L.starLight(1, 0.1).corePx);
  // Nothing is made up for a star with no luminosity or radius on record.
  assert.ok(Number.isFinite(L.starLight(null, null).sizePx));
});
