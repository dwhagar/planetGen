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
