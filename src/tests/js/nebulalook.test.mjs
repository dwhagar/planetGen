// tests/js/nebulalook.test.mjs -- static/nebulalook.js, the Galaxy Map's
// nebula colours (MAP.175): every family reaches 3:1 against the map's black
// background, and the contrast arithmetic is WCAG's.

import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";

import { importPage, STATIC_URL } from "./fakedom.mjs";

const N = await importPage("nebulalook.js");

const css = readFileSync(new URL("style.css", STATIC_URL), "utf8");

const rootBlock = css.slice(css.indexOf(":root {"), css.indexOf("}", css.indexOf(":root {")));
// Every map is drawn on black in either theme (`--map-bg`, Boss 2026-10-31).
const DARK = N.parseHex("#" + /--map-bg:\s*#([0-9a-fA-F]{6})\s*;/.exec(rootBlock)[1]);

test("the contrast arithmetic is WCAG's", () => {
  assert.ok(Math.abs(N.contrastRatio([0, 0, 0], [255, 255, 255]) - 21) < 1e-9);
  assert.equal(N.contrastRatio([10, 20, 30], [10, 20, 30]), 1);
  assert.deepEqual(N.compositeOver([200, 100, 0], 0.5, [0, 0, 100]), [100, 50, 50]);
  assert.deepEqual(N.parseHex("#1c1c24"), [28, 28, 36]);
});

test("the map itself is black", () => {
  assert.deepEqual(DARK, [0, 0, 0]);
});

test("every nebula family reaches 3:1 against the map's black background", () => {
  ["diffuse", "emission", "reflection", "planetary", "dark", "default"].forEach((family) => {
    assert.ok(N.NEBULA_LOOKS[family], `${family} has a look`);
    const ratio = N.lookContrast(N.NEBULA_LOOKS[family], DARK);
    assert.ok(ratio >= N.CONTRAST_TARGET, `${family} is ${ratio.toFixed(2)}:1`);
  });
});

test("the dark family was invisible before and is not now", () => {
  // The old near-black fill on the dark theme.
  assert.ok(N.lookContrast(["#1c1c24", 0.91, 0.56], DARK) < 1.3);
  assert.ok(N.lookContrast(N.NEBULA_LOOKS.dark, DARK) >= 3);
});

test("every look is a hex colour with a fill alpha above the rim's", () => {
  Object.values(N.NEBULA_LOOKS).forEach((look) => {
    assert.match(look[0], /^#[0-9a-f]{6}$/);
    assert.ok(look[1] > 0 && look[1] <= 1 && look[2] > 0 && look[2] <= look[1]);
  });
});
