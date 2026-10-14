// tests/js/nebulalook.test.mjs -- static/nebulalook.js, the Galaxy Map's
// nebula colours per theme (MAP.175): every family reaches 3:1 against the
// map's real background on both themes, and the contrast arithmetic is
// WCAG's.

import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";

import { importPage, STATIC_URL } from "./fakedom.mjs";

const N = await importPage("nebulalook.js");

// The map's background as the page paints it: --bg-subtle over the page's
// --bg, from style.css itself, so a change to the theme's colours fails here.
const css = readFileSync(new URL("style.css", STATIC_URL), "utf8");

function token(block, name) {
  const match = new RegExp("--" + name + ":\\s*#([0-9a-fA-F]{6})([0-9a-fA-F]{2})?\\s*;").exec(block);
  assert.ok(match, `${name} in the block`);
  return { rgb: N.parseHex("#" + match[1]), alpha: match[2] ? parseInt(match[2], 16) / 255 : 1 };
}

function background(block) {
  const bg = token(block, "bg");
  const subtle = token(block, "bg-subtle");
  return N.compositeOver(subtle.rgb, subtle.alpha, bg.rgb);
}

const rootBlock = css.slice(css.indexOf(":root {"), css.indexOf("}", css.indexOf(":root {")));
const darkStart = css.indexOf(':root[data-theme="dark"] {');
const darkBlock = css.slice(darkStart, css.indexOf("}", darkStart));
const LIGHT = background(rootBlock);
const DARK = background(darkBlock);

test("the contrast arithmetic is WCAG's", () => {
  assert.ok(Math.abs(N.contrastRatio([0, 0, 0], [255, 255, 255]) - 21) < 1e-9);
  assert.equal(N.contrastRatio([10, 20, 30], [10, 20, 30]), 1);
  assert.deepEqual(N.compositeOver([200, 100, 0], 0.5, [0, 0, 100]), [100, 50, 50]);
  assert.deepEqual(N.parseHex("#1c1c24"), [28, 28, 36]);
});

test("the two themes' backgrounds are the light and dark ones the map sits on", () => {
  assert.ok(N.luminance(LIGHT) > 0.7, String(LIGHT));
  assert.ok(N.luminance(DARK) < 0.03, String(DARK));
});

test("every nebula family reaches 3:1 against the map's background on both themes", () => {
  [["light", N.NEBULA_LOOKS_LIGHT_THEME, LIGHT], ["dark", N.NEBULA_LOOKS_DARK_THEME, DARK]].forEach(([name, looks, bg]) => {
    ["diffuse", "emission", "reflection", "planetary", "dark", "default"].forEach((family) => {
      assert.ok(looks[family], `${family} has a look on the ${name} theme`);
      const ratio = N.lookContrast(looks[family], bg);
      assert.ok(ratio >= N.CONTRAST_TARGET, `${family} on the ${name} theme is ${ratio.toFixed(2)}:1`);
    });
  });
});

test("the dark family was invisible before and is not now", () => {
  // The old near-black fill on the dark theme.
  assert.ok(N.lookContrast(["#1c1c24", 0.91, 0.56], DARK) < 1.3);
  assert.ok(N.lookContrast(N.NEBULA_LOOKS_DARK_THEME.dark, DARK) >= 3);
});

test("each theme gets its own set and both have the same families", () => {
  assert.equal(N.nebulaLooks(true), N.NEBULA_LOOKS_LIGHT_THEME);
  assert.equal(N.nebulaLooks(false), N.NEBULA_LOOKS_DARK_THEME);
  assert.deepEqual(Object.keys(N.NEBULA_LOOKS_LIGHT_THEME).sort(), Object.keys(N.NEBULA_LOOKS_DARK_THEME).sort());
  Object.values(N.NEBULA_LOOKS_LIGHT_THEME).concat(Object.values(N.NEBULA_LOOKS_DARK_THEME)).forEach((look) => {
    assert.match(look[0], /^#[0-9a-f]{6}$/);
    assert.ok(look[1] > 0 && look[1] <= 1 && look[2] > 0 && look[2] <= look[1]);
  });
});
