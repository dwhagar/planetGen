// tests/js/superbubble.test.mjs -- static/superbubble.js (MAP.236): how a superbubble is drawn, a faint fill and a crisp rim.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const B = await importPage("superbubble.js");

test("the fill is fainter than the old circles' 0.14 and the rim stands out from it", () => {
  const look = B.LOOK;
  assert.ok(look.fillAlpha > 0 && look.fillAlpha < 0.14, `${look.fillAlpha}`);
  assert.ok(look.rimAlpha >= 3 * look.fillAlpha, `${look.rimAlpha}`);
  assert.ok(look.rimAlpha <= 1);
});

test("a colour is three components between 0 and 1 as written to the screen", () => {
  assert.deepEqual(B.colorComponents("#ff8000"), [1, 128 / 255, 0]);
  assert.deepEqual(B.colorComponents("#000000"), [0, 0, 0]);
});

test("the shader draws the rim where the wall meets the line of sight, from the fill's alpha to the rim's", () => {
  assert.match(B.FRAGMENT_SHADER, /fwidth\(facing\)/);
  assert.match(B.FRAGMENT_SHADER, /mix\(fillAlpha, rimAlpha, rim\)/);
  assert.match(B.VERTEX_SHADER, /modelViewMatrix/);
  assert.ok(B.RIM_PX > 0);
});
