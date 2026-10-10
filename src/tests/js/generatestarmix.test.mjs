// tests/js/generatestarmix.test.mjs -- static/generatestarmix.js, the Generate
// page's star-mix total (ADM.45).

import { test } from "node:test";
import assert from "node:assert/strict";

import { FakeEvent, h, installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/");
const mod = await import(new URL("generatestarmix.js", STATIC_URL).href);

test("the total text says what to do", () => {
  assert.match(mod.totalText(100), /That is right/);
  assert.match(mod.totalText(100.004), /That is right/);
  assert.match(mod.totalText(110), /Lower one of them by 10 points/);
  assert.match(mod.totalText(95.5), /Raise one of them by 4.5 points/);
  assert.match(mod.totalText(NaN), /needs a number/);
});

test("the box keeps a running total as the shares change", () => {
  const inputs = ["69.7", "15.2", "15.1"].map((value) => h("input", { value }));
  const out = h("p", { "data-star-mix-total": "" });
  document.body.appendChild(h("fieldset", { "data-star-mix": "" }, [...inputs, out]));
  mod.init();
  assert.match(out.textContent, /That is right/);
  inputs[0].value = "80";
  inputs[0].dispatchEvent(new FakeEvent("input"));
  assert.match(out.textContent, /Lower one of them by 10\.3 points/);
});
