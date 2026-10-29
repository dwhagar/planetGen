// tests/js/poiprefs.test.mjs -- static/poiprefs.js (MAP.306): the point-of-interest choices a browser remembers.

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/galaxy");
const L = await import(new URL("maplayers.js", STATIC_URL).href);
const P = await import(new URL("poiprefs.js", STATIC_URL).href);

const memory = () => { const data = {}; return { getItem: (k) => (k in data ? data[k] : null), setItem: (k, v) => { data[k] = v; } }; };
const IDS = ["a", "b"];

function registry() {
  const reg = L.createLayerRegistry("galaxy");
  reg.register({ id: "a", label: "A", forward: true, urlKey: "hide", urlKeyOn: "show" });
  reg.register({ id: "b", label: "B", forward: true, defaultOn: false, urlKey: "hide", urlKeyOn: "show" });
  return reg;
}
const DEFAULTS = { a: "shown", b: "hidden" };

test("the states are saved and come back on a fresh map", () => {
  const storage = memory();
  const first = registry();
  first.setPulled("b", true);
  first.setOn("a", false);
  P.savePrefs(storage, first, IDS);
  const second = registry();
  assert.equal(P.applyPrefs(storage, second, IDS, DEFAULTS), true);
  assert.deepEqual(IDS.map((id) => second.stateOf(id)), ["hidden", "forward"]);
});

test("a map whose address chose a layer is left as the address says", () => {
  const storage = memory();
  const first = registry();
  first.setOn("a", false);
  P.savePrefs(storage, first, IDS);
  const second = registry();
  second.setOn("b", true);
  assert.equal(P.applyPrefs(storage, second, IDS, DEFAULTS), false);
  assert.deepEqual(IDS.map((id) => second.stateOf(id)), ["shown", "shown"]);
});

test("a store that is missing, empty, damaged or refuses writes is no trouble", () => {
  const reg = registry();
  assert.equal(P.applyPrefs(memory(), reg, IDS, DEFAULTS), false);
  const damaged = { getItem: () => "{nope", setItem: () => { throw new Error("full"); } };
  assert.deepEqual(P.loadPrefs(damaged), {});
  assert.doesNotThrow(() => P.savePrefs(damaged, reg, IDS));
  assert.deepEqual(P.loadPrefs({ getItem: () => JSON.stringify({ a: "forward", b: "sideways" }) }), { a: "forward" });
});
