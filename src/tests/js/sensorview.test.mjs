// tests/js/sensorview.test.mjs -- static/sensorview.js (UX.275): the painter table of the Sensor view.

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, importPage } from "./fakedom.mjs";

installDom();
const S = await importPage("sensorview.js");

test("the painters are the three that exist, each with the module that draws it", () => {
  assert.deepEqual(Object.keys(S.PAINTERS), ["sphere", "phenomenon", "nebula"]);
  assert.equal(S.PAINTERS.sphere.module, "bodyview.js");
  assert.equal(S.PAINTERS.phenomenon.module, "phenomenonrender.js");
  assert.equal(S.PAINTERS.nebula.module, "nebulaview.js");
});

test("a record is drawn by the painter it names, and by none when it names none", () => {
  assert.equal(S.painterOf({ v: 1, painter: "sphere", look: {} }).module, "bodyview.js");
  assert.equal(S.painterOf({ v: 1, painter: "phenomenon" }).id, "phenomenon");
  assert.equal(S.painterOf({ v: 1, painter: null, reason: "none yet" }), null);
  assert.equal(S.painterOf({ v: 1, painter: "hologram" }), null, "an unknown painter is none");
  assert.equal(S.painterOf({ v: 1, painter: "toString" }), null, "an inherited name is not a painter");
  assert.equal(S.painterOf({ v: S.SENSOR_VERSION + 1, painter: "sphere" }), null, "a newer record than this script reads");
  assert.equal(S.painterOf(null), null);
  assert.equal(S.painterOf("sphere"), null);
});
