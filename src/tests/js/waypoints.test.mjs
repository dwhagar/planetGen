// tests/js/waypoints.test.mjs -- static/waypoints.js, the stops picked for a course (NAV.49): what can be a waypoint, the list in the
// order picked, and the words on the Waypoint button.

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/galaxy");
const W = await import(new URL("waypoints.js", STATIC_URL).href);

test("a star system or a phenomenon with a NAV endpoint can be a waypoint", () => {
  assert.equal(W.isWaypoint({ name: "Sol", value: "system:00001000000000-4000000-000", kind: "system" }), true);
  assert.equal(W.isWaypoint({ name: "Veil", value: "nebula:3" }), true);
});

test("a map view, a sector, or anything without a NAV endpoint cannot", () => {
  assert.equal(W.isWaypoint({ name: "View", kind: "stage", value: "/galaxy?cam=1,2,3" }), false);
  assert.equal(W.isWaypoint({ name: "Sector", kind: "sector", value: "sector:12" }), false);
  assert.equal(W.isWaypoint({ name: "Page", value: "/system/12" }), false);
  assert.equal(W.isWaypoint(null), false);
  assert.equal(W.isWaypoint({ value: "system:1" }), false);
});

test("the button says to add, or to remove with the place in the course", () => {
  const entry = { name: "Sol", value: "system:1" };
  assert.equal(W.buttonWords(entry, 0), "Add waypoint");
  assert.equal(W.buttonWords(entry, 3), "Remove waypoint 3");
  assert.equal(W.buttonWords(null, 0), "Waypoint");
});

test("with no database named nothing is kept, and adding says so", async () => {
  await W.ready;
  assert.deepEqual(W.list(), []);
  const result = await W.add({ name: "Sol", value: "system:1" });
  assert.equal(result.ok, false);
});
