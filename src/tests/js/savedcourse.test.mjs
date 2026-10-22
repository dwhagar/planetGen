// tests/js/savedcourse.test.mjs -- static/savedcourse.js, the saved course record (NAV.17): both forms kept, the
// record's place in the saved-list store, and what opening it again says has changed.

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/nav");
const C = await import(new URL("savedcourse.js", STATIC_URL).href);

function answer(overrides) {
  return Object.assign({
    origin: { ref: "system:00001000000000-4000000-000", kind: "system", name: "Sol Fixture" },
    destination: { ref: "nebula:00001000000000-4000003-000", kind: "nebula", name: "Fixture Veil" },
    direct: { distance_ly: 12.345, bearing_deg: 45.2, mark_deg: 12.4, extra: "ignored" },
    route: {
      path: [1, 2, "phenomenon:nebula:3"],
      distance_ly: 13.5,
      hops: [{ distance_ly: 6, bearing_deg: 40, mark_deg: 10 }, { distance_ly: 7.5, bearing_deg: 50, mark_deg: 14 }],
      positions: { "1": [0, 0, 0] },
    },
  }, overrides || {});
}

const WHEN = new Date("2026-10-22T12:00:00Z");

test("a record keeps both forms, the two ends, a name and when", () => {
  const record = C.recordOf(answer(), "  Home run ", WHEN);
  assert.equal(record.kind, "course");
  assert.equal(record.name, "Home run");
  assert.equal(record.from, "system:00001000000000-4000000-000");
  assert.equal(record.to, "nebula:00001000000000-4000003-000");
  assert.equal(record.fromName, "Sol Fixture");
  assert.equal(record.saved, "2026-10-22T12:00:00.000Z");
  assert.deepEqual(record.direct, { distance_ly: 12.345, bearing_deg: 45.2, mark_deg: 12.4 });
  assert.deepEqual(record.route.path, ["1", "2", "phenomenon:nebula:3"]);
  assert.equal(record.route.distance_ly, 13.5);
  assert.deepEqual(record.route.hops[1], { distance_ly: 7.5, bearing_deg: 50, mark_deg: 14 });
  assert.equal(record.adjusted, null);
  assert.ok(!("positions" in record.route));
});

test("a course with no route keeps the direct line alone, and an unnamed one is named by its ends", () => {
  const record = C.recordOf(answer({ route: null }), "", WHEN);
  assert.equal(record.route, null);
  assert.equal(record.name, "Sol Fixture to Fixture Veil");
  assert.equal(record.direct.distance_ly, 12.345);
});

test("a record survives the store and a course is known by its two ends, in order", async () => {
  window.localStorage.clear();
  const store = C.courseStore(() => "galaxy_a");
  await store.ready;
  const record = C.recordOf(answer(), "Home run", WHEN);
  assert.equal((await store.add(record)).ok, true);
  assert.equal((await store.add(Object.assign({}, record, { name: "Same ends" }))).entry.name, "Home run", "same ends: already saved");
  const back = C.recordOf(answer({ origin: answer().destination, destination: answer().origin }), "Way back", WHEN);
  assert.equal((await store.add(back)).ok, true, "the other way round is another course");
  assert.equal(store.list().length, 2);
  assert.deepEqual(JSON.parse(window.localStorage.getItem("planetgen.courses.galaxy_a"))[0].route.path, ["1", "2", "phenomenon:nebula:3"]);
  assert.equal(await store.rename(C.courseKey(record), "Renamed"), true);
  assert.equal((await store.add({ name: "no ends" })).ok, false);
});

test("opening a course that has not changed says nothing", () => {
  const record = C.recordOf(answer(), "Home run", WHEN);
  assert.deepEqual(C.changesSince(record, answer()), []);
  const jitter = answer({ direct: { distance_ly: 12.3454, bearing_deg: 45.2004, mark_deg: 12.4 } });
  assert.deepEqual(C.changesSince(record, jitter), [], "a change below the tolerance is no change");
});

test("opening a course says what moved", () => {
  const record = C.recordOf(answer(), "Home run", WHEN);
  const moved = answer({
    direct: { distance_ly: 15, bearing_deg: 50, mark_deg: 12.4 },
    route: { path: [1, 5, "phenomenon:nebula:3"], distance_ly: 17.25, hops: [] },
  });
  const changes = C.changesSince(record, moved);
  assert.deepEqual(changes, [
    "The direct course is now 15 ly (it was 12.35 ly).",
    "The direct course heads 050 mark 012 (it was 045 mark 012).",
    "The route now passes through different systems (still 3 stops).",
    "The route is now 17.25 ly (it was 13.5 ly).",
  ]);
});

test("a route that gains a stop says how many", () => {
  const record = C.recordOf(answer(), "Home run", WHEN);
  const longer = answer({ route: { path: [1, 2, 4, "phenomenon:nebula:3"], distance_ly: 13.5, hops: [] } });
  assert.deepEqual(C.changesSince(record, longer), ["The route now passes through different systems (4 stops; it was 3)."]);
});

test("opening a course says when the route appears or goes", () => {
  const record = C.recordOf(answer(), "Home run", WHEN);
  assert.deepEqual(C.changesSince(record, answer({ route: null })), ["The route through known space is gone: no route is found now."]);
  const direct = C.recordOf(answer({ route: null }), "Direct only", WHEN);
  assert.deepEqual(C.changesSince(direct, answer()), ["A route through known space exists now (3 stops, 13.5 ly)."]);
});

test("a course whose end is gone says it can no longer be planned", () => {
  assert.equal(C.unavailableChange("no such system"), "The course can no longer be planned: no such system.");
  assert.equal(C.unavailableChange(""), "The course can no longer be planned: an end is gone.");
});

test("an adjusted path (NAV.6) that changes is reported once courses have one", () => {
  const record = C.recordOf(answer({ adjusted: { path: [[0, 0, 0], [1, 1, 1]] } }), "Bent", WHEN);
  assert.deepEqual(C.changesSince(record, answer({ adjusted: { path: [[0, 0, 0], [1, 1, 1]] } })), []);
  assert.deepEqual(C.changesSince(record, answer({ adjusted: { path: [[0, 0, 0], [2, 2, 2]] } })),
    ["The adjusted path around gravity wells has changed."]);
});
