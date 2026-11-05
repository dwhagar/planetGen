// tests/js/savedcourses.test.mjs -- static/savedcourses.js, the words and links of the saved courses list (NAV.4, NAV.18).

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/nav");
const S = await import(new URL("savedcourses.js", STATIC_URL).href);

const record = {
  name: "Home run", from: "system:00001000000000-4000000-000", to: "nebula:00001000000000-4000003-000",
  direct: { distance_ly: 12.345, bearing_deg: 45.2, mark_deg: 12.4 },
  route: { path: [1, 2, 3, 4], distance_ly: 13.5, hops: [] },
};

test("the two forms are direct and route", () => {
  assert.deepEqual(S.FORMS, ["direct", "route"]);
});

test("the direct form says the distance, the bearing and the mark", () => {
  const line = S.describeCourse(record, "direct");
  assert.match(line, /045 mark 012$/);
  assert.match(line, /12\.3 ly/);
});

test("a bearing below zero or at 360 wraps", () => {
  const line = S.describeCourse({ direct: { distance_ly: 1, bearing_deg: -10, mark_deg: 359.6 } }, "direct");
  assert.match(line, /350 mark 000$/);
});

test("the route form counts the stops between the ends", () => {
  assert.match(S.describeCourse(record, "route"), /, 2 stops between the ends$/);
  const one = Object.assign({}, record, { route: { path: [1, 2, 3], distance_ly: 5 } });
  assert.match(S.describeCourse(one, "route"), /, 1 stop between the ends$/);
});

test("a course without a form says so", () => {
  assert.equal(S.describeCourse(Object.assign({}, record, { route: null }), "route"), "No route through known space.");
  assert.equal(S.describeCourse(Object.assign({}, record, { direct: null }), "direct"), "No direct course.");
});

test("the links go to the NAV page and the Galaxy Map with the colons kept", () => {
  assert.equal(S.courseUrl(record), "/nav?from=system:00001000000000-4000000-000&to=nebula:00001000000000-4000003-000");
  assert.equal(S.mapUrl(record), "/galaxy?course=system:00001000000000-4000000-000,nebula:00001000000000-4000003-000");
});
