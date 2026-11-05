// tests/js/coursehops.test.mjs -- static/coursehops.js, the hops of the Galaxy Map's course that cross unknown space (NAV.36).

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const H = await import(new URL("coursehops.js", STATIC_URL).href);

test("known hops that follow each other make one run and each unknown hop is a run of its own", () => {
  const point = (hopUnknown) => ({ hopUnknown });
  assert.deepEqual(H.courseRuns([point(false), point(false), point(false)]), [{ from: 0, to: 2, unknown: false }]);
  assert.deepEqual(
    H.courseRuns([point(false), point(false), point(true), point(true), point(false), point(false)]),
    [{ from: 0, to: 1, unknown: false }, { from: 1, to: 2, unknown: true }, { from: 2, to: 3, unknown: true }, { from: 3, to: 5, unknown: false }]);
  assert.deepEqual(H.courseRuns([point(false), point(true)]), [{ from: 0, to: 1, unknown: true }]);
  assert.deepEqual(H.courseRuns([point(false)]), []);
});

test("a hop takes one glow sprite per few screen pixels, between two and the most", () => {
  assert.equal(H.beadCount(0), 2);
  assert.equal(H.beadCount(60), 10);
  assert.equal(H.beadCount(1e9), 240);
  assert.equal(H.beadCount(Infinity), 240);
});

test("the glow pulses between 0.3 and 1, and holds still under reduced motion", () => {
  for (let now = 0; now < 5000; now += 137) {
    const f = H.pulseFactor(now, false);
    assert.ok(f >= 0.29 && f <= 1.0001, String(f));
    assert.equal(H.pulseFactor(now, true), 1);
  }
  assert.notEqual(H.pulseFactor(0, false), H.pulseFactor(600, false));
});
