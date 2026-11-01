// tests/js/daycounter.test.mjs -- static/daycounter.js, the day-and-night counter of an object page (UX.243).

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const C = await import(new URL("daycounter.js", STATIC_URL).href);

const EPOCH = 946728000000;
const HOUR = C.MS_PER_HOUR;

test("the first moment of the count is the start of day 0", () => {
  const r = C.clockAt(EPOCH, EPOCH, 24, 12);
  assert.equal(r.day, 0);
  assert.equal(r.hoursIn, 0);
  assert.equal(r.isDay, true);
  assert.equal(r.hoursLeft, 12);
});

test("an Earth day: half a day in is the start of the night, 30 hours in is the next day", () => {
  assert.equal(C.clockAt(EPOCH + 12 * HOUR, EPOCH, 24, 12).isDay, false);
  const next = C.clockAt(EPOCH + 30 * HOUR, EPOCH, 24, 12);
  assert.equal(next.day, 1);
  assert.equal(next.hoursIn, 6);
  assert.equal(next.isDay, true);
  assert.equal(next.hoursLeft, 6);
});

test("the night counts down to the next sunrise", () => {
  const r = C.clockAt(EPOCH + 18 * HOUR, EPOCH, 24, 12);
  assert.equal(r.isDay, false);
  assert.equal(r.hoursLeft, 6);
});

test("a day of many hours counts at the real rate, not a compressed one", () => {
  const r = C.clockAt(EPOCH + 100 * HOUR, EPOCH, 5832, 2916);
  assert.equal(r.day, 0);
  assert.equal(r.hoursIn, 100);
});

test("a night not known leaves day or night unsaid, a night as long as the day is all night", () => {
  assert.equal(C.clockAt(EPOCH + HOUR, EPOCH, 24, null).isDay, null);
  assert.equal(C.clockAt(EPOCH + HOUR, EPOCH, 24, 24).isDay, false);
});

test("no usable day gives no clock", () => {
  assert.equal(C.clockAt(EPOCH, EPOCH, 0, 0), null);
  assert.equal(C.clockAt(EPOCH, EPOCH, NaN, 1), null);
  assert.equal(C.sentences(null, 24, 12), null);
});

test("a span is written with its units", () => {
  assert.equal(C.spanText(7.2092), "7 h 12 min 33 s");
  assert.equal(C.spanText(0), "0 h 00 min 00 s");
  assert.equal(C.spanText(-1), "–");
});

test("the sentences say the day number, the time in, and the phase with the time left", () => {
  const reading = C.clockAt(EPOCH + 30 * HOUR, EPOCH, 24, 12);
  const parts = C.sentences(reading, 24, 12);
  assert.match(parts.day, /^Day 2 since the count began$/);
  assert.match(parts.inDay, /^6 h 00 min 00 s into a day of 24 h 00 min 00 s$/);
  assert.match(parts.phase, /^It is day: 6 h 00 min 00 s of daylight left\.$/);
});

test("a body with no night says so", () => {
  const reading = C.clockAt(EPOCH + HOUR, EPOCH, 24, 0);
  assert.match(C.sentences(reading, 24, 0).phase, /^No night/);
});

test("start fills the three lines and the stop function ends the timer", () => {
  const nodes = {};
  const element = {
    dataset: { dayHours: "24", nightHours: "12", epochMs: String(EPOCH) },
    querySelector: (selector) => (nodes[selector] ??= { textContent: "" }),
  };
  let cleared = null;
  const stop = C.start(element, () => EPOCH + 30 * HOUR, () => 7, (timer) => { cleared = timer; });
  assert.match(nodes["[data-day-counter-day]"].textContent, /Day 2/);
  assert.match(nodes["[data-day-counter-in-day]"].textContent, /6 h 00 min 00 s into a day/);
  assert.match(nodes["[data-day-counter-phase]"].textContent, /It is day/);
  stop();
  assert.equal(cleared, 7);
});
