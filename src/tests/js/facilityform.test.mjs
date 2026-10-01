// tests/js/facilityform.test.mjs -- static/facilityform.js, the system
// page's "Place a facility" form (TEST.58): placement filters the hosts,
// the orbit slider shows only "in orbit" and its readout follows the
// slider on the same log scale and Kepler period as
// stellarObjects/facilities.py.

import { test } from "node:test";
import assert from "node:assert/strict";

import { FakeEvent, h, importPage, installDom, STATIC_URL } from "./fakedom.mjs";

const STEPS = 100;
const EARTH = { min: 6471, max: 924000, mass: 5.972e24 };

function form(placement, step) {
  const el = h("form", { class: "search-form facility-form" }, [
    h("select", { name: "placement", "data-facility-placement": "" }, [
      h("option", { value: "surface", selected: placement === "surface" }, "On the surface"),
      h("option", { value: "orbital", selected: placement === "orbital" }, "In orbit"),
      h("option", { value: "belt", selected: placement === "belt" }, "In a belt"),
    ]),
    h("select", { name: "host", "data-facility-host": "" }, [
      h("option", { value: "planet:1", "data-placements": "surface orbital", "data-min-km": EARTH.min, "data-max-km": EARTH.max, "data-mass-kg": EARTH.mass }, "Earth"),
      h("option", { value: "star:1", "data-placements": "orbital", "data-min-km": "700000", "data-max-km": "1e10", "data-mass-kg": "1.989e30" }, "Sun"),
      h("option", { value: "belt:1", "data-placements": "belt" }, "Main belt"),
    ]),
    h("div", { class: "facility-orbit", "data-facility-orbit": "", hidden: placement !== "orbital" }, [
      h("input", { type: "range", min: "0", max: String(STEPS), step: "1", value: String(step == null ? 50 : step) }),
      h("output", { "data-facility-orbit-readout": "" }),
    ]),
  ]);
  document.body.appendChild(el);
  return el;
}

const fmt = {};

async function load(placement, step) {
  installDom("http://localhost/system/4");
  if (!fmt.distance) {
    fmt.distance = (await import(new URL("distance.js", STATIC_URL).href)).formatDistanceKm;
    fmt.period = (await import(new URL("period.js", STATIC_URL).href)).formatDurationSeconds;
    fmt.speed = (await import(new URL("speed.js", STATIC_URL).href)).formatSpeedKms;
  }
  const el = form(placement, step);
  await importPage("facilityform.js");
  return {
    el,
    placement: el.querySelector("[data-facility-placement]"),
    host: el.querySelector("[data-facility-host]"),
    orbit: el.querySelector("[data-facility-orbit]"),
    slider: el.querySelector('input[type="range"]'),
    readout: el.querySelector("[data-facility-orbit-readout]"),
  };
}

function change(el, type) {
  el.dispatchEvent(new FakeEvent(type || "change"));
}

function expectedReadout(km, massKg) {
  const periodS = 2 * Math.PI * Math.sqrt((km * 1000) ** 3 / (6.6743e-11 * massKg));
  return `${fmt.distance(km)} from its host, one orbit every ${fmt.period(periodS)}, at ${fmt.speed((2 * Math.PI * km) / periodS)}`;
}

function shown(host) {
  return host.options.filter((o) => !o.hidden).map((o) => o.value);
}

test("each placement shows only the hosts that allow it", async () => {
  const f = await load("surface");
  assert.deepEqual(shown(f.host), ["planet:1"]);
  assert.equal(f.host.options[1].disabled, true, "a hidden host can't be sent either");
  f.placement.value = "orbital";
  change(f.placement);
  assert.deepEqual(shown(f.host), ["planet:1", "star:1"]);
  f.placement.value = "belt";
  change(f.placement);
  assert.deepEqual(shown(f.host), ["belt:1"]);
  assert.equal(f.host.value, "belt:1", "a host that no longer fits is replaced by the first that does");
});

test("the orbit slider shows only for an orbit", async () => {
  const f = await load("surface");
  assert.equal(f.orbit.hidden, true);
  f.placement.value = "orbital";
  change(f.placement);
  assert.equal(f.orbit.hidden, false);
  f.placement.value = "surface";
  change(f.placement);
  assert.equal(f.orbit.hidden, true);
});

test("the readout runs from just above the surface to the edge of the sphere of influence", async () => {
  const f = await load("orbital", 0);
  assert.equal(f.readout.textContent, expectedReadout(EARTH.min, EARTH.mass));
  assert.equal(f.slider.getAttribute("aria-valuetext"), fmt.distance(EARTH.min));
  f.slider.value = String(STEPS);
  change(f.slider, "input");
  assert.equal(f.readout.textContent, expectedReadout(EARTH.max, EARTH.mass));
  f.slider.value = String(STEPS / 2);
  change(f.slider, "input");
  // Log scale: the middle step is the geometric mean.
  assert.equal(f.readout.textContent, expectedReadout(Math.sqrt(EARTH.min * EARTH.max), EARTH.mass));
});

test("a geostationary distance reads one orbit a sidereal day", async () => {
  const f = await load("orbital", 0);
  // facilities.distance_from_step's inverse for 42,164 km.
  const step = (STEPS * Math.log(42164 / EARTH.min)) / Math.log(EARTH.max / EARTH.min);
  f.slider.value = String(step);
  change(f.slider, "input");
  assert.match(f.readout.textContent, /one orbit every 23\.9 hours, at 3\.07 km\/s/);
});

test("out-of-range slider values are clamped like distance_from_step", async () => {
  const f = await load("orbital", 0);
  f.slider.value = "-5";
  change(f.slider, "input");
  assert.equal(f.readout.textContent, expectedReadout(EARTH.min, EARTH.mass));
  f.slider.value = String(STEPS + 50);
  change(f.slider, "input");
  assert.equal(f.readout.textContent, expectedReadout(EARTH.max, EARTH.mass));
});

test("changing the host redraws the readout for that host", async () => {
  const f = await load("orbital", 0);
  f.host.value = "star:1";
  change(f.host);
  assert.equal(f.readout.textContent, expectedReadout(700000, 1.989e30));
});

test("a host without orbit limits leaves the readout empty", async () => {
  const f = await load("belt");
  assert.equal(f.readout.textContent, "");
});

test("a form missing a part is left alone", async () => {
  installDom("http://localhost/system/4");
  const el = form("orbital", 0);
  el.querySelector("[data-facility-orbit-readout]").remove();
  await importPage("facilityform.js");
  assert.equal(el.querySelector("[data-facility-host]").options[2].hidden, false);
});
