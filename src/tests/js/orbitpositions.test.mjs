// tests/js/orbitpositions.test.mjs -- static/orbitpositions.js (MAP.70), the
// browser twin of physics/body_positions.py, checked against Python's own
// positions for the same scene, and static/orbitclock.js, the time control.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

const P = await import(new URL("orbitpositions.js", STATIC_URL).href);
const C = await import(new URL("orbitclock.js", STATIC_URL).href);

const AU_KM = 149597870.7;

test("positionsAt matches the Python positions for every body at every time", () => {
  const { scene, cases } = F.orbitPositions;
  for (const c of cases) {
    const got = P.positionsAt(scene, c.years);
    assert.deepEqual(Object.keys(got).sort(), Object.keys(c.positions).sort());
    for (const ref of Object.keys(c.positions)) {
      c.positions[ref].forEach((expected, axis) => {
        const tolerance = Math.max(1, Math.abs(expected) * 1e-9);
        assert.ok(Math.abs(got[ref][axis] - expected) <= tolerance,
          `${ref} axis ${axis} at ${c.years} y: ${got[ref][axis]} vs ${expected}`);
      });
    }
  }
});

test("a planet keeps its radius and returns after a period", () => {
  const { scene } = F.orbitPositions;
  const start = P.positionsAt(scene, 0)["planet:1"];
  const later = P.positionsAt(scene, 1.8)["planet:1"];
  assert.ok(Math.abs(Math.hypot(...start) - 1.5 * AU_KM) < 1e-3 * AU_KM * 1e-6 + 1);
  assert.ok(Math.hypot(...start.map((v, i) => v - later[i])) < 1e4);
});

test("the clock starts at the present, runs in real time and speeds up", () => {
  const epoch = 1_000_000;
  const clock = C.createOrbitClock({ epochUnix: epoch, nowMs: (epoch + 86400) * 1000 });
  const start = clock.years();
  assert.ok(Math.abs(start - 86400 / C.SECONDS_PER_YEAR) < 1e-15);
  clock.tick((epoch + 86400) * 1000 + 1000);
  assert.ok(Math.abs(clock.years() - start - 1 / C.SECONDS_PER_YEAR) < 1e-15);
  clock.faster();
  clock.faster();
  clock.faster();
  const before = clock.years();
  clock.tick((epoch + 86400) * 1000 + 2000);
  assert.ok(Math.abs(clock.years() - before - 1) < 1e-12, "1 year per second at the 4th step");
});

test("pause freezes the clock and now jumps back to the present", () => {
  const epoch = 0;
  const t0 = 5_000_000;
  const clock = C.createOrbitClock({ epochUnix: epoch, nowMs: t0 });
  clock.faster();
  clock.tick(t0 + 1000);
  clock.pause();
  const frozen = clock.years();
  clock.tick(t0 + 5000);
  assert.equal(clock.years(), frozen);
  clock.now(t0 + 6000);
  assert.ok(Math.abs(clock.years() - (t0 + 6000) / 1000 / C.SECONDS_PER_YEAR) < 1e-15);
  assert.equal(clock.rateIndex(), 0);
  assert.equal(clock.playing(), true);
});

test("a comet's orbit is described as bound with its period, or unbound and not returning (GEN.182)", () => {
  const kepler = { perihelion_distance_km: 2 * AU_KM, eccentricity: 0.5, period_years: 12.3456 };
  const closed = Object.fromEntries(P.cometOrbitFields({ type: "elliptical", kepler }));
  assert.equal(closed.Orbit, "Bound, returns every 12.3 years");
  assert.equal(closed.Perihelion, "2.00 AU");
  const open = Object.fromEntries(P.cometOrbitFields({ type: "parabolic", kepler: { ...kepler, eccentricity: 0.9987, period_years: null } }));
  assert.match(open.Orbit, /^Unbound.*does not return/);
  assert.equal(open.Eccentricity, "0.9987");
});

test("a parabolic comet is gone once its pass is over, an elliptical one never is (GEN.182)", () => {
  const kepler = {
    perihelion_distance_km: 1 * AU_KM, eccentricity: 0.999, inclination_deg: 0, arg_periapsis_deg: 0,
    ascending_node_deg: 0, primary_mass_solar: 1, mean_anomaly_deg: 10, period_years: 1000,
    parabolic_mean_anomaly: 0,
  };
  const comet = (ref, type) => ({ ref, orbit: { around: "barycenter", type, kepler } });
  const scene = { stars: [], planets: [], comets: [comet("comet:1", "parabolic"), comet("comet:2", "elliptical")] };
  const near = P.relativeAt(scene, 0);
  assert.equal(near["comet:1"].gone, false);
  const later = P.relativeAt(scene, 1e6);
  assert.equal(later["comet:1"].gone, true);
  assert.equal(later["comet:2"].gone, false);
});
