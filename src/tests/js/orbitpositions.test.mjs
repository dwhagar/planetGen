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

test("positionsAt places the stars of a wide pair and of a triple where the tree's eccentric orbits put them (MAP.255)", () => {
  for (const [name, { scene, cases }] of Object.entries(F.orbitPositions.trees)) {
    for (const c of cases) {
      const got = P.positionsAt(scene, c.years);
      assert.deepEqual(Object.keys(got).sort(), Object.keys(c.positions).sort(), name);
      for (const ref of Object.keys(c.positions)) {
        c.positions[ref].forEach((expected, axis) => {
          const tolerance = Math.max(1, Math.abs(expected) * 1e-9);
          assert.ok(Math.abs(got[ref][axis] - expected) <= tolerance,
            `${name}: ${ref} axis ${axis} at ${c.years} y: ${got[ref][axis]} vs ${expected}`);
        });
      }
    }
  }
});

test("a wide pair's stars are closest at periapsis and farthest half a period later (MAP.255)", () => {
  const { scene } = F.orbitPositions.trees.wide;
  const node = scene.tree.nodes[0];
  const periodsToPeriapsis = -node.mean_anomaly_deg / 360 * node.period_years;
  const apart = (years) => {
    const at = P.positionsAt(scene, years);
    return Math.hypot(...at["star:2"].map((c, axis) => c - at["star:1"][axis]));
  };
  assert.ok(Math.abs(apart(periodsToPeriapsis) - 30 * 0.4 * AU_KM) < 1e-6 * AU_KM);
  assert.ok(Math.abs(apart(periodsToPeriapsis + node.period_years / 2) - 30 * 1.6 * AU_KM) < 1e-6 * AU_KM);
});

test("relativeAt hangs each star from the centre of its node and each node from its parent's, as virtual points (MAP.213)", () => {
  const { scene } = F.orbitPositions.trees.triple;
  for (const years of [0, 1.8, 400]) {
    const at = P.relativeAt(scene, years);
    assert.deepEqual(Object.keys(at).filter((ref) => at[ref].virtual).sort(), ["node:0", "node:1"]);
    assert.deepEqual([at["node:1"].around, at["node:0"].around], ["barycenter", "node:1"]);
    assert.deepEqual(at["node:1"].rel, [0, 0, 0]);
    assert.deepEqual([at["star:1"].around, at["star:2"].around, at["star:3"].around], ["node:0", "node:0", "node:1"]);
    // The chain sums to positionsAt's positions, which are Python's.
    const placed = P.positionsAt(scene, years);
    const chain = (ref) => {
      const { around, rel } = at[ref];
      return around === "barycenter" ? rel : chain(around).map((c, axis) => c + rel[axis]);
    };
    for (const ref of ["star:1", "star:2", "star:3"]) {
      chain(ref).forEach((c, axis) => assert.ok(Math.abs(c - placed[ref][axis]) <= Math.max(1, Math.abs(c) * 1e-9), `${ref} ${axis}`));
    }
    assert.equal(Object.keys(placed).some((ref) => ref.startsWith("node:")), false);
  }
});

test("a node's paths are its two children's ellipses about its centre, each at the other's share of the separation (MAP.213)", () => {
  const { scene } = F.orbitPositions.trees.triple;
  const paths = P.nodePaths(scene, 64);
  assert.deepEqual(paths.map((path) => [path.ref, path.around]).sort(), [
    ["node:0", "node:1"], ["star:1", "node:0"], ["star:2", "node:0"], ["star:3", "node:1"]]);
  const node = scene.tree.nodes[0];
  const star1 = paths.find((path) => path.ref === "star:1");
  const star2 = paths.find((path) => path.ref === "star:2");
  assert.equal(star1.points.length, 65);
  // At periapsis the pair is a(1 - e) apart, star 2 on 0.25 of it (a quarter of the pair's mass is star 2).
  const gap = Math.hypot(...star2.points[0].map((c, axis) => c - star1.points[0][axis]));
  assert.ok(Math.abs(gap - node.separation_km * (1 - node.eccentricity)) < 1e-6 * node.separation_km);
  assert.ok(Math.abs(Math.hypot(...star2.points[0]) - 0.75 * gap) < 1e-6 * gap);
  const at = P.relativeAt(scene, 0);
  const near = (point, target) => Math.hypot(...point.map((c, axis) => c - target[axis]));
  assert.ok(star2.points.some((point) => near(point, at["star:2"].rel) < 0.2 * node.separation_km));
});

test("a scene with stars and no tree of orbits is refused (MAP.255)", () => {
  assert.throws(() => P.relativeAt({ stars: [{ ref: "star:1", orbit: null }], planets: [], comets: [] }, 0), /no tree of orbits/);
});

test("wrappedTurnDeg matches kepler.wrapped_turn_deg, whole orbits dropped before the division", () => {
  for (const c of F.orbitPositions.turns) {
    assert.ok(Math.abs(P.wrappedTurnDeg(c.elapsed, c.period) - c.degrees) <= 1e-12 * 360,
      `${c.elapsed} over ${c.period}: ${P.wrappedTurnDeg(c.elapsed, c.period)} vs ${c.degrees}`);
  }
});

test("stumpff matches kepler.stumpff_c2_c3 inside and outside the series limit", () => {
  for (const c of F.orbitPositions.stumpff) {
    const got = P.stumpff(c.z);
    c.c.forEach((expected, n) => {
      assert.ok(Math.abs(got[n] - expected) <= 1e-15 * Math.abs(expected) + 1e-300, `c${n + 2}(${c.z}): ${got[n]} vs ${expected}`);
    });
  }
});

test("universalStep matches kepler.universal_step for every size and sign of step", () => {
  for (const c of F.orbitPositions.steps) {
    const { position, velocity } = P.universalStep(c.position, c.velocity, c.mu, c.dt);
    c.position_after.forEach((expected, axis) => {
      assert.ok(Math.abs(position[axis] - expected) <= Math.max(1, Math.abs(expected) * 1e-9),
        `position ${axis} after ${c.dt} s: ${position[axis]} vs ${expected}`);
    });
    c.velocity_after.forEach((expected, axis) => {
      assert.ok(Math.abs(velocity[axis] - expected) <= Math.max(1e-6, Math.abs(expected) * 1e-9),
        `velocity ${axis} after ${c.dt} s: ${velocity[axis]} vs ${expected}`);
    });
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
  assert.equal(closed.Perihelion, "2 AU");
  const open = Object.fromEntries(P.cometOrbitFields({ type: "parabolic", kepler: { ...kepler, eccentricity: 0.9987, period_years: null } }));
  assert.match(open.Orbit, /^Unbound.*does not return/);
  assert.equal(open.Eccentricity, "0.9987 (a ratio: 0 is a circle, 1 is a parabola)");
  assert.equal(closed.Eccentricity, "0.500 (a ratio: 0 is a circle, 1 is a parabola)");
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
