// html/static/orbitpositions.js
//
// Where every body of a system is at any time (MAP.70): the browser twin of
// planetgen/physics/body_positions.py, over the scene
// GET /api/systems/<id>/scene returns. positionsAt(scene, years) moves each
// body `years` past the scene's epoch along its stored orbit and answers
// {ref: [x, y, z]} in kilometres: planets and moons by their stored state
// vectors moved with the universal-variable Kepler step (GEN.289), a pair's
// second star on a circular orbit (phase plus 360 degrees a period), comets on
// Kepler or Barker orbits, a close pair either side of the barycenter by the
// other's mass fraction. tests/test_js_unit.py checks the two copies against
// each other.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { formatDistanceKm, formatPeriodYears } = await import(`./unitladder.js${VERSION_QUERY}`);

const AU_KM = 149597870.7;
const TWO_PI = 2 * Math.PI;
const MU_AU3_YR2_PER_SOLAR_MASS = 4 * Math.PI * Math.PI;
const rad = (deg) => (deg * Math.PI) / 180;
const deg = (radians) => (radians * 180) / Math.PI;

// A circular orbit's position (orbits.orbital_position_au), same units as distance.
export function orbitalPosition(distance, inclinationDeg, ascendingNodeDeg, phaseDeg) {
  const u = rad(phaseDeg);
  const i = rad(inclinationDeg);
  const node = rad(ascendingNodeDeg);
  const cosU = Math.cos(u);
  const sinU = Math.sin(u);
  const cosI = Math.cos(i);
  const cosNode = Math.cos(node);
  const sinNode = Math.sin(node);
  return [
    distance * (cosNode * cosU - sinNode * sinU * cosI),
    distance * (sinNode * cosU + cosNode * sinU * cosI),
    distance * sinU * Math.sin(i),
  ];
}

// Kepler's equation M = E - e sin E (kepler.solve_eccentric_anomaly).
export function solveEccentricAnomaly(meanAnomaly, eccentricity, tolerance = 1e-12) {
  const m = ((meanAnomaly % TWO_PI) + TWO_PI) % TWO_PI;
  let e = eccentricity < 0.8 ? m : Math.PI;
  for (let n = 0; n < 100; n += 1) {
    const delta = (e - eccentricity * Math.sin(e) - m) / (1 - eccentricity * Math.cos(e));
    e -= delta;
    if (Math.abs(delta) < tolerance) return e;
  }
  let lo = 0;
  let hi = TWO_PI;
  while (hi - lo > tolerance) {
    const mid = (lo + hi) / 2;
    if (mid === lo || mid === hi) break;
    if (mid - eccentricity * Math.sin(mid) - m < 0) lo = mid;
    else hi = mid;
  }
  return (lo + hi) / 2;
}

const cubeRoot = (value) => Math.sign(value) * Math.pow(Math.abs(value), 1 / 3);

// Barker's equation D^3 + 3D = 3 Mp, D = tan(true anomaly / 2).
export function solveBarker(parabolicMeanAnomaly) {
  const w = 1.5 * parabolicMeanAnomaly;
  const s = Math.abs(w) < 1e150 ? Math.sqrt(w * w + 1) : Math.abs(w);
  return cubeRoot(w + s) + cubeRoot(w - s);
}

function relativeCircularKm(orbit, years) {
  const phase = orbit.phase_deg + wrappedTurnDeg(years, orbit.period_years);
  const p = orbitalPosition(orbit.distance_km / AU_KM, orbit.inclination_deg, orbit.ascending_node_deg, ((phase % 360) + 360) % 360);
  return [p[0] * AU_KM, p[1] * AU_KM, p[2] * AU_KM];
}

// The degrees an orbit of `period` turns in `elapsed`, whole orbits dropped first so a long step keeps its fraction
// (kepler.wrapped_turn_deg; `%` is C's fmod, exact).
export function wrappedTurnDeg(elapsed, period) {
  return period ? 360 * ((elapsed % period) / period) : 0;
}

// The Stumpff functions c2 and c3 (kepler.stumpff_c2_c3): ten terms of the series below |z| = 1, the closed forms beyond.
const STUMPFF_SERIES_LIMIT = 1;
const C2_COEFFICIENTS = [];
const C3_COEFFICIENTS = [];
for (let k = 0, f = 2; k < 10; k += 1, f += 2) {
  let factorial2 = 1;
  for (let n = 2; n <= f; n += 1) factorial2 *= n;
  C2_COEFFICIENTS.push(1 / factorial2);
  C3_COEFFICIENTS.push(1 / (factorial2 * (f + 1)));
}

export function stumpff(z) {
  if (Math.abs(z) < STUMPFF_SERIES_LIMIT) {
    let c2 = 0;
    let c3 = 0;
    for (let k = C2_COEFFICIENTS.length - 1; k >= 0; k -= 1) {
      c2 = c2 * -z + C2_COEFFICIENTS[k];
      c3 = c3 * -z + C3_COEFFICIENTS[k];
    }
    return [c2, c3];
  }
  if (z > 0) {
    const root = Math.sqrt(z);
    return [(1 - Math.cos(root)) / z, (root - Math.sin(root)) / root ** 3];
  }
  const root = Math.sqrt(-z);
  return [(Math.cosh(root) - 1) / -z, (Math.sinh(root) - root) / root ** 3];
}

const UNIVERSAL_TOLERANCE = 1e-13;
const UNIVERSAL_MAX_ITERATIONS = 60;

// kepler.universal_step: a two-body state `dt` on, whatever the conic, in any
// consistent units. Halley's method on the universal anomaly kept inside a
// bracket that holds the root, bisecting when a step leaves it; a whole number
// of orbits is dropped first. Returns {position, velocity} (arrays of three)
// and throws on a bad input or no convergence, as the Python does.
export function universalStep(position, velocity, mu, dt) {
  if (!(mu > 0) || !Number.isFinite(mu)) throw new Error(`universalStep: mu must be positive, got ${mu}`);
  if (![...position, ...velocity, dt].every(Number.isFinite)) throw new Error("universalStep: the state and step must be finite");
  if (dt === 0) return { position: [...position], velocity: [...velocity] };
  const r0 = Math.hypot(...position);
  if (r0 === 0) throw new Error("universalStep: the body sits on its primary");
  const v0Sq = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
  const sqrtMu = Math.sqrt(mu);
  const s0 = (position[0] * velocity[0] + position[1] * velocity[1] + position[2] * velocity[2]) / sqrtMu;
  const alpha = 2 / r0 - v0Sq / mu;
  let step = dt;

  const residual = (chi) => {
    const z = alpha * chi * chi;
    const [c2, c3] = stumpff(z);
    const f = s0 * chi * chi * c2 + (1 - alpha * r0) * chi ** 3 * c3 + r0 * chi - sqrtMu * step;
    const df = s0 * chi * (1 - z * c3) + (1 - alpha * r0) * chi * chi * c2 + r0;
    const d2f = s0 * (1 - z * c2) + (1 - alpha * r0) * chi * (1 - z * c3);
    return [f, df, d2f];
  };

  let chi;
  if (alpha > 1e-12 / r0) {
    chi = sqrtMu * step * alpha;
    const period = (2 * Math.PI) / (sqrtMu * alpha ** 1.5);
    if (Math.abs(step) > period) {
      step %= period;
      chi = sqrtMu * step * alpha;
    }
  } else if (alpha < -1e-12 / r0) {
    const a = 1 / alpha;
    const sign = step > 0 ? 1 : -1;
    const inner = (-2 * mu * alpha * step) / (s0 * sqrtMu + sign * Math.sqrt(-mu * a) * (1 - r0 * alpha));
    chi = inner > 0 ? sign * Math.sqrt(-a) * Math.log(inner) : sign * Math.sqrt(-a);
  } else {
    chi = (sqrtMu * step) / r0;
  }

  const guess = Math.abs(chi) || 1;
  let [low, high] = step > 0 ? [0, guess] : [-guess, 0];
  for (let n = 0; n < 200; n += 1) {
    if (step > 0 && residual(high)[0] < 0) [low, high] = [high, high * 2];
    else if (step < 0 && residual(low)[0] > 0) [low, high] = [low * 2, low];
    else break;
  }
  if (!(low <= chi && chi <= high)) chi = 0.5 * (low + high);

  let converged = false;
  for (let n = 0; n < UNIVERSAL_MAX_ITERATIONS; n += 1) {
    const [f, df, d2f] = residual(chi);
    if (f > 0) high = chi;
    else low = chi;
    const denominator = 2 * df * df - f * d2f;
    let next = chi - (denominator !== 0 ? (2 * f * df) / denominator : f / df);
    if (!(low < next && next < high)) next = 0.5 * (low + high);
    const done = Math.abs(next - chi) <= UNIVERSAL_TOLERANCE * Math.max(1, Math.abs(next))
      || high - low <= 4e-16 * Math.max(1, Math.abs(next));
    chi = next;
    if (done) { converged = true; break; }
  }
  if (!converged) throw new Error(`universalStep: no convergence in ${UNIVERSAL_MAX_ITERATIONS} iterations`);

  const z = alpha * chi * chi;
  const [c2, c3] = stumpff(z);
  const f = 1 - (chi * chi * c2) / r0;
  const g = step - (chi ** 3 * c3) / sqrtMu;
  const r = position.map((a, i) => f * a + g * velocity[i]);
  const rLength = Math.hypot(...r);
  const gDot = 1 - (chi * chi * c2) / rLength;
  const fDot = (sqrtMu * chi * (z * c3 - 1)) / (rLength * r0);
  return { position: r, velocity: position.map((a, i) => fDot * a + gDot * velocity[i]) };
}

const SECONDS_PER_YEAR = 365.25 * 86400;

// A planet's or moon's position relative to what it goes round, km (GEN.289):
// its stored state vector moved from the epoch it is true at to `years` past
// the scene's epoch (body_positions._relative_state_km). The pull is
// mu = (2 pi / period)^2 distance^3, in km and seconds.
function relativeStateKm(orbit, years, sceneEpochUnix) {
  const offset = orbit.epoch_unix == null || sceneEpochUnix == null ? 0 : sceneEpochUnix - orbit.epoch_unix;
  const omega = TWO_PI / (orbit.period_years * SECONDS_PER_YEAR);
  const mu = omega * omega * orbit.distance_km ** 3;
  return universalStep(orbit.position_km, orbit.velocity_kms, mu, years * SECONDS_PER_YEAR + offset).position;
}

function cometRelativeKm(orbit, years) {
  const k = orbit.kepler;
  const q = k.perihelion_distance_km / AU_KM;
  let trueAnomaly;
  let distance;
  if (orbit.type === "elliptical") {
    const anomaly = rad(k.mean_anomaly_deg + wrappedTurnDeg(years, k.period_years));
    const a = q / (1 - k.eccentricity);
    const ecc = solveEccentricAnomaly(anomaly, k.eccentricity);
    trueAnomaly = 2 * Math.atan2(
      Math.sqrt(1 + k.eccentricity) * Math.sin(ecc / 2),
      Math.sqrt(1 - k.eccentricity) * Math.cos(ecc / 2));
    distance = a * (1 - k.eccentricity * Math.cos(ecc));
  } else {
    const mu = MU_AU3_YR2_PER_SOLAR_MASS * k.primary_mass_solar;
    const anomaly = k.parabolic_mean_anomaly + Math.sqrt(mu / (2 * q ** 3)) * years;
    const d = solveBarker(anomaly);
    trueAnomaly = 2 * Math.atan(d);
    distance = q * (1 + d * d);
  }
  const latitude = deg(rad(k.arg_periapsis_deg) + trueAnomaly) % 360;
  const p = orbitalPosition(distance, k.inclination_deg, k.ascending_node_deg, latitude);
  return [p[0] * AU_KM, p[1] * AU_KM, p[2] * AU_KM];
}

// GEN.182: a comet on an open (parabolic) orbit makes one pass and leaves. The
// drawn path reaches out to this distance; a comet farther out than that is
// off its path (not yet arrived, or gone for good) and is not shown.
export const COMET_PATH_LIMIT_AU = 2000;

const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];

// Where every body is, `years` after the scene's epoch, as a position
// relative to what it goes round: {ref: {around, rel: [x, y, z] km}}.
// `around` is "barycenter" or another body's ref; a close pair's stars sit
// either side of the barycenter. Belts have no single point and are left out.
export function relativeAt(scene, years) {
  const out = {};
  for (const star of scene.stars) {
    if (!star.orbit) out[star.ref] = { around: "barycenter", rel: [0, 0, 0] };
  }
  for (const star of scene.stars) {
    if (!star.orbit) continue;
    const relative = relativeCircularKm(star.orbit, years);
    if (star.orbit.around === "barycenter") {
      const f = star.orbit.secondary_mass_fraction;
      out[star.ref] = { around: "barycenter", rel: relative.map((c) => c * (1 - f)) };
      for (const other of scene.stars) {
        if (!other.orbit) out[other.ref] = { around: "barycenter", rel: relative.map((c) => -c * f) };
      }
    } else {
      out[star.ref] = { around: star.orbit.around, rel: relative };
    }
  }
  for (const planet of scene.planets) {
    out[planet.ref] = { around: planet.orbit.around, rel: relativeStateKm(planet.orbit, years, scene.epoch_unix) };
    for (const moon of planet.moons) {
      out[moon.ref] = { around: planet.ref, rel: relativeStateKm(moon.orbit, years, scene.epoch_unix) };
    }
  }
  for (const comet of scene.comets) {
    const rel = cometRelativeKm(comet.orbit, years);
    const gone = comet.orbit.type !== "elliptical" && Math.hypot(rel[0], rel[1], rel[2]) > COMET_PATH_LIMIT_AU * AU_KM;
    out[comet.ref] = { around: comet.orbit.around, rel: rel, gone: gone };
  }
  return out;
}

// {ref: [x, y, z]} km, `years` after the scene's epoch (negative: before),
// in the scene's frame.
export function positionsAt(scene, years) {
  const relative = relativeAt(scene, years);
  const out = {};
  const resolve = (ref) => {
    if (out[ref]) return out[ref];
    const { around, rel } = relative[ref];
    out[ref] = around === "barycenter" ? rel : add(resolve(around), rel);
    return out[ref];
  };
  for (const ref of Object.keys(relative)) resolve(ref);
  return out;
}

// The points of one orbit's path, relative to what it goes round, in km: a
// circular orbit's circle, a comet's ellipse or the near part of its
// parabola (out to `maxAu`). `samples` points, the first repeated at the
// end of a closed path.
export function orbitPath(orbit, samples = 128, maxAu = COMET_PATH_LIMIT_AU) {
  const out = [];
  if (!orbit.kepler) {
    const r = orbit.distance_km / AU_KM;
    for (let n = 0; n <= samples; n += 1) {
      const p = orbitalPosition(r, orbit.inclination_deg, orbit.ascending_node_deg, (360 * n) / samples);
      out.push([p[0] * AU_KM, p[1] * AU_KM, p[2] * AU_KM]);
    }
    return out;
  }
  const k = orbit.kepler;
  const q = k.perihelion_distance_km / AU_KM;
  const e = orbit.type === "elliptical" ? k.eccentricity : 1;
  const semiLatus = q * (1 + e);
  let limit = Math.PI;
  if (e >= 1) limit = Math.acos(Math.max(-1, Math.min(1, semiLatus / maxAu - 1))) || 0.1;
  for (let n = 0; n <= samples; n += 1) {
    const nu = -limit + (2 * limit * n) / samples;
    const r = semiLatus / (1 + e * Math.cos(nu));
    const p = orbitalPosition(r, k.inclination_deg, k.ascending_node_deg, deg(rad(k.arg_periapsis_deg) + nu));
    out.push([p[0] * AU_KM, p[1] * AU_KM, p[2] * AU_KM]);
  }
  return out;
}

// What the info panel says about a comet's orbit (GEN.182): a closed orbit
// returns every period; a parabolic one is unbound, passes the star once and
// does not come back, so its drawn path is open at both ends. `orbit` is a
// comet's scene `orbit`; returns [label, text] pairs.
export function cometOrbitFields(orbit) {
  const k = orbit.kepler;
  const fields = [];
  if (orbit.type === "elliptical") {
    const period = k.period_years;
    fields.push(["Orbit", "Bound, returns every " + formatPeriodYears(period)]);
  } else {
    fields.push(["Orbit", "Unbound: passes the star once and does not return"]);
  }
  fields.push(["Perihelion", formatDistanceKm(k.perihelion_distance_km)]);
  fields.push(["Eccentricity", k.eccentricity.toFixed(orbit.type === "elliptical" ? 3 : 4)]);
  return fields;
}
