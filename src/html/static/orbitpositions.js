// html/static/orbitpositions.js
//
// Where every body of a system is at any time (MAP.70): the browser twin of
// planetgen/physics/body_positions.py, over the scene
// GET /api/systems/<id>/scene returns. positionsAt(scene, years) moves each
// body `years` past the scene's epoch along its stored orbit and answers
// {ref: [x, y, z]} in kilometres: planets, moons and a pair's second star on
// circular orbits (phase plus 360 degrees a period), comets on Kepler or
// Barker orbits, a close pair either side of the barycenter by the other's
// mass fraction. tests/test_js_unit.py checks the two copies against each other.

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
export function solveEccentricAnomaly(meanAnomaly, eccentricity, tolerance = 1e-10) {
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
  const phase = orbit.phase_deg + (orbit.period_years ? (360 * years) / orbit.period_years : 0);
  const p = orbitalPosition(orbit.distance_km / AU_KM, orbit.inclination_deg, orbit.ascending_node_deg, ((phase % 360) + 360) % 360);
  return [p[0] * AU_KM, p[1] * AU_KM, p[2] * AU_KM];
}

function cometRelativeKm(orbit, years) {
  const k = orbit.kepler;
  const q = k.perihelion_distance_km / AU_KM;
  let trueAnomaly;
  let distance;
  if (orbit.type === "elliptical") {
    const anomaly = rad(k.mean_anomaly_deg + (360 * years) / k.period_years);
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

const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];

// {ref: [x, y, z]} km, `years` after the scene's epoch (negative: before).
// Belts have no single point and are left out.
export function positionsAt(scene, years) {
  const out = {};
  for (const star of scene.stars) {
    if (!star.orbit) out[star.ref] = [0, 0, 0];
  }
  for (const star of scene.stars) {
    if (!star.orbit) continue;
    const relative = relativeCircularKm(star.orbit, years);
    if (star.orbit.around === "barycenter") {
      const f = star.orbit.secondary_mass_fraction;
      out[star.ref] = relative.map((c) => c * (1 - f));
      for (const other of scene.stars) {
        if (!other.orbit) out[other.ref] = relative.map((c) => -c * f);
      }
    } else {
      out[star.ref] = add(out[star.orbit.around], relative);
    }
  }
  const origin = (around) => (around === "barycenter" ? [0, 0, 0] : out[around]);
  for (const planet of scene.planets) {
    out[planet.ref] = add(origin(planet.orbit.around), relativeCircularKm(planet.orbit, years));
    for (const moon of planet.moons) {
      out[moon.ref] = add(out[planet.ref], relativeCircularKm(moon.orbit, years));
    }
  }
  for (const comet of scene.comets) {
    out[comet.ref] = add(origin(comet.orbit.around), cometRelativeKm(comet.orbit, years));
  }
  return out;
}
