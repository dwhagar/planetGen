// html/static/systemscale.js
//
// How the 3D system view maps real kilometres to what is drawn (MAP.71).
// True scale makes planets invisible dots, so the view has two modes:
//   "true"        every distance and radius in km over UNIT_KM: nothing
//                 shrunk or grown, the honest picture;
//   "compressed"  distances on the System Map's log scale and bodies
//                 enlarged (also on a log scale), so every orbit and
//                 every body stays visible; a moon's orbit stays outside
//                 its planet's drawn size.
// A layout answers three things for a scene: where a body whose position
// relative to its parent is `rel` km is drawn (`place`), how large it is
// drawn (`radiusOf`), and the note saying which mode is shown (`note`).
// Plain data in, plain data out: no three.js here.

export const MODE_TRUE = "true";
export const MODE_COMPRESSED = "compressed";
export const UNIT_KM = 1e6; // true scale: one world unit is a million km.

const SUN_RADIUS_KM = 695700;
const STAR_RADIUS_UNITS = 5; // a sun-sized star, compressed
const COMPRESSED_INNER_UNITS = 14; // the innermost orbit around the middle
const COMPRESSED_SPAN_UNITS = 150; // inner to outer orbit
const MOON_GAP_UNITS = 0.8;
const MOON_SPAN_UNITS = 6;

const NOTES = {
  [MODE_TRUE]: "True scale: distances and sizes are real, so most bodies are too small to see until you zoom in.",
  [MODE_COMPRESSED]: "Compressed scale: distances are on a log scale and bodies are drawn larger than life.",
};

const length = (v) => Math.hypot(v[0], v[1], v[2]);

// Which distances a parent's children sit at: orbit radii of planets and
// belts (inner and outer edge), a comet's perihelion, a second star's
// separation. Moons and planets group by `around`.
function groupDistances(scene) {
  const groups = {};
  const add = (around, km) => {
    if (km > 0) (groups[around] = groups[around] || []).push(km);
  };
  for (const star of scene.stars) {
    if (star.orbit) add(star.orbit.around, star.orbit.distance_km);
  }
  for (const planet of scene.planets) {
    add(planet.orbit.around, planet.orbit.distance_km);
    for (const moon of planet.moons) add(planet.ref, moon.orbit.distance_km);
  }
  for (const belt of scene.belts) {
    add(belt.around, belt.inner_km);
    add(belt.around, belt.outer_km);
  }
  for (const comet of scene.comets) add(comet.orbit.around, comet.orbit.kepler.perihelion_distance_km);
  return groups;
}

export function radiusUnitsCompressed(kind, radiusKm) {
  if (kind === "star") return STAR_RADIUS_UNITS * Math.pow(Math.max(radiusKm, 1) / SUN_RADIUS_KM, 0.35);
  if (kind === "comet") return 0.18;
  if (kind === "moon") return 0.12 + 0.12 * Math.log10(1 + radiusKm / 100);
  return 0.35 + 0.45 * Math.log10(1 + radiusKm / 500); // planet
}

// A layout of `scene` in `mode`.
export function createLayout(scene, mode) {
  const compressed = mode === MODE_COMPRESSED;
  const bodies = {};
  for (const star of scene.stars) bodies[star.ref] = star;
  for (const planet of scene.planets) {
    bodies[planet.ref] = planet;
    for (const moon of planet.moons) bodies[moon.ref] = moon;
  }
  for (const comet of scene.comets) bodies[comet.ref] = comet;

  const radii = {};
  for (const ref of Object.keys(bodies)) {
    const b = bodies[ref];
    radii[ref] = compressed ? radiusUnitsCompressed(b.kind, b.radius_km) : b.radius_km / UNIT_KM;
  }

  // Per parent: the compressed map r(km) -> units. Planets and the like
  // around the middle (or a star) spread over the span; a planet's moons
  // sit just outside its drawn size.
  const distances = groupDistances(scene);
  const maps = {};
  for (const around of Object.keys(distances)) {
    const list = distances[around];
    const lo = Math.min(...list);
    const hi = Math.max(...list);
    const isPlanet = around !== "barycenter" && !around.startsWith("star:");
    const inner = isPlanet ? radii[around] + MOON_GAP_UNITS : COMPRESSED_INNER_UNITS;
    const span = isPlanet ? MOON_SPAN_UNITS : COMPRESSED_SPAN_UNITS;
    const logSpan = Math.log(hi / lo);
    maps[around] = (km) => {
      if (km <= 0) return 0;
      if (km <= lo) return inner * (km / lo);
      return inner + (logSpan > 0 ? (span * Math.log(km / lo)) / logSpan : 0);
    };
  }
  const fallback = (km) => km / UNIT_KM;

  return {
    mode: mode,
    note: NOTES[compressed ? MODE_COMPRESSED : MODE_TRUE],
    radiusOf: (ref) => radii[ref],
    // The drawn offset of a body `rel` km from `around`, in world units.
    place: function (around, rel) {
      if (!compressed) return [rel[0] / UNIT_KM, rel[1] / UNIT_KM, rel[2] / UNIT_KM];
      const d = length(rel);
      if (d === 0) return [0, 0, 0];
      const map = maps[around] || fallback;
      const scale = map(d) / d;
      return [rel[0] * scale, rel[1] * scale, rel[2] * scale];
    },
    // A distance from the middle drawn (a heliopause, a belt's edge).
    distance: function (around, km) {
      if (!compressed) return km / UNIT_KM;
      return (maps[around] || fallback)(km);
    },
  };
}

// World positions of every body (the offsets summed up their parents) from
// relativeAt()'s answer: {ref: [x, y, z]} in world units.
export function layoutPositions(layout, relative) {
  const out = {};
  const resolve = (ref) => {
    if (out[ref]) return out[ref];
    const { around, rel } = relative[ref];
    const offset = layout.place(around, rel);
    if (around === "barycenter") {
      out[ref] = offset;
    } else {
      const base = resolve(around);
      out[ref] = [base[0] + offset[0], base[1] + offset[1], base[2] + offset[2]];
    }
    return out[ref];
  };
  for (const ref of Object.keys(relative)) resolve(ref);
  return out;
}
