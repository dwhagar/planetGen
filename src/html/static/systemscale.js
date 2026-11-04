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

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { convert } = await import(`./unitregistry.js${VERSION_QUERY}`);

export const MODE_TRUE = "true";
export const MODE_COMPRESSED = "compressed";
export const UNIT_KM = 1e6; // true scale: one world unit is a million km.

const SUN_RADIUS_KM = convert(1, "Rsun", "km");
const STAR_RADIUS_UNITS = 5; // a sun-sized star, compressed
const COMPRESSED_INNER_UNITS = 14; // the innermost orbit around the middle
const COMPRESSED_SPAN_UNITS = 150; // inner to outer orbit
const MOON_GAP_UNITS = 0.8;
const MOON_SPAN_UNITS = 6;
const MAX_CLEARED_ECCENTRICITY = 0.8; // the most eccentric node orbit drawn wide enough to clear at periapsis
const NODE_MARGIN = 1.15; // room between the two sides of a tree node's orbit, as a share of what each side needs

const NOTES = {
  [MODE_TRUE]: "True scale: distances and sizes are real, so most bodies are too small to see until you zoom in.",
  [MODE_COMPRESSED]: "Compressed scale: distances are on a log scale and bodies are drawn larger than life.",
};

const nodeRef = (nodeNo) => `node:${nodeNo}`; // orbitpositions.js's `nodeRef`: the ref of a tree node's centre of mass

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
    // A close pair's separation sits among the orbits round the middle; a wide pair's stars are placed by the tree (MAP.213).
    if (star.orbit && star.orbit.around === "barycenter") add(star.orbit.around, star.orbit.distance_km);
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

// Which parents have something going round them: planets, belts and comets round a star.
function orbitedRefs(scene) {
  const out = new Set();
  for (const planet of scene.planets) out.add(planet.orbit.around);
  for (const belt of scene.belts || []) out.add(belt.around);
  for (const comet of scene.comets) out.add(comet.orbit.around);
  return out;
}

// The drawn scale of each non-group tree node, compressed: Map(nodeRef -> units per km). Children are sized before their
// parents (a node's reach is how far its farthest star, with what goes round it, lies from its centre).
function nodeScales(scene, radii) {
  const tree = scene.tree;
  const orbited = orbitedRefs(scene);
  const byParent = new Map();
  for (const member of tree.members) {
    const list = byParent.get(member.parent_node_no) || [];
    list[member.slot] = { reach: orbited.has(member.ref) ? COMPRESSED_INNER_UNITS + COMPRESSED_SPAN_UNITS : radii[member.ref] || 0,
      mass: member.mass_solar };
    byParent.set(member.parent_node_no, list);
  }
  const scales = new Map();
  const reaches = new Map();
  const children = new Map();
  for (const node of tree.nodes) {
    if (node.parent_node_no == null) continue;
    if (!children.has(node.parent_node_no)) children.set(node.parent_node_no, []);
    children.get(node.parent_node_no).push(node);
  }
  const size = (node) => {
    for (const child of children.get(node.node_no) || []) size(child);
    const sides = byParent.get(node.node_no) || [];
    for (const child of children.get(node.node_no) || []) sides[child.slot] = { reach: reaches.get(child.node_no), mass: child.mass_solar };
    if (node.kind === "group" || !sides[0] || !sides[1]) {
      reaches.set(node.node_no, Math.max(0, ...sides.map((side) => (side ? side.reach : 0))));
      return;
    }
    // (the two sides are closest at periapsis, a(1 - e) apart, so an eccentric orbit is drawn wider by that much, up to e = 0.8)
    const periapsisShare = 1 - Math.min(node.eccentricity || 0, MAX_CLEARED_ECCENTRICITY);
    const units = Math.max(COMPRESSED_INNER_UNITS, (NODE_MARGIN * (sides[0].reach + sides[1].reach)) / periapsisShare);
    const total = sides[0].mass + sides[1].mass;
    let reach = Math.max(units * (sides[1].mass / total) + sides[0].reach, units * (sides[0].mass / total) + sides[1].reach);
    // A pair with bodies going round its centre (a close pair below a further star's orbit) is as wide as the span they cover
    // beyond its stars' farthest reach from the centre, a(1 + e).
    if (orbited.has(nodeRef(node.node_no))) reach = Math.max(reach, units * (1 + (node.eccentricity || 0)) + COMPRESSED_SPAN_UNITS);
    reaches.set(node.node_no, reach);
    if (node.separation_km > 0) scales.set(nodeRef(node.node_no), units / node.separation_km);
  };
  for (const node of tree.nodes) {
    if (node.parent_node_no == null) size(node);
  }
  return scales;
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
  // MAP.213: a tree node's orbit is drawn on its own scale, so a 10,000 AU outer pair and a 0.1 AU inner pair both read: the
  // node's separation is as many units as its two sides need to clear each other (a star with bodies going round it needs
  // their whole span, a star alone its drawn size, a sub-pair the reach of its own orbit), and the orbit is that, scaled
  // straight, so its ellipse keeps its shape.
  // The bodies going round a node's centre (`around` is its `nodeRef`) start outside the node's own orbit, where the stars
  // of the node are never found: beyond a(1 + e) the map compresses them over the span like any other parent's.
  if (compressed && scene.tree) {
    const nodes = new Map(scene.tree.nodes.map((node) => [nodeRef(node.node_no), node]));
    for (const [ref, scale] of nodeScales(scene, radii)) {
      const edgeKm = nodes.get(ref).separation_km * (1 + (nodes.get(ref).eccentricity || 0));
      const list = distances[ref];
      const hi = list ? Math.max(...list) : 0;
      maps[ref] = hi > edgeKm
        ? (km) => (km <= edgeKm ? km * scale : edgeKm * scale + (COMPRESSED_SPAN_UNITS * Math.log(km / edgeKm)) / Math.log(hi / edgeKm))
        : (km) => km * scale;
    }
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
