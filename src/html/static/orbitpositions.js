// html/static/orbitpositions.js
//
// Where every body of a system is at any time (MAP.70): the browser twin of
// planetgen/physics/body_positions.py, over the scene
// GET /api/systems/<id>/scene returns. positionsAt(scene, years) moves each
// body `years` past the scene's epoch along its stored orbit and answers
// {ref: [x, y, z]} in kilometres: planets and moons by their stored state
// vectors moved with the universal-variable Kepler step (GEN.289), a pair's
// stars where the scene's tree of orbits puts them (GEN.311, MAP.255: each node
// the Kepler orbit, eccentric when it is, of its slot-1 child about its slot-0
// child, the two children either side of the node's centre by the other's
// mass), comets on Kepler or Barker orbits. tests/test_js_unit.py checks the
// two copies against each other.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { formatDistanceKm, formatPeriodYears } = await import(`./unitladder.js${VERSION_QUERY}`);
const { convert } = await import(`./unitregistry.js${VERSION_QUERY}`);

const AU_KM = convert(1, "au", "km");
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

const SECONDS_PER_YEAR = convert(1, "yr", "s");

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

// The position on an elliptical orbit of perihelion distance `q` and eccentricity `e` at mean anomaly `meanAnomaly`
// (radians), in the units of `q` (kepler.comet_orbital_state, elliptical).
function ellipticalPosition(q, e, meanAnomaly, inclinationDeg, ascendingNodeDeg, argPeriapsisDeg) {
  const a = q / (1 - e);
  const ecc = solveEccentricAnomaly(meanAnomaly, e);
  const trueAnomaly = 2 * Math.atan2(Math.sqrt(1 + e) * Math.sin(ecc / 2), Math.sqrt(1 - e) * Math.cos(ecc / 2));
  const distance = a * (1 - e * Math.cos(ecc));
  const latitude = deg(rad(argPeriapsisDeg) + trueAnomaly) % 360;
  return orbitalPosition(distance, inclinationDeg, ascendingNodeDeg, latitude);
}

function cometRelativeKm(orbit, years) {
  const k = orbit.kepler;
  const q = k.perihelion_distance_km / AU_KM;
  let p;
  if (orbit.type === "elliptical") {
    const anomaly = rad(k.mean_anomaly_deg + wrappedTurnDeg(years, k.period_years));
    p = ellipticalPosition(q, k.eccentricity, anomaly, k.inclination_deg, k.ascending_node_deg, k.arg_periapsis_deg);
  } else {
    const mu = MU_AU3_YR2_PER_SOLAR_MASS * k.primary_mass_solar;
    const anomaly = k.parabolic_mean_anomaly + Math.sqrt(mu / (2 * q ** 3)) * years;
    const d = solveBarker(anomaly);
    const latitude = deg(rad(k.arg_periapsis_deg) + 2 * Math.atan(d)) % 360;
    p = orbitalPosition(q * (1 + d * d), k.inclination_deg, k.ascending_node_deg, latitude);
  }
  return [p[0] * AU_KM, p[1] * AU_KM, p[2] * AU_KM];
}

// The slot-1 child of a tree node relative to its slot-0 child at `unix`, km, in the system's frame
// (system_tree.node_relative_au). The mean anomaly runs from the node's own epoch, a whole number of orbits dropped first.
function nodeRelativeKm(node, unix) {
  const anomaly = node.period_years > 0
    ? node.mean_anomaly_deg + wrappedTurnDeg(unix - node.epoch_unix, node.period_years * SECONDS_PER_YEAR)
    : node.mean_anomaly_deg;
  const periapsis = (node.separation_km / AU_KM) * (1 - node.eccentricity);
  const p = ellipticalPosition(periapsis, node.eccentricity, rad(((anomaly % 360) + 360) % 360),
    node.inclination_deg, node.ascending_node_deg, node.arg_periapsis_deg);
  return [p[0] * AU_KM, p[1] * AU_KM, p[2] * AU_KM];
}

// GEN.182: a comet on an open (parabolic) orbit makes one pass and leaves. The
// drawn path reaches out to this distance; a comet farther out than that is
// off its path (not yet arrived, or gone for good) and is not shown.
export const COMET_PATH_LIMIT_AU = 2000;

const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];

// The scene's tree of orbits as the walk needs it, built once per scene: each node's children by slot (a member or a node),
// the one root, and the ref of the member at the origin (system_tree.member_offsets_au).
const TREES = new WeakMap();

function treeLayout(tree) {
  let layout = TREES.get(tree);
  if (layout) return layout;
  const children = new Map();
  const slotsOf = (nodeNo) => {
    if (!children.has(nodeNo)) children.set(nodeNo, new Map());
    return children.get(nodeNo);
  };
  for (const member of tree.members) slotsOf(member.parent_node_no).set(member.slot, { member });
  for (const node of tree.nodes) {
    if (node.parent_node_no != null) slotsOf(node.parent_node_no).set(node.slot, { node });
  }
  const roots = tree.nodes.filter((node) => node.parent_node_no == null);
  if (tree.nodes.length && roots.length !== 1) throw new Error(`a tree has exactly one root node, found ${roots.length}`);
  layout = { children, root: roots[0] || null };
  TREES.set(tree, layout);
  return layout;
}

// The node centres and the offset of each star from its node's centre, `years` after the scene's epoch, km, in the
// tree's own frame (the root's centre of mass at the origin): {centres: Map(node_no -> [x, y, z]), offsets: {ref: [x, y, z]},
// shares: Map(node_no -> {low, high})}. The tree's nodes hold their phase at their own epoch; the scene's epoch is the time
// zero of `years`, the first node's when the scene has none (body_positions._star_positions_km,
// system_tree.member_offsets_au). A group node's members sit at the node's centre.
function treeAt(scene, years) {
  const tree = scene.tree;
  if (!tree) throw new Error("the scene carries no tree of orbits for its stars");
  const { children, root } = treeLayout(tree);
  const centres = new Map();
  const offsets = {};
  const shares = new Map();
  if (!root) {
    for (const member of tree.members) offsets[member.ref] = [0, 0, 0];
    return { centres: centres, offsets: offsets, shares: shares };
  }
  const start = scene.epoch_unix == null ? tree.nodes[0].epoch_unix : scene.epoch_unix;
  const unix = start + years * SECONDS_PER_YEAR;
  centres.set(root.node_no, [0, 0, 0]);
  const pending = [root];
  while (pending.length) {
    const node = pending.pop();
    const slots = children.get(node.node_no) || new Map();
    if (node.kind === "group") {
      for (const child of slots.values()) {
        if (child.node) {
          centres.set(child.node.node_no, centres.get(node.node_no));
          pending.push(child.node);
        } else {
          offsets[child.member.ref] = [0, 0, 0];
        }
      }
      continue;
    }
    const low = slots.get(0);
    const high = slots.get(1);
    if (!low || !high) throw new Error(`node ${node.node_no} needs a child in slot 0 and in slot 1`);
    const massLow = (low.member || low.node).mass_solar;
    const massHigh = (high.member || high.node).mass_solar;
    const relative = nodeRelativeKm(node, unix);
    const total = massLow + massHigh;
    shares.set(node.node_no, { low: -massHigh / total, high: massLow / total, relative: relative });
    for (const [child, share] of [[low, -massHigh / total], [high, massLow / total]]) {
      const offset = relative.map((c) => c * share);
      if (child.node) {
        centres.set(child.node.node_no, add(centres.get(node.node_no), offset));
        pending.push(child.node);
      } else {
        offsets[child.member.ref] = offset;
      }
    }
  }
  const missing = tree.members.filter((member) => !offsets[member.ref]).map((member) => member.member_no);
  if (missing.length) throw new Error(`members ${missing} hang from no node of the tree`);
  return { centres: centres, offsets: offsets, shares: shares };
}

// The ref a tree node's centre of mass goes by in relativeAt()'s answer (it is a point, not a body).
export const nodeRef = (nodeNo) => `node:${nodeNo}`;

// Every star's position `years` after the scene's epoch, km, in the scene's frame: {ref: [x, y, z]}, the first star of a
// wide pair at the origin (body_positions._star_positions_km).
function starPositionsKm(scene, years) {
  const { centres, offsets } = treeAt(scene, years);
  const parents = new Map(scene.tree.members.map((member) => [member.ref, member.parent_node_no]));
  const out = {};
  for (const ref of Object.keys(offsets)) {
    const parent = parents.get(ref);
    out[ref] = centres.has(parent) ? add(centres.get(parent), offsets[ref]) : offsets[ref];
  }
  const shift = out[scene.tree.origin] || [0, 0, 0];
  for (const ref of Object.keys(out)) out[ref] = out[ref].map((c, axis) => c - shift[axis]);
  return out;
}

// Where every body is, `years` after the scene's epoch, as a position
// relative to what it goes round: {ref: {around, rel: [x, y, z] km}}.
// `around` is "barycenter" or another body's ref. The stars hang in the tree
// of orbits (GEN.311): each is placed relative to the centre of mass of the
// node it hangs from, and each node's centre relative to its parent node's,
// the root's at the barycenter. A node's centre is a point, not a body: its
// entry is `virtual` and its ref is `nodeRef(node_no)`. A group node's centre
// is its parent's. Belts have no single point and are left out.
export function relativeAt(scene, years) {
  const out = {};
  if (scene.stars.length) {
    const { centres, offsets } = treeAt(scene, years);
    const parentOf = new Map(scene.tree.nodes.map((node) => [node.node_no, node.parent_node_no]));
    for (const [nodeNo, centre] of centres) {
      const parent = parentOf.get(nodeNo);
      const parentCentre = parent == null ? [0, 0, 0] : centres.get(parent);
      out[nodeRef(nodeNo)] = {
        around: parent == null ? "barycenter" : nodeRef(parent), rel: centre.map((c, axis) => c - parentCentre[axis]), virtual: true,
      };
    }
    const memberNode = new Map(scene.tree.members.map((member) => [member.ref, member.parent_node_no]));
    for (const star of scene.stars) {
      const node = memberNode.get(star.ref);
      out[star.ref] = { around: centres.has(node) ? nodeRef(node) : "barycenter", rel: offsets[star.ref] };
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
// in the scene's frame: the stars where the tree puts them, the first star of
// a wide pair at the origin, every other body at the sum of its offsets up
// what it goes round. Tree node centres are not bodies and are left out.
export function positionsAt(scene, years) {
  const relative = relativeAt(scene, years);
  const out = scene.stars.length ? starPositionsKm(scene, years) : {};
  const resolve = (ref) => {
    if (out[ref]) return out[ref];
    const { around, rel } = relative[ref];
    out[ref] = around === "barycenter" ? rel : add(resolve(around), rel);
    return out[ref];
  };
  for (const ref of Object.keys(relative)) {
    if (!relative[ref].virtual) resolve(ref);
  }
  return out;
}

// The paths the nodes of the scene's tree trace in their own frames (MAP.213), for a view to draw: one per child of each
// node, {ref, around, points}, the child's ref (a star's, or `nodeRef` of a sub-pair's centre), the node centre it goes round
// (`nodeRef`) and its path about that centre in km, `samples` points from periapsis, the first repeated at the end. A node's
// child on the slot-0 side sits opposite the slot-1 child, each at the other's share of the separation. A group node has no
// orbit and no paths.
export function nodePaths(scene, samples = 128) {
  const tree = scene.tree;
  const { children } = treeLayout(tree);
  const out = [];
  for (const node of tree.nodes) {
    if (node.kind === "group") continue;
    const slots = children.get(node.node_no);
    const low = slots.get(0);
    const high = slots.get(1);
    const massLow = (low.member || low.node).mass_solar;
    const massHigh = (high.member || high.node).mass_solar;
    const total = massLow + massHigh;
    const path = nodeOrbitPath(node, samples);
    for (const [child, share] of [[low, -massHigh / total], [high, massLow / total]]) {
      out.push({
        ref: child.member ? child.member.ref : nodeRef(child.node.node_no), around: nodeRef(node.node_no),
        points: path.map((point) => point.map((c) => c * share)),
      });
    }
  }
  return out;
}

// The points of one tree node's relative orbit, km: the slot-1 child about the slot-0 child (nodeRelativeKm over a whole
// turn), `samples` points from periapsis, the first repeated at the end.
function nodeOrbitPath(node, samples) {
  const e = node.eccentricity || 0;
  const semiLatus = (node.separation_km / AU_KM) * (1 - e * e);
  const out = [];
  for (let n = 0; n <= samples; n += 1) {
    const nu = (2 * Math.PI * n) / samples;
    const p = orbitalPosition(semiLatus / (1 + e * Math.cos(nu)), node.inclination_deg, node.ascending_node_deg, deg(rad(node.arg_periapsis_deg) + nu));
    out.push([p[0] * AU_KM, p[1] * AU_KM, p[2] * AU_KM]);
  }
  return out;
}

// The points of one orbit's path, relative to what it goes round, in km: a
// planet's or moon's ellipse (a circle when its eccentricity is 0), a
// comet's ellipse or the near part of its parabola (out to `maxAu`).
// `samples` points, the first repeated at the end of a closed path.
export function orbitPath(orbit, samples = 128, maxAu = COMET_PATH_LIMIT_AU) {
  const out = [];
  if (!orbit.kepler) {
    const e = orbit.eccentricity || 0;
    const semiLatus = (orbit.distance_km / AU_KM) * (1 - e * e);
    const periapsis = rad(orbit.arg_periapsis_deg || 0);
    for (let n = 0; n <= samples; n += 1) {
      const nu = (2 * Math.PI * n) / samples;
      const p = orbitalPosition(semiLatus / (1 + e * Math.cos(nu)), orbit.inclination_deg, orbit.ascending_node_deg, deg(periapsis + nu));
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
  fields.push(["Eccentricity", formatEccentricity(k.eccentricity, orbit.type)]);
  return fields;
}

// An eccentricity has no unit, so it says what the number is and where it sits on its scale (UX.210).
export function formatEccentricity(eccentricity, type) {
  return eccentricity.toFixed(type === "elliptical" ? 3 : 4) + " (a ratio: 0 is a circle, 1 is a parabola)";
}
