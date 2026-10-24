// html/static/sectorreach.js
//
// Which sector's own scene the Galaxy Map has up while the camera is still in a window of
// sectors, and how much of it shows (MAP.152, docs/design/interface/fly-through-view-distance.md
// section 4.7): a sector's scene (every star, phenomenon and body of it, galaxysector.js)
// fades in as the camera comes within `FADE_IN_EDGES` sector edges of the sector's centre,
// is whole by `FULL_EDGES`, and is dropped past `UNLOAD_EDGES`. The gap between the two
// ends is the hysteresis: a camera hovering near the edge of the range does not load and
// drop the scene again and again, and a scene already loaded stays put while the camera
// wanders a little further than it was fetched at.
//
// Pure functions of numbers (no DOM, no three.js), so tests/js/sectorreach.test.mjs holds
// them to the spec. `galaxystageview.js` feeds them the camera and the sectors in view.

// The camera's distance, in sector edges, from a sector's centre ...
export const FADE_IN_EDGES = 1.5; // ... inside which its scene starts to fade in,
export const FULL_EDGES = 1; // ... inside which it is whole,
export const UNLOAD_EDGES = 3; // ... and past which it is dropped.

// One scene is up at a time, so a sector takes over from the one up only when it is clearly the
// nearer: its distance under this share of the other's. Two sectors the camera is about
// equally far from do not trade places on every frame.
export const TAKE_OVER_SHARE = 0.85;

// How long (ms) a sector must stay in range before its scene is asked for, so a camera
// passing through does not fetch the sectors it flies by.
export const DWELL_MS = 250;

function smoothstep(lo, hi, x) {
  const t = Math.min(1, Math.max(0, (x - lo) / (hi - lo)));
  return t * t * (3 - 2 * t);
}

// How much of a sector's scene shows with the camera `distance` pc from its centre, `edge` pc
// across: 0 from FADE_IN_EDGES out, 1 from FULL_EDGES in, smooth between.
export function sectorShare(distance, edge) {
  if (!(edge > 0)) return 0;
  return 1 - smoothstep(FULL_EDGES, FADE_IN_EDGES, distance / edge);
}

function distanceBetween(a, b) {
  return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

// The sector whose scene the camera at `eye` ([x, y, z] pc) wants up, or null.
//   candidates  [{key, center: [x, y, z]}]: the sectors in view that have a scene to show;
//   edge        a sector's edge, pc;
//   current     the key of the sector whose scene is up now, or null.
// Returns {key, distance, share}: `share` is how much of it shows (0 for a scene kept up
// between FADE_IN_EDGES and UNLOAD_EDGES, loaded but not yet drawn).
export function chooseSector(eye, candidates, edge, current) {
  let nearest = null;
  let kept = null;
  candidates.forEach(function (candidate) {
    const distance = distanceBetween(eye, candidate.center);
    const found = { key: candidate.key, distance: distance, share: sectorShare(distance, edge) };
    if (candidate.key === current) kept = found;
    if (distance < FADE_IN_EDGES * edge && (!nearest || distance < nearest.distance)) nearest = found;
  });
  if (kept && kept.distance > UNLOAD_EDGES * edge) kept = null;
  if (!kept) return nearest;
  if (nearest && nearest.key !== kept.key && nearest.distance < TAKE_OVER_SHARE * kept.distance) return nearest;
  return kept;
}
