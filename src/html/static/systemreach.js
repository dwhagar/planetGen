// html/static/systemreach.js
//
// Which star's system the Galaxy Map has up while the camera is in an opened sector, and how much of it shows
// (MAP.152, docs/design/fly-through-view-distance.md sections 4.4 and 4.7): a system's disc becomes drawable when it
// subtends about `LOAD_PX` pixels, fades in to `FULL_PX`, and is dropped below `UNLOAD_PX`. The gap between the
// two ends is the hysteresis: a camera hovering near the size does not load and drop the system again and
// again. A system's real size is not known until its scene is loaded, so until then it is `NOMINAL_RADIUS_PC`
// (50 AU: about 3,000 AU, 0.014 pc, from a 1000 px map at 50 degrees of view).
//
// Pure functions of numbers (no DOM, no three.js), so tests/js/systemreach.test.mjs holds them to the spec.
// `galaxystageview.js` feeds them the camera and the stars of the open sector.

export const LOAD_PX = 20; // the system's disc is this many pixels across: it starts to fade in
export const FULL_PX = 40; // ... and is whole
export const UNLOAD_PX = 10; // ... and below this it is dropped
export const NOMINAL_RADIUS_PC = 50 / 206264.806; // 50 AU
// A system takes over from the one up only when its disc is clearly the larger.
export const TAKE_OVER_SHARE = 0.85;
// How long (ms) a system must stay in range before its scene is asked for.
export const DWELL_MS = 250;

function smoothstep(lo, hi, x) {
  const t = Math.min(1, Math.max(0, (x - lo) / (hi - lo)));
  return t * t * (3 - 2 * t);
}

// How many pixels across a system of `radiusPc` looks from `distance` pc, on a map `heightPx` tall with a vertical
// field of view of `fovDeg`.
export function discPixels(distance, radiusPc, heightPx, fovDeg) {
  if (!(distance > 0) || !(radiusPc > 0)) return distance > 0 ? 0 : Infinity;
  return (radiusPc * heightPx) / (distance * Math.tan((fovDeg * Math.PI) / 360));
}

// The camera distance (pc) at which a system of `radiusPc` is `px` pixels across.
export function distanceForPixels(px, radiusPc, heightPx, fovDeg) {
  return (radiusPc * heightPx) / (px * Math.tan((fovDeg * Math.PI) / 360));
}

// How much of a system shows when its disc is `px` pixels across: 0 at LOAD_PX, 1 at FULL_PX.
export function systemShare(px) {
  return smoothstep(LOAD_PX, FULL_PX, px);
}

// The system whose scene the camera at `eye` ([x, y, z] pc) wants up, or null.
//   candidates  [{key, center: [x, y, z], radiusPc (optional: the real one once its scene is known)}];
//   view        {heightPx, fovDeg};
//   current     the key of the system whose scene is up now, or null.
// Returns {key, px, share}: `share` is 0 for a scene kept up between UNLOAD_PX and LOAD_PX.
export function chooseSystem(eye, candidates, view, current) {
  let nearest = null;
  let kept = null;
  candidates.forEach(function (candidate) {
    const distance = Math.hypot(eye[0] - candidate.center[0], eye[1] - candidate.center[1], eye[2] - candidate.center[2]);
    const px = discPixels(distance, candidate.radiusPc || NOMINAL_RADIUS_PC, view.heightPx, view.fovDeg);
    const found = { key: candidate.key, px: px, share: systemShare(px) };
    if (candidate.key === current) kept = found;
    if (px >= LOAD_PX && (!nearest || px > nearest.px)) nearest = found;
  });
  if (kept && kept.px < UNLOAD_PX) kept = null;
  if (!kept) return nearest;
  if (nearest && nearest.key !== kept.key && nearest.px * TAKE_OVER_SHARE > kept.px) return nearest;
  return kept;
}
