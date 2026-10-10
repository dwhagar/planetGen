// html/static/nearfield.js
//
// MAP.149: the near field of the fly-through map. What stands close to the
// camera dissolves (a depth fade), a soft see-through tube thins what stands
// between the camera and the thing looked at (the focus), and the container
// the camera is inside is drawn from the inside. docs/design/
// fly-through-view-distance.md (4.2 and 4.6) has the rule; this file is its
// arithmetic, free of the page, so the shaders (which repeat `nearField` in
// GLSL, built from the same constants below: NEAR_FIELD_GLSL), the picker
// (galaxystageview.js, galaxymap3d.js) and the tests share one definition.
//
// Everything is worked out in VIEW SPACE: the camera at the origin looking
// down -z at the focus, which sits at (0, 0, -D), D the camera-to-focus
// distance (the orbit radius). So a point's depth is -z, and the tube round
// the line from the camera to the focus is just the distance from the z axis.
//
//     depth fade   smoothstep(Z0 D, Z1 D, depth)
//     tube         ahead = smoothstep(0, Rt, D - depth)      (nearer than the focus)
//                  tube  = 1 - (1 - A_MIN) (1 - smoothstep(Rt, 1.6 Rt, rho)) ahead
//     Rt           max(focus radius, 0.12 D)
//
// `ahead` ramps the tube in over a distance Rt in front of the focus, so
// flat slabs lying in the focus plane are untouched and nothing pops as
// the focus plane is crossed. Picking uses the same number at the point it
// hits: only what is visible enough can be picked, so picking agrees with
// drawing.

export const DEPTH_FADE_START = 0.15;
export const DEPTH_FADE_END = 0.5;
export const TUBE_MIN_ALPHA = 0.08;
export const TUBE_RADIUS_FLOOR = 0.12;
export const TUBE_EDGE = 1.6;
// A near-field share below this can't be picked (it is drawn too faint to
// have been seen), so the pick falls through to what is behind.
export const PICK_MIN_VISIBLE = 0.35;
// The faces of the container the camera is in, seen from inside: at this
// share of the block's own opacity (its edges at full).
export const INSIDE_FACE_ALPHA = 0.25;

// Prominence of a region (4.2): the one looked at and the container the
// camera is in at (nearly) full strength, the rest faint, thinning with
// their distance from the focus.
export const CONTAINER_PROMINENCE = 0.8;
export const CONTEXT_MIN = 0.08;
export const CONTEXT_MAX = 0.3;
// How much shallower (magnitudes) a faint region's star limit is than the
// focus's (MAP.148 uses limitingMagnitude).
export const CONTEXT_MAGNITUDE_DROP = 2.5;

export function smoothstep(edge0, edge1, x) {
  const t = Math.min(1, Math.max(0, (x - edge0) / (edge1 - edge0)));
  return t * t * (3 - 2 * t);
}

// 0 within DEPTH_FADE_START of the focus distance D from the camera, 1
// from DEPTH_FADE_END of it on.
export function depthFade(depth, D) {
  if (!(D > 0)) return 1;
  return smoothstep(DEPTH_FADE_START * D, DEPTH_FADE_END * D, depth);
}

// The tube's radius for focus distance D and a focus `radius` wide (0 when
// unknown).
export function tubeRadius(D, radius) {
  return Math.max(radius || 0, TUBE_RADIUS_FLOOR * D);
}

// 1 outside the tube or at and behind the focus, TUBE_MIN_ALPHA inside it
// in front of the focus; (x, y, depth) a point in view space.
export function tubeFade(x, y, depth, D, radius) {
  if (!(D > 0)) return 1;
  const rt = tubeRadius(D, radius);
  const ahead = smoothstep(0, rt, D - depth);
  const inside = 1 - smoothstep(rt, TUBE_EDGE * rt, Math.hypot(x, y));
  return 1 - (1 - TUBE_MIN_ALPHA) * inside * ahead;
}

// The near field's share of a point at (x, y, z) in view space (z < 0 in
// front of the camera).
export function nearField(x, y, z, D, radius) {
  return depthFade(-z, D) * tubeFade(x, y, -z, D, radius);
}

// A world point in view space, from a camera's matrixWorldInverse
// elements (column-major, as three.js keeps them).
export function viewSpace(e, x, y, z) {
  return [
    e[0] * x + e[4] * y + e[8] * z + e[12],
    e[1] * x + e[5] * y + e[9] * z + e[13],
    e[2] * x + e[6] * y + e[10] * z + e[14],
  ];
}

// The near field's share at a world point, for the picker.
export function nearFieldAtWorld(e, x, y, z, D, radius) {
  const v = viewSpace(e, x, y, z);
  return nearField(v[0], v[1], v[2], D, radius);
}

// Whether a near-field share is enough to be picked.
export function visibleEnough(share) {
  return share >= PICK_MIN_VISIBLE;
}

// The opacity of a faint region `distance` from the focus, `lambda` being
// about one region width: a0 exp(-d / lambda), kept between CONTEXT_MIN
// (so its ring stays visible) and CONTEXT_MAX.
export function contextOpacity(distance, lambda) {
  const share = lambda > 0 ? Math.exp(-distance / lambda) : 0;
  return Math.min(CONTEXT_MAX, Math.max(CONTEXT_MIN, share));
}

// A region's prominence p: 1 for the one looked at, CONTAINER_PROMINENCE
// for the container the camera is in, else contextOpacity.
export function prominence(kind, distance, lambda) {
  if (kind === "focus") return 1;
  if (kind === "container") return CONTAINER_PROMINENCE;
  return contextOpacity(distance, lambda);
}

// The limiting magnitude of a region of prominence p, given the focus's.
export function limitingMagnitude(focusLimit, p) {
  return focusLimit - CONTEXT_MAGNITUDE_DROP * (1 - p);
}

// Whether a camera at (x, y, z) is inside a block with bounds {r0, r1, t0,
// t1, z0, z1} (polar around the galaxy's axis, as the stage API gives them).
export function blockContains(b, x, y, z) {
  if (z < b.z0 || z >= b.z1) return false;
  const r = Math.hypot(x, y);
  if (r < b.r0 || r >= b.r1) return false;
  if (r === 0) return true;
  const t = ((Math.atan2(y, x) % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
  return t >= b.t0 && t < b.t1;
}

// A block's center in world space ([x, y, z]), the middle of its bounds,
// as galaxyprisms.cellCoordinates works it (the block shader compares it
// with its vertices' `prismCenter` to tell which block is the container).
export function blockCenter(b) {
  const r = (b.r0 + b.r1) / 2;
  const t = (b.t0 + b.t1) / 2;
  return [r * Math.cos(t), r * Math.sin(t), (b.z0 + b.z1) / 2];
}

// The same near-field arithmetic in GLSL: `float nearField(vec3 view, float D,
// float radius)`, the shaders' twin of nearField above. Built from the
// constants, so the two can't drift.
export const NEAR_FIELD_GLSL = [
  "float nfStep(float a, float b, float x) {",
  "  float t = clamp((x - a) / (b - a), 0.0, 1.0);",
  "  return t * t * (3.0 - 2.0 * t);",
  "}",
  "float nearField(vec3 view, float D, float radius) {",
  "  if (D <= 0.0) return 1.0;",
  "  float depth = -view.z;",
  "  float rt = max(radius, " + TUBE_RADIUS_FLOOR.toFixed(4) + " * D);",
  "  float ahead = nfStep(0.0, rt, D - depth);",
  "  float inside = 1.0 - nfStep(rt, " + TUBE_EDGE.toFixed(4) + " * rt, length(view.xy));",
  "  float tube = 1.0 - (1.0 - " + TUBE_MIN_ALPHA.toFixed(4) + ") * inside * ahead;",
  "  return nfStep(" + DEPTH_FADE_START.toFixed(4) + " * D, " + DEPTH_FADE_END.toFixed(4) + " * D, depth) * tube;",
  "}",
].join("\n");
