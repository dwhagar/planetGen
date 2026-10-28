// html/static/flytofit.js
//
// Where the camera goes to look at one thing so that it fills the window (MAP.319, docs/design/interface/fly-through-view-distance.md):
// the thing is a bounding sphere (a star system, a planet and its moons, a nebula, a region), and `fitDistance` is how
// far from its middle the camera sits for that sphere to fill the smaller side of the view. The distance is the sphere's
// radius over the sine of half the smaller field of view (the sphere touches the edge of the view along that side),
// plus a margin of `FIT_MARGIN`, never less than for a sphere of `options.minRadius`. The viewing direction is the
// map's to keep; the zoom policy of the map clamps the result to what the view allows (MAP.152).
//
// Pure functions of numbers (no DOM, no three.js), so tests/js/flytofit.test.mjs holds them to the spec. The Galaxy
// Map (galaxystageview.js) and the 3D system view (systemview3d.js) feed them the sphere.

// The share added to the distance so the target does not touch the window's edge.
export const FIT_MARGIN = 0.08;
// The smallest sphere (pc) the Galaxy Map frames a target as, whatever it measures: a nebula or region smaller than
// this is looked at as a view of this radius.
export const MIN_VIEW_RADIUS_PC = 2;

// Half of the smaller of the two fields of view, in radians: `vfovDeg` is the camera's vertical field of view and
// `aspect` the window's width over its height (the horizontal field follows from the two).
export function halfFieldRad(vfovDeg, aspect) {
  const vertical = Math.tan((vfovDeg * Math.PI) / 360);
  const horizontal = vertical * (aspect > 0 ? aspect : 1);
  return Math.atan(Math.min(vertical, horizontal));
}

// The distance (in the units of `radius`) from the sphere's middle to the camera for the sphere to fill the view.
//   options.margin     the share added (default FIT_MARGIN)
//   options.minRadius  a sphere smaller than this is taken as this size (default 0: the map's zoom limit clamps)
export function fitDistance(radius, vfovDeg, aspect, options) {
  const margin = options && options.margin !== undefined ? options.margin : FIT_MARGIN;
  const least = options && options.minRadius > 0 ? options.minRadius : 0;
  const size = Math.max(Number.isFinite(radius) ? radius : 0, least);
  return (size / Math.sin(halfFieldRad(vfovDeg, aspect))) * (1 + margin);
}

// How far the farthest part of `spheres` reaches from `centre`: spheres are {x, y, z, r} (r the sphere's own radius,
// 0 for a point). Spheres without a position are skipped.
export function reachFrom(centre, spheres) {
  let far = 0;
  spheres.forEach(function (s) {
    if (!Number.isFinite(s.x) || !Number.isFinite(s.y) || !Number.isFinite(s.z)) return;
    far = Math.max(far, Math.hypot(s.x - centre[0], s.y - centre[1], s.z - centre[2]) + (s.r > 0 ? s.r : 0));
  });
  return far;
}

// The farthest any vertex of a flat position list ([x, y, z, x, y, z, ...], an orbit's line) lies from the origin
// of its own frame: how far the orbit reaches from what it goes round, its apoapsis for an ellipse.
export function reachOfPositions(flat) {
  let far = 0;
  for (let n = 0; n + 2 < flat.length; n += 3) {
    far = Math.max(far, Math.hypot(flat[n], flat[n + 1], flat[n + 2]));
  }
  return far;
}
