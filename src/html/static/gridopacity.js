// html/static/gridopacity.js
//
// MAP.188: the translucent blocks of a drill-down stage add up along the
// line of sight. A ray through n cells meets n faces, so the picture's
// opacity is 1 - (1 - a)^n, which at deep zoom (a free camera looking through
// many cells) makes the cube grid nearly opaque. Boss (MAP.128): "never under
// any frame be more opaque than 50%". So the faces' (and lines') alpha is
// scaled by a share that keeps the worst ray at or under GRID_ACCUMULATION_CAP,
// worked out from the stage's cells and the way the camera looks. Free of
// the page, so the shader's uniform and the tests share one definition.

// The most the blocks' faces and lines may add up to along any ray.
export const GRID_ACCUMULATION_CAP = 0.5;

function median(values) {
  const sorted = values.slice().sort(function (a, b) { return a - b; });
  return sorted[Math.floor(sorted.length / 2)];
}

// The stage's cells (each {r0, r1, c0, c1, t0, t1, opacity}, spherical about the
// galaxy's center: radii, colatitudes and longitudes) as the numbers the
// estimate needs: the extent of the whole stage along each axis, the typical
// cell edge (the cube root of its volume, the cells being near cubes: a shell
// is as deep as a sector is wide) and the most opaque cell. null for no cells.
export function cellStats(cells) {
  if (!cells || !cells.length) return null;
  const edges = [];
  const low = [Infinity, Infinity, Infinity];
  const high = [-Infinity, -Infinity, -Infinity];
  let maxOpacity = 0;
  cells.forEach(function (c) {
    const r = (c.r0 + c.r1) / 2;
    const colat = (c.c0 + c.c1) / 2;
    const t = (c.t0 + c.t1) / 2;
    const plane = r * Math.sin(colat);
    const at = [plane * Math.cos(t), plane * Math.sin(t), r * Math.cos(colat)];
    for (let k = 0; k < 3; k++) {
      low[k] = Math.min(low[k], at[k]);
      high[k] = Math.max(high[k], at[k]);
    }
    const volume = ((c.r1 * c.r1 * c.r1 - c.r0 * c.r0 * c.r0) / 3) * (Math.cos(c.c0) - Math.cos(c.c1)) * (c.t1 - c.t0);
    edges.push(Math.cbrt(Math.max(volume, 1e-27)));
    maxOpacity = Math.max(maxOpacity, c.opacity || 0);
  });
  const edge = median(edges);
  return {
    edge: edge,
    extent: [high[0] - low[0] + edge, high[1] - low[1] + edge, high[2] - low[2] + edge],
    maxOpacity: maxOpacity,
  };
}

// About how many cells the longest ray through the stage crosses, looking
// along `direction` (a unit vector, any sign): its path inside the stage's
// box times the cells per unit length along it (a ray crosses a face of a
// cubic cell for every |component| of its direction per edge). At least 1.
export function crossings(stats, direction) {
  if (!stats) return 1;
  let path = Infinity;
  let rate = 0;
  for (let k = 0; k < 3; k++) {
    const d = Math.abs(direction[k]);
    if (d > 1e-9) path = Math.min(path, stats.extent[k] / d);
    rate += d / stats.edge;
  }
  if (!isFinite(path)) return 1;
  return Math.max(1, path * rate);
}

// What n faces of opacity `a` add up to.
export function accumulated(n, a) {
  return 1 - Math.pow(1 - Math.min(1, Math.max(0, a)), n);
}

// The share (0..1) to multiply every face's and line's alpha by, so n of the
// most opaque (`maxOpacity`) add up to at most GRID_ACCUMULATION_CAP.
export function gridScale(n, maxOpacity) {
  if (!(maxOpacity > 0) || !(n > 1)) return 1;
  const each = 1 - Math.pow(1 - GRID_ACCUMULATION_CAP, 1 / n);
  return Math.min(1, each / maxOpacity);
}
