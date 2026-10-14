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

// The stage's cells (each {r0, r1, t0, t1, z0, z1, opacity}, polar about the
// galaxy's axis) as the numbers the estimate needs: the extent of the whole
// stage and the typical cell edge, horizontally and vertically (world
// units), and the most opaque cell. null for no cells.
export function cellStats(cells) {
  if (!cells || !cells.length) return null;
  const across = [];
  const high = [];
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity, minZ = Infinity, maxZ = -Infinity;
  let maxOpacity = 0;
  cells.forEach(function (c) {
    const r = (c.r0 + c.r1) / 2;
    const t = (c.t0 + c.t1) / 2;
    const x = r * Math.cos(t);
    const y = r * Math.sin(t);
    const z = (c.z0 + c.z1) / 2;
    minX = Math.min(minX, x); maxX = Math.max(maxX, x);
    minY = Math.min(minY, y); maxY = Math.max(maxY, y);
    minZ = Math.min(minZ, z); maxZ = Math.max(maxZ, z);
    across.push(Math.sqrt(Math.max(c.r1 - c.r0, 1e-9) * Math.max(r * (c.t1 - c.t0), 1e-9)));
    high.push(Math.max(c.z1 - c.z0, 1e-9));
    maxOpacity = Math.max(maxOpacity, c.opacity || 0);
  });
  const edgeH = median(across);
  const edgeV = median(high);
  return {
    edgeH: edgeH, edgeV: edgeV,
    extentH: Math.max(maxX - minX, maxY - minY) + edgeH,
    extentV: maxZ - minZ + edgeV,
    maxOpacity: maxOpacity,
  };
}

// About how many cells the longest ray through the stage crosses, looking
// along `direction` (a unit vector, any sign): its path inside the stage's
// box times the cells per unit length along it. At least 1.
export function crossings(stats, direction) {
  if (!stats) return 1;
  const h = Math.hypot(direction[0], direction[1]);
  const v = Math.abs(direction[2]);
  let path = Infinity;
  if (h > 1e-9) path = Math.min(path, stats.extentH / h);
  if (v > 1e-9) path = Math.min(path, stats.extentV / v);
  if (!isFinite(path)) return 1;
  return Math.max(1, path * (h / stats.edgeH + v / stats.edgeV));
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
