// html/static/camerarelative.js
//
// MAP.240: browser buffers that hold positions as float32 lose precision far from the origin (the spacing of float32 at 262 kpc
// is about 0.016 pc, far above what a close view shows), so a shape built from absolute coordinates jitters and snaps. The
// remedy is the one the star buffer uses, for shapes that are one object each: the vertices are made relative to an origin of
// the object's own, in double precision on the CPU, and the object carries the origin as its position. three.js multiplies the
// model and view matrices in double too, so what the GPU gets is the small offset from the camera.
//
// Pure arithmetic, free of three.js, so a test can hold it. An origin is a plain `[x, y, z]`; a list of points is a flat
// `[x, y, z, x, y, z, ...]` array of doubles.

// How far from its origin a buffer's values may be, in parsecs: float32 keeps about 7 digits, so 10 kpc leaves 1e-3 pc (200 AU),
// a tenth of a pixel at the closest zoom of a map layer (0.01 pc a pixel). The browser test holds every layer in a window to it;
// the one view of the whole galaxy draws its blocks as a single part that spans more, where a pixel is 100 pc or more.
export const SAFE_EXTENT_PC = 10000;

// The middle of the box that holds the points: the origin that keeps a single object's offsets smallest.
export function centreOf(flat) {
  if (!flat.length) return [0, 0, 0];
  const low = [Infinity, Infinity, Infinity];
  const high = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < flat.length; i += 3) {
    for (let axis = 0; axis < 3; axis++) {
      const value = flat[i + axis];
      if (value < low[axis]) low[axis] = value;
      if (value > high[axis]) high[axis] = value;
    }
  }
  return [(low[0] + high[0]) / 2, (low[1] + high[1]) / 2, (low[2] + high[2]) / 2];
}

// The points as float32 offsets from `origin`, each difference worked out in double before it is rounded.
export function offsetsFrom(flat, origin) {
  const out = new Float32Array(flat.length);
  for (let i = 0; i < flat.length; i += 3) {
    out[i] = flat[i] - origin[0];
    out[i + 1] = flat[i + 1] - origin[1];
    out[i + 2] = flat[i + 2] - origin[2];
  }
  return out;
}

// One object's buffer: `{ origin, positions }`, the origin the middle of the points (or the one given).
export function relativeBuffer(flat, origin) {
  const at = origin || centreOf(flat);
  return { origin: at, positions: offsetsFrom(flat, at) };
}

// Points spread over a wide space as a few objects, each with an origin of its own: the points are sorted into cubes `cellPc`
// a side, and each occupied cube becomes `{ origin, positions, indices }` (`indices`: which of the points are in it, in the
// order they were given, so a parallel list such as colours can be split alike). The origin is the cube's middle, so every
// offset is at most half the cube's diagonal.
export function groupByCell(flat, cellPc) {
  const cubes = new Map();
  for (let i = 0, n = 0; i < flat.length; i += 3, n++) {
    const key = [0, 1, 2].map((axis) => Math.floor(flat[i + axis] / cellPc)).join(",");
    let cube = cubes.get(key);
    if (!cube) {
      cube = { cell: key.split(",").map(Number), points: [], indices: [] };
      cubes.set(key, cube);
    }
    cube.points.push(flat[i], flat[i + 1], flat[i + 2]);
    cube.indices.push(n);
  }
  return Array.from(cubes.values()).map((cube) => {
    const origin = cube.cell.map((c) => (c + 0.5) * cellPc);
    return { origin, positions: offsetsFrom(cube.points, origin), indices: cube.indices };
  });
}

// A polyline through `stops` (a list of `[x, y, z]`) as straight pieces no longer than `maxPc`, each `{ origin, positions }` with
// its origin at its own middle and the two end points as float32 offsets from it: a course across the galaxy is then as exact
// as one across a sector. Collinear pieces draw the same line.
export function linePieces(stops, maxPc) {
  const pieces = [];
  for (let leg = 0; leg + 1 < stops.length; leg++) {
    const from = stops[leg];
    const to = stops[leg + 1];
    const length = Math.hypot(to[0] - from[0], to[1] - from[1], to[2] - from[2]);
    const count = Math.max(1, Math.ceil(length / maxPc));
    for (let n = 0; n < count; n++) {
      const a = [0, 1, 2].map((axis) => from[axis] + (to[axis] - from[axis]) * (n / count));
      const b = [0, 1, 2].map((axis) => from[axis] + (to[axis] - from[axis]) * ((n + 1) / count));
      pieces.push(relativeBuffer([a[0], a[1], a[2], b[0], b[1], b[2]]));
    }
  }
  return pieces;
}

// The biggest absolute value in a buffer (what the browser test bounds), 0 for none.
export function extentOf(array) {
  let most = 0;
  for (let i = 0; i < array.length; i++) {
    const value = Math.abs(array[i]);
    if (value > most) most = value;
  }
  return most;
}
