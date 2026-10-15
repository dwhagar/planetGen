// html/static/galaxytilecut.js
//
// MAP.152: the tiles the Galaxy Map fetches for a camera, chosen by
// distance. The same arithmetic as planetgen.galaxy.viewport's
// `distance_cut_tiles` and `fetch_cut_tiles` (tests/test_viewport.py and
// tests/js/galaxytilecut.test.mjs hold the two equal): start from the root
// cube and split a tile into its eight halves while its edge is more than
// `k` times the distance from the camera to its box, so near space is fine
// and far space coarse, and keep the tiles inside the cone the camera looks
// along. The fetch set is the cuts at `k * FADE_BAND` and `k / FADE_BAND`
// together, the tiles static/starfade.js fades in and out between.

// Equal to planetgen.galaxy.viewport.DISTANCE_CUT_FADE_BAND and static/starfade.js's
// (a test holds them equal); not imported, since a module may not import a sibling
// without its version query.
export const FADE_BAND = 1.15;
export const DISTANCE_CUT_EDGES = 2;
export const DISTANCE_CUT_HALF_ANGLE_DEG = 70;
export const DISTANCE_CUT_MAX_TILES = 512;
const SQRT3_HALF = Math.sqrt(3) / 2;

// The edge of a level-`level` tile.
export function tileEdgeOf(level, rootEdge) {
  return rootEdge / Math.pow(2, level);
}

// A tile key's box, `[lo, edge]`.
export function tileBoxOf(key, rootEdge) {
  const parts = key.split("/").map(Number);
  const edge = tileEdgeOf(parts[0], rootEdge);
  const origin = -rootEdge / 2;
  return { lo: [origin + parts[1] * edge, origin + parts[2] * edge, origin + parts[3] * edge], edge };
}

// The distance cut. `config`: { rootEdge, finestLevel (the finest level the
// cut itself reaches), bounds: [[x, y, z], [x, y, z]] or null, k,
// halfAngleDeg, maxTiles }. Returns tile keys, nearest first.
export function distanceCutTiles(camera, look, config) {
  const rootEdge = config.rootEdge;
  const k = config.k || DISTANCE_CUT_EDGES;
  const half = rootEdge / 2;
  const bounds = config.bounds || [[-half, -half, -half], [half, half, half]];
  const halfAngleDeg = config.halfAngleDeg == null ? DISTANCE_CUT_HALF_ANGLE_DEG : config.halfAngleDeg;
  const maxTiles = config.maxTiles || DISTANCE_CUT_MAX_TILES;
  let axis = [look[0] - camera[0], look[1] - camera[1], look[2] - camera[2]];
  const length = Math.sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
  const cone = halfAngleDeg < 180 && length > 0;
  if (cone) {
    axis = axis.map((v) => v / length);
  }
  const halfAngle = (halfAngleDeg * Math.PI) / 180;
  const found = [];
  const pending = [[0, 0, 0, 0]];
  while (pending.length) {
    const [level, ix, iy, iz] = pending.pop();
    const edge = tileEdgeOf(level, rootEdge);
    const lo = [-half + ix * edge, -half + iy * edge, -half + iz * edge];
    let outside = false;
    for (let a = 0; a < 3; a++) {
      if (lo[a] + edge <= bounds[0][a] || lo[a] > bounds[1][a]) outside = true;
    }
    if (outside) continue;
    let distanceSq = 0;
    for (let a = 0; a < 3; a++) {
      const nearest = Math.min(Math.max(camera[a], lo[a]), lo[a] + edge);
      distanceSq += (camera[a] - nearest) * (camera[a] - nearest);
    }
    if (cone && distanceSq > 0) {
      const offset = [0, 1, 2].map((a) => lo[a] + edge / 2 - camera[a]);
      const reach = Math.sqrt(offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]);
      const radius = edge * SQRT3_HALF;
      if (reach > radius) {
        const along = (offset[0] * axis[0] + offset[1] * axis[1] + offset[2] * axis[2]) / reach;
        if (Math.acos(Math.max(-1, Math.min(1, along))) > halfAngle + Math.asin(radius / reach)) continue;
      }
    }
    if (level < config.finestLevel && edge * edge > k * k * distanceSq) {
      for (let dx = 0; dx < 2; dx++) {
        for (let dy = 0; dy < 2; dy++) {
          for (let dz = 0; dz < 2; dz++) {
            pending.push([level + 1, 2 * ix + dx, 2 * iy + dy, 2 * iz + dz]);
          }
        }
      }
    } else {
      found.push([distanceSq, `${level}/${ix}/${iy}/${iz}`]);
    }
  }
  found.sort((p, q) => (p[0] - q[0]) || (p[1] < q[1] ? -1 : p[1] > q[1] ? 1 : 0));
  return found.slice(0, maxTiles).map((entry) => entry[1]);
}

// The tiles to fetch: the cuts at `k * FADE_BAND` and `k / FADE_BAND`, each once, nearest first.
export function fetchCutTiles(camera, look, config) {
  const k = config.k || DISTANCE_CUT_EDGES;
  const keys = new Map();
  [FADE_BAND, 1 / FADE_BAND].forEach((factor) => {
    distanceCutTiles(camera, look, Object.assign({}, config, { k: k * factor })).forEach((key) => keys.set(key, true));
  });
  const distance = (key) => {
    const box = tileBoxOf(key, config.rootEdge);
    let sum = 0;
    for (let a = 0; a < 3; a++) {
      const gap = Math.max(box.lo[a] - camera[a], camera[a] - (box.lo[a] + box.edge), 0);
      sum += gap * gap;
    }
    return sum;
  };
  const ordered = Array.from(keys.keys()).map((key) => [distance(key), key]);
  ordered.sort((p, q) => (p[0] - q[0]) || (p[1] < q[1] ? -1 : p[1] > q[1] ? 1 : 0));
  return ordered.slice(0, config.maxTiles || DISTANCE_CUT_MAX_TILES).map((entry) => entry[1]);
}
