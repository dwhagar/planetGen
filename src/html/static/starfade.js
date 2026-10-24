// MAP.153 and MAP.152: stars fade in and out as the camera nears or leaves
// the tile that lists them. The map fetches the distance cut of tiles
// (static/galaxytilecut.js: near ones small, far ones big) and hands over
// from a tile to its eight halves over a band of distance, so nothing pops.
//
// Every tile T has a "leaving" share from the camera's distance to its box,
//
//     u = CUT_EDGES * distance / edge,   L(T) = smoothstep((u - 1/B) / (B - 1/B))
//
// (B = FADE_BAND): 1 while the camera is at least B edges' worth away, 0 once it
// is within 1/B of that, so a tile gives way to its halves as the camera comes
// in. A tile is shown for what its ancestors have left of it,
//
//     share(T) = L(T) * product of (1 - L(A)) over T's ancestors A,
//
// so across a hand-over the shares of a tile and of its halves add up to 1: a star
// listed in both stays whole, one listed in one comes in or goes out smoothly.
// The finest tile of the cut has nothing finer to give way to (L = 1 there). The
// finest "detail" tiles round the target fade by distance alone, whole within
// DETAIL_FULL_PC and gone at DETAIL_GONE_PC. docs/design/interface/fly-through-view-distance.md
// sections 4.5 and 4.7.
//
// A star's `fade` is its tile: `{ parts: [level, ix, iy, iz], detail }`, or
// `null` for a star in no tile (carried over while its tile loads, or a
// landmark), which is whole.

// Equal to planetgen.galaxy.viewport.DISTANCE_CUT_EDGES and DISTANCE_CUT_FADE_BAND
// (a test holds them equal).
export const CUT_EDGES = 3;
export const FADE_BAND = 1.15;
export const DETAIL_FULL_PC = 60;
export const DETAIL_GONE_PC = 120;

export function smoothstep(x) {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
}

// How far `point` is from the box of the tile `parts`.
export function boxDistance(parts, point, rootEdge) {
  const edge = rootEdge / 2 ** parts[0];
  let sum = 0;
  for (let axis = 0; axis < 3; axis++) {
    const lo = -rootEdge / 2 + parts[axis + 1] * edge;
    const gap = Math.max(lo - point[axis], point[axis] - (lo + edge), 0);
    sum += gap * gap;
  }
  return Math.sqrt(sum);
}

// The leaving share of the tile `parts` with the camera `distance` from its box.
export function leavingShare(distance, edge) {
  const u = CUT_EDGES * distance / edge;
  return smoothstep((u - 1 / FADE_BAND) / (FADE_BAND - 1 / FADE_BAND));
}

// A tile's share (0..1) with the camera at `camera`. `rootEdge` is the root
// tile's edge and `finestLevel` the finest level the cut itself reaches.
export function tileShare(parts, camera, rootEdge, finestLevel) {
  let share = 1;
  let [level, ix, iy, iz] = parts;
  while (level > 0 && share > 0) {
    level -= 1;
    ix >>= 1;
    iy >>= 1;
    iz >>= 1;
    const edge = rootEdge / 2 ** level;
    share *= 1 - leavingShare(boxDistance([level, ix, iy, iz], camera, rootEdge), edge);
  }
  if (share > 0 && parts[0] < finestLevel) {
    share *= leavingShare(boxDistance(parts, camera, rootEdge), rootEdge / 2 ** parts[0]);
  }
  return share;
}

// A detail tile's share: whole while the camera is within DETAIL_FULL_PC.
export function detailShare(distance) {
  return smoothstep((DETAIL_GONE_PC - distance) / (DETAIL_GONE_PC - DETAIL_FULL_PC));
}

// A star's opacity from the tile that lists it (`fade`, see above), with the
// camera at `camera` = [x, y, z].
export function starShare(fade, camera, rootEdge, finestLevel) {
  if (!fade) {
    return 1;
  }
  return fade.detail ? detailShare(boxDistance(fade.parts, camera, rootEdge))
    : tileShare(fade.parts, camera, rootEdge, finestLevel);
}

// A star that is in no tile (carried over, or a landmark) is whole.
export const POINT_FADE = null;

// MAP.155: a cloud sprite fades in between CLOUD_RAMP_PX[0] and
// CLOUD_RAMP_PX[1] pixels of radius instead of switching on at 2 px.
export const CLOUD_RAMP_PX = [1, 4];

export function cloudSizeShare(radiusPx) {
  return smoothstep((radiusPx - CLOUD_RAMP_PX[0]) / (CLOUD_RAMP_PX[1] - CLOUD_RAMP_PX[0]));
}

// The sum of `starShare` over `fades` for a camera at `camera`.
export function opacitySum(fades, camera, rootEdge, finestLevel) {
  let total = 0;
  for (const fade of fades) {
    total += starShare(fade, camera, rootEdge, finestLevel);
  }
  return total;
}
