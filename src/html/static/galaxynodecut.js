// html/static/galaxynodecut.js
//
// MAP.261: the sector-block nodes the Galaxy Map needs for a camera, and how much of each to draw. The same
// arithmetic as planetgen.galaxy.nodecut (tests/test_nodecut.py and tests/js/galaxynodecut.test.mjs hold the two
// equal; docs/design-drafts single-space-map-zoom section 2.2): walk the block ladder from its three root blocks and
// keep the blocks inside the view cone, going on to a block's children while its projected size is at least
// FETCH_FROM_PX. A block with children hands its place over to them over HANDOVER_PX of projected size, so a node is
// drawn for what its ancestors have left of it and a node and its children always add up to one.

const VERSION_QUERY = new URL(import.meta.url).search;
const {
  DRILL_LEVELS, drillBlockBounds, drillChildren, drillParent, drillWedgeCount, formatDrillKey,
} = await import(`./galaxyprisms.js${VERSION_QUERY}`);

// Equal to planetgen.galaxy.nodecut's (a test holds them equal).
export const HANDOVER_PX = [70, 140];
export const FETCH_FROM_PX = 56;
export const CONE_HALF_ANGLE_DEG = 50;
export const MAX_NODES = 1200;
export const DEFAULT_FOCAL_PX = 800 / (2 * Math.tan((25 * Math.PI) / 180));
const SQRT3_HALF = Math.sqrt(3) / 2;
const TOP = DRILL_LEVELS[0];

export function roots() {
  const count = drillWedgeCount(TOP, 0, 0);
  const found = [];
  for (let wedge = 0; wedge < count; wedge++) found.push({ m: TOP, shell: 0, band: 0, wedge });
  return found;
}

// A block's bounding box, [lo, hi] in parsecs (planetgen.galaxy.nodecut.block_box): the extremes of
// radius * (sin(colat) cos(lon), sin(colat) sin(lon), cos(colat)) over its bounds, widened a thousandth of a parsec.
export function blockBox(block, edgePc) {
  const b = drillBlockBounds(block, edgePc);
  const cosNorth = Math.cos(b.c0);
  const cosSouth = Math.cos(b.c1);
  const zLow = Math.min(b.r0 * cosSouth, b.r1 * cosSouth);
  const zHigh = Math.max(b.r0 * cosNorth, b.r1 * cosNorth);
  const sinLow = Math.min(Math.sin(b.c0), Math.sin(b.c1));
  const sinHigh = b.c0 <= Math.PI / 2 && Math.PI / 2 <= b.c1 ? 1 : Math.max(Math.sin(b.c0), Math.sin(b.c1));
  const rhoLow = b.r0 * sinLow;
  const rhoHigh = b.r1 * sinHigh;
  const angles = [b.t0, b.t1];
  let turn = Math.ceil(b.t0 / (Math.PI / 2)) * (Math.PI / 2);
  while (turn <= b.t1) {
    angles.push(turn);
    turn += Math.PI / 2;
  }
  const xs = [];
  const ys = [];
  angles.forEach((angle) => [rhoLow, rhoHigh].forEach((rho) => {
    xs.push(rho * Math.cos(angle));
    ys.push(rho * Math.sin(angle));
  }));
  const margin = 1e-3;
  return [
    [Math.min(...xs) - margin, Math.min(...ys) - margin, zLow - margin],
    [Math.max(...xs) + margin, Math.max(...ys) + margin, zHigh + margin],
  ];
}

const geometries = new Map();

// A block's cut geometry, { lo, hi, size } (size the longest edge of the box), kept per sector edge.
export function geometry(block, edgePc) {
  const key = formatDrillKey(block) + "@" + edgePc;
  let found = geometries.get(key);
  if (found === undefined) {
    const box = blockBox(block, edgePc);
    found = { lo: box[0], hi: box[1], size: Math.max(box[1][0] - box[0][0], box[1][1] - box[0][1], box[1][2] - box[0][2]) };
    geometries.set(key, found);
  }
  return found;
}

export function smoothstep(x) {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
}

function distanceSqTo(geo, camera) {
  let sum = 0;
  for (let a = 0; a < 3; a++) {
    const nearest = Math.min(Math.max(camera[a], geo.lo[a]), geo.hi[a]);
    sum += (camera[a] - nearest) * (camera[a] - nearest);
  }
  return sum;
}

// The block's size on screen, pixels; infinite with the camera inside its box.
export function projectedSize(block, camera, edgePc, focalPx) {
  const geo = geometry(block, edgePc);
  const distanceSq = distanceSqTo(geo, camera);
  return distanceSq <= 0 ? Infinity : (geo.size * focalPx) / Math.sqrt(distanceSq);
}

// The share a node gives its children at projected size `sizePx`; 0 for a sector.
export function handover(sizePx, block) {
  if (block.m === 1) return 0;
  return smoothstep((sizePx - HANDOVER_PX[0]) / (HANDOVER_PX[1] - HANDOVER_PX[0]));
}

// A weight calculator for one camera: `weight(block)` is the share of the block's own content to draw (what its
// ancestors have given it less what it gives its children). A block gives children only when the cut expanded it
// (`expanded`, a Set of keys; without one every block with children does), so a node whose children were not
// fetched stays whole. The hand-over of each block is worked out once.
export function weightCalculator(camera, edgePc, focalPx, expanded) {
  const given = new Map();
  const giveOf = (block) => {
    const key = formatDrillKey(block);
    let share = given.get(key);
    if (share === undefined) {
      share = expanded && !expanded.has(key) ? 0 : handover(projectedSize(block, camera, edgePc, focalPx), block);
      given.set(key, share);
    }
    return share;
  };
  // What a block's children inherit: its own inheritance times what it gives.
  const handed = new Map();
  const handedDown = (block) => {
    const key = formatDrillKey(block);
    let found = handed.get(key);
    if (found === undefined) {
      const parent = drillParent(block);
      const inherited = parent === null ? 1 : handedDown(parent);
      found = inherited > 0 ? inherited * giveOf(block) : 0;
      handed.set(key, found);
    }
    return found;
  };
  const weights = new Map();
  return (block) => {
    const key = formatDrillKey(block);
    let found = weights.get(key);
    if (found === undefined) {
      const parent = drillParent(block);
      const inherited = parent === null ? 1 : handedDown(parent);
      found = inherited > 0 ? inherited * (1 - giveOf(block)) : 0;
      weights.set(key, found);
    }
    return found;
  };
}

// The weight of one block for a camera (planetgen.galaxy.nodecut.node_weight).
export function nodeWeight(block, camera, edgePc, focalPx, expanded) {
  return weightCalculator(camera, edgePc, focalPx, expanded)(block);
}

// The blocks of a block's children, kept per block (the ladder's counts are rounded floats, costly to redo).
const childLists = new Map();
export function childrenOf(block) {
  const key = formatDrillKey(block);
  let found = childLists.get(key);
  if (found === undefined) {
    if (childLists.size > 40000) childLists.clear();
    found = [];
    drillChildren(block).forEach((group) => group.blocks.forEach((child) => found.push(child)));
    childLists.set(key, found);
  }
  return found;
}

// A binary max-heap of open entries by projected size (ties: the smaller key first), as an array.
function above(p, q) {
  return p.size > q.size || (p.size === q.size && p.key < q.key);
}

function pushOpen(heap, entry) {
  heap.push(entry);
  let i = heap.length - 1;
  while (i > 0) {
    const up = (i - 1) >> 1;
    if (!above(heap[i], heap[up])) break;
    [heap[i], heap[up]] = [heap[up], heap[i]];
    i = up;
  }
}

function popLargest(heap) {
  const top = heap[0];
  const last = heap.pop();
  if (heap.length) {
    heap[0] = last;
    let i = 0;
    for (;;) {
      const left = 2 * i + 1;
      const right = left + 1;
      let best = i;
      if (left < heap.length && above(heap[left], heap[best])) best = left;
      if (right < heap.length && above(heap[right], heap[best])) best = right;
      if (best === i) break;
      [heap[i], heap[best]] = [heap[best], heap[i]];
      i = best;
    }
  }
  return top;
}

// A value to nine significant digits (planetgen.galaxy.nodecut._quantized: the order of the cut must not follow the last bits).
const quantized = (value) => (Number.isFinite(value) ? Number(value.toPrecision(9)) : value);

// The cut. `config`: { edgePc, focalPx, bounds: [[x, y, z], [x, y, z]] or null, halfAngleDeg, maxNodes, fetchFromPx }.
// Best first: the roots, then the block with the biggest projected size that is worth opening (its size at least
// fetchFromPx), until opening one more would pass maxNodes. Returns { keys (nearest first, ties by key), expanded
// (a Set of the keys opened) }.
export function nodeCut(camera, look, config) {
  const edgePc = config.edgePc;
  const focalPx = config.focalPx || DEFAULT_FOCAL_PX;
  const bounds = config.bounds || null;
  const halfAngleDeg = config.halfAngleDeg == null ? CONE_HALF_ANGLE_DEG : config.halfAngleDeg;
  const maxNodes = config.maxNodes || MAX_NODES;
  const fetchFrom = config.fetchFromPx == null ? FETCH_FROM_PX : config.fetchFromPx;
  let axis = [look[0] - camera[0], look[1] - camera[1], look[2] - camera[2]];
  const length = Math.sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
  const cone = halfAngleDeg < 180 && length > 0;
  if (cone) axis = axis.map((v) => v / length);
  const halfAngle = (halfAngleDeg * Math.PI) / 180;
  const found = [];
  const open = [];
  // The block's entry if it is kept (inside the bounds and the cone), else null.
  const consider = (block) => {
    const geo = geometry(block, edgePc);
    if (bounds) {
      for (let a = 0; a < 3; a++) {
        if (geo.hi[a] <= bounds[0][a] || geo.lo[a] > bounds[1][a]) return null;
      }
    }
    const distanceSq = distanceSqTo(geo, camera);
    if (cone && distanceSq > 0) {
      const offset = [0, 1, 2].map((a) => (geo.lo[a] + geo.hi[a]) / 2 - camera[a]);
      const reach = Math.sqrt(offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]);
      const radius = geo.size * SQRT3_HALF;
      if (reach > radius) {
        const along = (offset[0] * axis[0] + offset[1] * axis[1] + offset[2] * axis[2]) / reach;
        if (Math.acos(Math.max(-1, Math.min(1, along))) > halfAngle + Math.asin(radius / reach)) return null;
      }
    }
    const size = distanceSq <= 0 ? Infinity : (geo.size * focalPx) / Math.sqrt(distanceSq);
    return { block, key: formatDrillKey(block), distanceSq: quantized(distanceSq), size: quantized(size) };
  };
  const add = (entry) => {
    found.push(entry);
    if (entry.block.m > 1 && entry.size >= fetchFrom) pushOpen(open, entry);
  };
  roots().forEach((block) => {
    const entry = consider(block);
    if (entry) add(entry);
  });
  const expanded = new Set();
  while (open.length) {
    const entry = popLargest(open);
    const children = [];
    childrenOf(entry.block).forEach((child) => {
      const kept = consider(child);
      if (kept) children.push(kept);
    });
    if (found.length + children.length > maxNodes) break;
    expanded.add(entry.key);
    children.forEach(add);
  }
  found.sort((p, q) => (p.distanceSq - q.distanceSq) || (p.key < q.key ? -1 : p.key > q.key ? 1 : 0));
  return { keys: found.map((entry) => entry.key), expanded };
}

// The blocks of `cut.expanded` that may hand over to their children: those whose children in the cut are all in hand
// (`loaded(key)`). The page draws with these, so a block whose children are still on their way stays whole instead of
// leaving a hole.
export function settledExpanded(cut, loaded) {
  const kept = new Set(cut.keys);
  const settled = new Set();
  cut.expanded.forEach((key) => {
    const ready = childrenOf(parseNodeKey(key)).every((child) => {
      const childKey = formatDrillKey(child);
      return !kept.has(childKey) || loaded(childKey);
    });
    if (ready) settled.add(key);
  });
  return settled;
}

// A block as its node key `m.shell.band.wedge`.
export function formatNodeKey(block) {
  return formatDrillKey(block);
}

// A node key `m.shell.band.wedge` as a block.
export function parseNodeKey(key) {
  const parts = String(key).split(".").map(Number);
  return { m: parts[0], shell: parts[1], band: parts[2], wedge: parts[3] };
}
