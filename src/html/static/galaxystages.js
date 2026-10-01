// html/static/galaxystages.js
//
// The Galaxy Map drill-down's stages (docs/design/galaxy-drilldown-
// navigation.md, sections 4, 5 and 8.1), without any drawing: which
// blocks a stage holds and how many sectors each can hold, stage URLs,
// the breadcrumb, labels, and the camera path of a flight into a block.
// galaxymap3d.js draws and drives them. Like galaxyprisms.js it has no
// three.js import, so it runs under plain node for tests
// (src/tests/test_galaxystages.py).
//
// A stage is {at, slab}: `at` is the container (null for the galaxy, or a
// drill block {m, ring, wedge, slab} of 243, 27 or 3 sectors a side) and
// `slab` the child slab pulled out of it (null in the 3D view). So:
//
//   stage 1  {at: null, slab: null}   the galaxy, 3D, level-243 blocks
//   stage 2  {at: null, slab: s}      one level-243 slab, top-down
//   stage 3  {at: 243-block, slab: null}  its level-27 children, 3D
//   stage 4  {at: 243-block, slab: s}     one slab of them, top-down
//   stages 5-6 the same for a level-27 block (level-3 children),
//   stages 7-8 for a level-3 block, whose children are sectors and whose
//   child slabs are sector layers.

const VERSION_QUERY = new URL(import.meta.url).search;

const {
  drillBlockBounds, drillChainOf, drillChildren, drillParent, drillSlabs, drillSlotRange, drillWedgeCount,
  formatDrillKey, parseDrillKey, ringSectorCount, sectorAddressAt, sectorDrawable,
} = await import(`./galaxyprisms.js${VERSION_QUERY}`);

// The flight's zoom/pan trade-off (van Wijk and Nuij's rho), and its
// duration per unit of path length, clamped (section 5.3).
export const FLIGHT_RHO = 1.4;
export const FLIGHT_SPEED = 1.1;
export const FLIGHT_SECONDS = [0.5, 1.6];
// The 3D stages' tilt from straight down, and the pull-out's length.
export const STAGE_TILT_DEG = 35;
export const PULL_OUT_MS = 600;
// A top-down view fits the slab's footprint with this much margin.
export const FIT_MARGIN = 1.1;

// --- The galaxy's outline -----------------------------------------------------

// The outline the browser can work out for itself: for each sector ring,
// the highest |layer| drawn (-1 for none), from the analytic density bound
// (sectorDrawable). A ring's drawn layers are always -L..L: the bound only
// falls away from the plane. Binary search per ring, so a few thousand
// rings cost tens of thousands of density evaluations, once.
export function galaxyOutline(edgePc, shape, galaxyRadius) {
  const maxRing = Math.max(0, Math.ceil(galaxyRadius / edgePc) - 1);
  const extents = new Int32Array(maxRing + 1);
  // No ring reaches past the galaxy's radius in z either.
  const layerCap = Math.ceil(galaxyRadius / edgePc);
  for (let ring = 0; ring <= maxRing; ring++) {
    if (!shape) {
      extents[ring] = layerCap;
      continue;
    }
    if (!sectorDrawable(ring, 0, edgePc, shape)) {
      extents[ring] = -1;
      continue;
    }
    let lo = 0;
    let hi = layerCap;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (sectorDrawable(ring, mid, edgePc, shape)) lo = mid;
      else hi = mid - 1;
    }
    extents[ring] = lo;
  }
  let lastRing = -1;
  let maxLayer = -1;
  for (let ring = 0; ring <= maxRing; ring++) {
    if (extents[ring] >= 0) {
      lastRing = ring;
      maxLayer = Math.max(maxLayer, extents[ring]);
    }
  }
  return { extents: extents, maxRing: lastRing, maxLayer: maxLayer };
}

// How many sectors drill block `block` can hold under `outline`: per
// member ring, its slots in the block's wedge times its drawn layers in
// the block's layer range.
export function drillBlockTotal(block, outline) {
  const m = block.m;
  const half = (m - 1) / 2;
  const layerFirst = m === 1 ? block.slab : block.slab * m - half;
  const layerLast = m === 1 ? block.slab : block.slab * m + half;
  const ringFirst = m === 1 ? block.ring : block.ring * m;
  const ringLast = m === 1 ? block.ring : block.ring * m + m - 1;
  let total = 0;
  for (let i = ringFirst; i <= Math.min(ringLast, outline.extents.length - 1); i++) {
    const reach = outline.extents[i];
    if (reach < 0) continue;
    const layers = Math.min(layerLast, reach) - Math.max(layerFirst, -reach) + 1;
    if (layers <= 0) continue;
    const slots = drillSlotRange(block, i);
    total += Math.max(0, slots.last - slots.first + 1) * layers;
  }
  return total;
}

// --- Stages -------------------------------------------------------------------

export function stageNumber(stage) {
  const depth = stage.at ? { 243: 1, 27: 2, 3: 3 }[stage.at.m] : 0;
  return 1 + 2 * depth + (stage.slab == null ? 0 : 1);
}

export function isTopDown(stage) {
  return stage.slab != null;
}

export function sameBlock(a, b) {
  if (!a || !b) return a === b;
  return a.m === b.m && a.ring === b.ring && a.wedge === b.wedge && a.slab === b.slab;
}

export function sameStage(a, b) {
  return sameBlock(a.at, b.at) && (a.slab == null ? b.slab == null : a.slab === b.slab);
}

// The container's children grouped by child slab, lowest first, each
// block with its `bounds` and `total` (sectors it can hold); blocks that
// can hold none are left out, and so are slabs left empty. Without an
// outline (no shape), every child counts and `total` is null.
export function stageChildren(at, outline, edgePc) {
  const groups = drillChildren(at, outline.maxRing, outline.maxLayer);
  const out = [];
  groups.forEach(function (group) {
    const blocks = [];
    group.blocks.forEach(function (block) {
      const total = outline.shapeless ? null : drillBlockTotal(block, outline);
      if (total === 0) return;
      block.bounds = drillBlockBounds(block, edgePc);
      block.total = total;
      blocks.push(block);
    });
    if (blocks.length) out.push({ slab: group.slab, blocks: blocks });
  });
  return out;
}

// The stage one step up: a top-down view goes back to its container's 3D
// view; a 3D view to the top-down slab its container sits in.
export function parentStage(stage) {
  if (stage.slab != null) return { at: stage.at, slab: null };
  if (!stage.at) return null;
  return { at: drillParent(stage.at), slab: stage.at.slab };
}

// Every stage from the galaxy down to `stage`, top first.
export function stageChain(stage) {
  const chain = [];
  for (let s = stage; s; s = parentStage(s)) chain.unshift(s);
  return chain;
}

// The stage a sector is picked at (stage 8: its level-3 block's layer).
export function sectorStage(ring, layer, slot) {
  const chain = drillChainOf(ring, layer, slot);
  return { at: chain[2], slab: layer };
}

// --- Designations -------------------------------------------------------------

// provisional_sector_designation's hex code: ring << 33 | (layer + 4096)
// << 20 | slot.
export function sectorDesignation(ring, layer, slot) {
  return ((BigInt(ring) << 33n) | (BigInt(layer + 4096) << 20n) | BigInt(slot)).toString(16).toUpperCase();
}

// sectorDesignation's inverse, or null for anything that isn't one (a
// real slot of a real ring).
export function parseSectorDesignation(text) {
  const hex = String(text == null ? "" : text).trim();
  if (!/^[0-9A-Fa-f]{1,16}$/.test(hex)) return null;
  const packed = BigInt("0x" + hex);
  const ring = Number(packed >> 33n);
  const layer = Number((packed >> 20n) & 0x1FFFn) - 4096;
  const slot = Number(packed & 0xFFFFFn);
  if (slot >= ringSectorCount(ring)) return null;
  return { ring: ring, layer: layer, slot: slot };
}

// --- The address bar (section 9.3) -------------------------------------------

// What the address bar's text asks for:
// - {sector: {ring, layer, slot}} for a designation ("1AA5FDEC0C99"),
//   "312/-3/1042" or "ring 312 layer -3 slot 1042";
// - {sector, point: [x, y, z]} for "x, y, z" in pc (with or without
//   "pc"), the sector holding that point;
// - {name} for anything else (looked up by the server);
// - {problem} for an address that can't be a sector (a slot past its
//   ring's end, say), or empty text.
// Designations are 9 to 16 hex digits (the layer field alone puts one
// past 8), so a short name made of hex letters ("Bead") stays a name.
export function parseAddress(text, edgePc) {
  const raw = String(text == null ? "" : text).trim();
  if (!raw) return { problem: "Type a designation, ring/layer/slot, x, y, z in pc, or a name." };
  const number = "(-?\\d+)";
  let match = raw.match(new RegExp("^" + number + "\\s*/\\s*" + number + "\\s*/\\s*" + number + "$"))
    || raw.match(new RegExp("^ring\\s+" + number + "[\\s,]+layer\\s+" + number + "[\\s,]+slot\\s+" + number + "$", "i"));
  if (match) {
    const sector = { ring: Number(match[1]), layer: Number(match[2]), slot: Number(match[3]) };
    if (sector.ring < 0 || sector.slot < 0 || sector.slot >= ringSectorCount(sector.ring)) {
      return { problem: "Ring " + sector.ring + " has no slot " + sector.slot + "." };
    }
    return { sector: sector };
  }
  const real = "(-?\\d+(?:\\.\\d+)?)";
  match = raw.match(new RegExp("^\\(?\\s*" + real + "\\s*[,\\s]\\s*" + real + "\\s*[,\\s]\\s*" + real + "\\s*\\)?\\s*(?:pc)?$", "i"));
  if (match) {
    const point = [Number(match[1]), Number(match[2]), Number(match[3])];
    return { sector: sectorAddressAt(point[0], point[1], point[2], edgePc), point: point };
  }
  if (/^[0-9A-Fa-f]{9,16}$/.test(raw)) {
    const sector = parseSectorDesignation(raw);
    if (!sector) return { problem: "There is no sector " + raw.toUpperCase() + "." };
    return { sector: sector };
  }
  return { name: raw };
}

// --- URLs ---------------------------------------------------------------------

// A stage URL's query (section 8.1): "" for the galaxy, "?slab=s",
// "?at=m.ring.wedge.slab" or "?at=...&slab=s".
export function stageQuery(stage) {
  const parts = [];
  if (stage.at) parts.push("at=" + formatDrillKey(stage.at));
  if (stage.slab != null) parts.push("slab=" + stage.slab);
  return parts.length ? "?" + parts.join("&") : "";
}

// The stage a URL's query asks for: {stage, sector, problem}. `sector`
// ({ring, layer, slot}) is set by ?sector=<designation>, which opens its
// stage 8. Anything malformed opens the galaxy, with `problem` saying why.
// Whether the asked-for block or slab holds anything is checked by the
// caller against the outline (validStage).
export function parseStageQuery(search) {
  const params = new URLSearchParams(search || "");
  const galaxy = { at: null, slab: null };
  const designation = params.get("sector");
  if (designation != null) {
    const sector = parseSectorDesignation(designation);
    if (!sector) return { stage: galaxy, sector: null, problem: "There is no sector " + designation + "." };
    return { stage: sectorStage(sector.ring, sector.layer, sector.slot), sector: sector, problem: null };
  }
  let at = null;
  if (params.has("at")) {
    at = parseDrillKey(params.get("at"));
    if (!at) return { stage: galaxy, sector: null, problem: "There is no block " + params.get("at") + "." };
  }
  let slab = null;
  if (params.has("slab")) {
    const raw = params.get("slab");
    if (!/^-?\d+$/.test(raw)) return { stage: galaxy, sector: null, problem: "There is no slab " + raw + "." };
    slab = Number(raw);
  }
  return { stage: { at: at, slab: slab }, sector: null, problem: null };
}

// Whether `stage` holds anything under `outline`: its container can hold
// sectors (or is the galaxy) and its slab is one of the container's
// non-empty child slabs. Returns a problem string, or null when fine.
export function validStage(stage, outline, edgePc) {
  if (stage.at && !outline.shapeless && drillBlockTotal(stage.at, outline) === 0) {
    return "Block " + blockLabel(stage.at) + " is outside the galaxy.";
  }
  if (stage.slab != null) {
    const slabs = stageChildren(stage.at, outline, edgePc).map(function (g) { return g.slab; });
    if (slabs.indexOf(stage.slab) < 0) {
      return (stage.at && stage.at.m === 3 ? "Layer " : "Slab ") + stage.slab + " is empty here.";
    }
  }
  return null;
}

// --- Labels -------------------------------------------------------------------

export function formatInt(n) {
  return Math.round(n).toLocaleString("en-US");
}

// "7·14": a block's ring and wedge; a sector's "1,705·-20·3,225"
// (ring·layer·slot).
export function blockLabel(block) {
  if (block.m === 1) return [formatInt(block.ring), block.slab, formatInt(block.wedge)].join("·");
  return block.ring + "·" + block.wedge;
}

// What a child slab of `at` is called: a sector layer inside a level-3
// block, a slab otherwise.
export function slabNoun(at) {
  return at && at.m === 3 ? "Layer" : "Slab";
}

// The sector layers child slab `slab` of container `at` covers.
export function slabLayers(at, slab) {
  const childM = at ? at.m / (at.m === 3 ? 3 : 9) : 243;
  const half = (childM - 1) / 2;
  return { first: slab * childM - half, last: slab * childM + half, childM: childM };
}

// The breadcrumb for `stage`: one crumb per stage from the galaxy down,
// {stage, label}. A crumb's label names what that stage shows.
export function crumbs(stage) {
  return stageChain(stage).map(function (s, n, chain) {
    let label;
    if (!s.at && s.slab == null) label = "Galaxy";
    else if (s.slab == null) label = "Block " + blockLabel(s.at);
    else label = slabNoun(s.at) + " " + s.slab;
    return { stage: s, label: label, last: n === chain.length - 1 };
  });
}

// The stages beside a crumb's (its siblings, for the crumb's menu): the
// other non-empty slabs of the same container, or the blocks next to
// this one in the same parent slab. [{stage, label}].
export function crumbSiblings(stage, outline, edgePc) {
  if (stage.slab != null) {
    return stageChildren(stage.at, outline, edgePc).map(function (g) {
      return { stage: { at: stage.at, slab: g.slab }, label: slabNoun(stage.at) + " " + g.slab };
    });
  }
  if (!stage.at) return [];
  const parent = drillParent(stage.at);
  const group = stageChildren(parent, outline, edgePc).find(function (g) { return g.slab === stage.at.slab; });
  return (group ? group.blocks : []).map(function (block) {
    return { stage: { at: stripBlock(block), slab: null }, label: "Block " + blockLabel(block) };
  });
}

// A block without the bounds/total stageChildren adds.
export function stripBlock(block) {
  return { m: block.m, ring: block.ring, wedge: block.wedge, slab: block.slab };
}

// Degrees counterclockwise from +X, as the wedge labels read.
export function bearingRange(bounds) {
  const deg = function (rad) { return (((rad * 180) / Math.PI) % 360 + 360) % 360; };
  return [deg(bounds.t0), deg(bounds.t1) || 360];
}

// --- Framing and the flight ---------------------------------------------------

// The footprint of `blocks` (their bounds) in the plane: its center (the
// middle of the corners' bounding box) and the radius that reaches every
// corner from it, sampling each block's arcs.
export function footprint(blocks) {
  let minX = Infinity;
  let minY = Infinity;
  let maxX = -Infinity;
  let maxY = -Infinity;
  const points = [];
  blocks.forEach(function (block) {
    const b = block.bounds;
    const steps = Math.max(1, Math.ceil((b.t1 - b.t0) / (Math.PI / 32)));
    for (let k = 0; k <= steps; k++) {
      const t = b.t0 + ((b.t1 - b.t0) * k) / steps;
      [b.r0, b.r1].forEach(function (r) {
        points.push([r * Math.cos(t), r * Math.sin(t)]);
      });
    }
  });
  points.forEach(function (p) {
    minX = Math.min(minX, p[0]);
    minY = Math.min(minY, p[1]);
    maxX = Math.max(maxX, p[0]);
    maxY = Math.max(maxY, p[1]);
  });
  const center = [(minX + maxX) / 2, (minY + maxY) / 2];
  let radius = 0;
  points.forEach(function (p) {
    radius = Math.max(radius, Math.hypot(p[0] - center[0], p[1] - center[1]));
  });
  return { center: center, radius: radius };
}

// Van Wijk and Nuij's smooth zoom-and-pan from view (c0, w0) to (c1, w1):
// c a 2D center, w the visible width. Returns {S, at(s) -> {center, w}}
// over path length s in [0, S] (section 5.3).
export function flightPath(c0, w0, c1, w1, rho) {
  rho = rho || FLIGHT_RHO;
  const dx = c1[0] - c0[0];
  const dy = c1[1] - c0[1];
  const u = Math.hypot(dx, dy);
  if (u < 1e-9 * Math.max(w0, w1)) {
    const S = Math.abs(Math.log(w1 / w0)) / rho;
    const k = w1 >= w0 ? 1 : -1;
    return {
      S: S,
      at: function (s) {
        return { center: [c0[0], c0[1]], w: w0 * Math.exp(k * rho * s) };
      },
    };
  }
  const rho2 = rho * rho;
  const rho4 = rho2 * rho2;
  const b0 = (w1 * w1 - w0 * w0 + rho4 * u * u) / (2 * w0 * rho2 * u);
  const b1 = (w1 * w1 - w0 * w0 - rho4 * u * u) / (2 * w1 * rho2 * u);
  const r0 = Math.log(-b0 + Math.sqrt(b0 * b0 + 1));
  const r1 = Math.log(-b1 + Math.sqrt(b1 * b1 + 1));
  const S = (r1 - r0) / rho;
  return {
    S: S,
    at: function (s) {
      const w = (w0 * Math.cosh(r0)) / Math.cosh(rho * s + r0);
      const d = (w0 / rho2) * (Math.cosh(r0) * Math.tanh(rho * s + r0) - Math.sinh(r0));
      return { center: [c0[0] + (dx * d) / u, c0[1] + (dy * d) / u], w: w };
    },
  };
}

// A flight's duration in milliseconds for path length S.
export function flightMs(S) {
  return 1000 * Math.min(FLIGHT_SECONDS[1], Math.max(FLIGHT_SECONDS[0], S / FLIGHT_SPEED));
}

export function easeInOut(t) {
  return t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2;
}

export { drillBlockBounds, drillSlabs, drillWedgeCount };
