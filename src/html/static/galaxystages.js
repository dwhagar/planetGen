// html/static/galaxystages.js
//
// The Galaxy Map's stages (docs/design/galaxy-drilldown-navigation.md, as
// MAP.150 and MAP.151 changed it; docs/design/drilldown-region-sizes.md),
// without any drawing: which blocks a stage holds and how many sectors each
// can hold, stage URLs, the breadcrumb, labels, and the camera path of a
// flight into a block. galaxymap3d.js draws and drives them. Like
// galaxyprisms.js it has no three.js import, so it runs under plain node
// for tests (src/tests/test_galaxystages.py).
//
// A stage is {win}: see "Stages" below.

const VERSION_QUERY = new URL(import.meta.url).search;
const { formatNumber } = await import(`./numberformat.js${VERSION_QUERY}`);

const {
  DEFAULT_WINDOW, WINDOW_WIDTHS, drillBlockBounds, drillChainOf, drillChildren, drillParent, drillSlotRange,
  drillWedgeCount, formatDrillKey, ringSectorCount, sectorAddressAt, sectorDrawable, windowCells,
} = await import(`./galaxyprisms.js${VERSION_QUERY}`);

// The flight's zoom/pan trade-off (van Wijk and Nuij's rho), and its
// duration per unit of path length, clamped (section 5.3).
export const FLIGHT_RHO = 1.4;
export const FLIGHT_SPEED = 1.1;
export const FLIGHT_SECONDS = [0.5, 1.6];
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
//
// A stage is {win}: `win` is null for the galaxy (every level-243 block of it,
// a click on one opens the window one size down) or a window {centre, g}: `g`
// cells across around `centre`, a cell of any level (a block of 243, 27 or 3
// sectors a side, or a sector as {m: 1, ring, wedge: slot, slab: layer}).
// There is nothing to pick but a cell: a click on one takes the window one
// level down, centred on the cell clicked (a sector is opened instead), and
// Up takes it one level up, centred on the cell holding the centre (the
// galaxy above level 243). A window's data is /galaxy/window
// (planetgen/galaxy/window.py); the galaxy's is /galaxy/stage.

const TWO_PI = 2 * Math.PI;
const galaxyCache = new WeakMap();

export function galaxyStage() {
  return { win: null };
}

export function sameBlock(a, b) {
  if (!a || !b) return a === b;
  return a.m === b.m && a.ring === b.ring && a.wedge === b.wedge && a.slab === b.slab;
}

function sameWindow(a, b) {
  if (!a || !b) return !a && !b;
  return a.g === b.g && sameBlock(a.centre, b.centre);
}

export function sameStage(a, b) {
  return sameWindow(a.win, b.win);
}

// What a stage's data is read for: its window, or null (the galaxy).
export function holder(stage) {
  return stage.win;
}

// A block without the bounds/total a view adds.
export function stripBlock(block) {
  return { m: block.m, ring: block.ring, wedge: block.wedge, slab: block.slab };
}

function distinct(values) {
  return Array.from(new Set(values)).sort(function (p, q) { return p - q; });
}

// The galaxy's own blocks (the level-243 ones that can hold a sector), each
// with its `bounds` and `total` (sectors it can hold); without an outline (no
// shape), every block counts and `total` is null. Worked out once per outline.
export function galaxyBlocks(outline, edgePc) {
  let blocks = galaxyCache.get(outline);
  if (!blocks) {
    blocks = [];
    drillChildren(null, outline.maxRing, outline.maxLayer).forEach(function (group) {
      group.blocks.forEach(function (block) {
        const total = outline.shapeless ? null : drillBlockTotal(block, outline);
        if (total === 0) return;
        block.bounds = drillBlockBounds(block, edgePc);
        block.total = total;
        blocks.push(block);
      });
    });
    galaxyCache.set(outline, blocks);
  }
  return blocks;
}

// The window's blocks with their bounds and sector totals (those that can
// hold none are left out), as galaxyBlocks makes the galaxy's.
export function windowView(win, outline, edgePc) {
  const blocks = [];
  windowCells(win.centre, win.g).forEach(function (block) {
    const total = outline.shapeless ? null : drillBlockTotal(block, outline);
    if (total === 0) return;
    block.bounds = drillBlockBounds(block, edgePc);
    block.total = total;
    blocks.push(block);
  });
  // Bearings are measured within half a turn of the centre's, so a window
  // across bearing 0 doesn't span the whole turn.
  const mid = ((win.centre.wedge + 0.5) / drillWedgeCount(win.centre.m, win.centre.ring)) * TWO_PI;
  return { win: win, blocks: blocks, a0: mid - Math.PI, a1: mid + Math.PI };
}

// The choices of a view: each block is one, a click on it (windowInto) opens
// the window one size down.
function cellOptions(blocks) {
  return blocks.map(function (block) {
    return { pick: { kind: "cell", block: stripBlock(block) }, blocks: [block], a0: block.bounds.t0, a1: block.bounds.t1 };
  });
}

// Everything about `stage`: {stage (as given), view: {win, blocks, a0, a1},
// kind ("cell"), options (one per block), problem (when the window holds
// nothing of the galaxy)}.
export function resolveStage(stage, outline, edgePc) {
  if (!stage.win) {
    const blocks = galaxyBlocks(outline, edgePc);
    return { stage: stage, view: { win: null, blocks: blocks, a0: 0, a1: TWO_PI }, kind: "cell", options: cellOptions(blocks), problem: null };
  }
  const view = windowView(stage.win, outline, edgePc);
  if (!view.blocks.length) {
    return { stage: stage, problem: "Nothing of the galaxy lies around " + blockLabel(stage.win.centre) + "." };
  }
  return { stage: stage, view: view, kind: "cell", options: cellOptions(view.blocks), problem: null };
}

// The cell one size down that a click on `block` opens a window round: the
// one under `at` (a point {x, y} of the plane in parsecs, where the click
// landed, kept inside the block), or holding the block's middle without
// one. Null for a sector.
export function windowChildCentre(block, at, edgePc) {
  if (block.m === 1) return null;
  const childM = block.m === 3 ? 1 : block.m / 9;
  let radius = (block.ring + 0.5) * block.m;
  let turn = (block.wedge + 0.5) / drillWedgeCount(block.m, block.ring);
  if (at && edgePc) {
    const b = block.bounds || drillBlockBounds(block, edgePc);
    const inside = function (value, low, high) {
      const margin = (high - low) * 1e-6;
      return Math.min(high - margin, Math.max(low + margin, value));
    };
    radius = inside(Math.hypot(at.x, at.y), b.r0, b.r1) / (edgePc);
    // The bearing counted forward from the block's first, a click just
    // outside its sides going to the nearer.
    const span = b.t1 - b.t0;
    let along = (((Math.atan2(at.y, at.x) - b.t0) % TWO_PI) + TWO_PI) % TWO_PI;
    if (along > span) along = along - span < TWO_PI - along ? span : 0;
    turn = ((b.t0 + inside(along, 0, span)) / TWO_PI + 1) % 1;
  }
  const ring = Math.floor(radius / childM);
  const wedge = Math.min(drillWedgeCount(childM, ring) - 1, Math.floor(turn * drillWedgeCount(childM, ring)));
  return { m: childM, ring: ring, wedge: wedge, slab: (block.slab * block.m) / childM };
}

export function windowStage(centre, g) {
  return { win: { centre: stripBlock(centre), g: g || DEFAULT_WINDOW } };
}

// The window stage a click on `block` (at the point `at`, see
// windowChildCentre) opens, or null on a sector.
export function windowInto(stage, block, at, edgePc) {
  const centre = windowChildCentre(block, at, edgePc);
  return centre ? windowStage(centre, stage.win ? stage.win.g : DEFAULT_WINDOW) : null;
}

// The steps a visitor takes, galaxy first, to reach `stage` (resolveStage
// results): the window at each level on the way down to its centre.
export function stageSteps(stage, outline, edgePc) {
  const steps = [resolveStage(galaxyStage(), outline, edgePc)];
  if (!stage.win) return steps;
  const chain = [];
  for (let centre = stage.win.centre; centre; centre = drillParent(centre)) chain.push(centre);
  chain.reverse().forEach(function (centre) {
    const next = resolveStage(windowStage(centre, stage.win.g), outline, edgePc);
    if (!next.problem) steps.push(next);
  });
  return steps;
}

// The stage one step back from `stage` (the last different one on the way to
// it), or null at the galaxy.
export function parentStage(stage, outline, edgePc) {
  const steps = stageSteps(stage, outline, edgePc);
  for (let n = steps.length - 1; n >= 0; n--) {
    if (!sameStage(steps[n].stage, stage)) return steps[n].stage;
  }
  return null;
}

// Whether sector (ring, layer, slot) can hold anything under `outline`.
function sectorInGalaxy(ring, layer, slot, outline) {
  if (ring < 0 || slot < 0 || slot >= ringSectorCount(ring)) return false;
  if (outline.shapeless) return true;
  return ring < outline.extents.length && Math.abs(layer) <= outline.extents[ring];
}

// The stage that shows a sector among its neighbours (MAP.26): the window of
// sectors centred on it. Null when the sector is outside the galaxy.
export function sectorStage(ring, layer, slot, outline) {
  if (!sectorInGalaxy(ring, layer, slot, outline)) return null;
  return windowStage({ m: 1, ring: ring, wedge: slot, slab: layer }, DEFAULT_WINDOW);
}

// The smallest stage showing every one of `sectors` ({ring, layer, slot}),
// for the NAV course overlay (section 9.4): the window of the finest level
// whose DEFAULT_WINDOW cells around the course's middle hold the cell of
// every stop. An empty list, or a course too wide for a window of 243-sector
// blocks, gives the galaxy.
export function courseStage(sectors, outline, edgePc) {
  const galaxy = galaxyStage();
  if (!sectors || !sectors.length) return galaxy;
  const stops = sectors.filter(function (s) { return sectorInGalaxy(s.ring, s.layer, s.slot, outline); });
  if (stops.length !== sectors.length) return galaxy;
  const middle = [0, 0, 0];
  stops.forEach(function (s) {
    const b = drillBlockBounds({ m: 1, ring: s.ring, wedge: s.slot, slab: s.layer }, edgePc);
    const r = (b.r0 + b.r1) / 2;
    const t = (b.t0 + b.t1) / 2;
    middle[0] += (r * Math.cos(t)) / stops.length;
    middle[1] += (r * Math.sin(t)) / stops.length;
    middle[2] += ((b.z0 + b.z1) / 2) / stops.length;
  });
  const at = sectorAddressAt(middle[0], middle[1], middle[2], edgePc);
  const around = drillChainOf(at.ring, at.layer, at.slot).reverse();
  const chains = stops.map(function (s) { return drillChainOf(s.ring, s.layer, s.slot).reverse(); });
  for (let level = 0; level < around.length; level++) {
    const wanted = new Set(windowCells(around[level], DEFAULT_WINDOW).map(formatDrillKey));
    if (chains.every(function (chain) { return wanted.has(formatDrillKey(chain[level])); })) {
      return windowStage(around[level], DEFAULT_WINDOW);
    }
  }
  return galaxy;
}

// "m.ring.wedge.slab" of any level, a sector (m = 1) included, or null.
export function parseCellKey(key) {
  const parts = String(key == null ? "" : key).split(".");
  if (parts.length !== 4 || !parts.every(function (p) { return /^-?\d+$/.test(p); })) return null;
  const v = parts.map(Number);
  const cell = { m: v[0], ring: v[1], wedge: v[2], slab: v[3] };
  if ([243, 27, 3, 1].indexOf(cell.m) < 0 || cell.ring < 0 || !Number.isSafeInteger(cell.ring)) return null;
  if (cell.wedge < 0 || cell.wedge >= drillWedgeCount(cell.m, cell.ring)) return null;
  return cell;
}

export function windowLabel(win) {
  return (win.centre.m === 1 ? "Sector " : "Block ") + blockLabel(win.centre) + " (" + win.centre.m + ")";
}

// The boundary of the area `blocks` cover in the plane, as edges that
// all lie on the blocks' own sides (so a highlight starts and ends on
// real wedge lines, MAP.52): {radial: [{t, r0, r1}], circles: [{r, t0,
// t1}]}, bearings in radians measured within half a turn of `mid`. A
// side two blocks share is left out: radial edges at the ends of each
// ring's runs of blocks, and along each ring boundary the bearings
// covered on one side of it but not the other.
export function outlineEdges(blocks, mid) {
  const eps = 1e-9;
  const rings = new Map();
  blocks.forEach(function (block) {
    const b = block.bounds;
    const key = b.r0 + "/" + b.r1;
    if (!rings.has(key)) rings.set(key, { r0: b.r0, r1: b.r1, spans: [] });
    const t0 = mid + ((((b.t0 - mid + Math.PI) % TWO_PI) + TWO_PI) % TWO_PI) - Math.PI;
    rings.get(key).spans.push([t0, t0 + (b.t1 - b.t0)]);
  });
  const radial = [];
  // Per radius: the bearings covered just inside it and just outside it.
  const atRadius = new Map();
  function side(r, which, spans) {
    if (!atRadius.has(r)) atRadius.set(r, { inside: [], outside: [] });
    Array.prototype.push.apply(atRadius.get(r)[which], spans);
  }
  rings.forEach(function (ring) {
    const runs = mergeSpans(ring.spans, eps);
    runs.forEach(function (run) {
      if (run[1] - run[0] >= TWO_PI - eps) return;
      radial.push({ t: run[0], r0: ring.r0, r1: ring.r1 });
      radial.push({ t: run[1], r0: ring.r0, r1: ring.r1 });
    });
    side(ring.r1, "inside", runs);
    side(ring.r0, "outside", runs);
  });
  const circles = [];
  atRadius.forEach(function (sides, r) {
    if (r <= eps) return;
    spanXor(mergeSpans(sides.inside, eps), mergeSpans(sides.outside, eps), eps).forEach(function (span) {
      circles.push({ r: r, t0: span[0], t1: span[1] });
    });
  });
  return { radial: radial, circles: circles };
}

// [t0, t1] spans joined where they touch or overlap, in order.
function mergeSpans(spans, eps) {
  const sorted = spans.slice().sort(function (p, q) { return p[0] - q[0]; });
  const out = [];
  sorted.forEach(function (span) {
    const last = out[out.length - 1];
    if (last && span[0] <= last[1] + eps) last[1] = Math.max(last[1], span[1]);
    else out.push([span[0], span[1]]);
  });
  return out;
}

// The bearings in exactly one of two merged span lists, merged.
function spanXor(a, b, eps) {
  const cuts = distinct([].concat.apply([], a.concat(b)));
  const inside = function (list, t) { return list.some(function (s) { return t > s[0] && t < s[1]; }); };
  const out = [];
  for (let k = 0; k + 1 < cuts.length; k++) {
    const t = (cuts[k] + cuts[k + 1]) / 2;
    if (cuts[k + 1] - cuts[k] > eps && inside(a, t) !== inside(b, t)) out.push([cuts[k], cuts[k + 1]]);
  }
  return mergeSpans(out, eps);
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

// A stage URL's query (section 8.1): "" for the galaxy, else "?win=" with the
// centre's "m.ring.wedge.slab" and "&g=" the cells across (when not the default).
export function stageQuery(stage) {
  if (!stage.win) return "";
  const parts = ["win=" + formatDrillKey(stage.win.centre)];
  if (stage.win.g !== DEFAULT_WINDOW) parts.push("g=" + stage.win.g);
  return "?" + parts.join("&");
}

// The printed object ID of a system or body in a link (API.23: hex digits in up to three dashed groups,
// `0008000000-4000000-00D`; see objectref.js), and a body's `<kind>:<id>` reference.
const ID_PART = /^[0-9A-Fa-f]+(?:-[0-9A-Fa-f]+){0,2}$/;
const OBJECT_REF = /^(?:star|planet|moon|belt|comet):[0-9A-Fa-f]+(?:-[0-9A-Fa-f]+){0,2}$/;

// The stage a URL's query asks for: {stage, sector, open, problem}. `sector`
// ({ring, layer, slot}) is set by ?sector=<designation>, whose stage the
// caller works out (sectorStage needs the outline); `open` by &open=1 with
// it, the sector opened in place (MAP.66) and not only selected, and `system`
// (&system=<id>, with open, and &object=<ref> a body of it) the system of it
// opened in place too (MAP.125). ?win= (and &g=) names a window. Anything
// malformed, and the arcs, slabs and blocks of the links before MAP.150,
// opens the galaxy; `problem` says why for the first. Whether a window holds
// anything is the caller's to check against the outline (resolveStage).
export function parseStageQuery(search) {
  const params = new URLSearchParams(search || "");
  const galaxy = galaxyStage();
  const designation = params.get("sector");
  if (designation != null) {
    const sector = parseSectorDesignation(designation);
    if (!sector) return { stage: galaxy, sector: null, problem: "There is no sector " + designation + "." };
    // The ID stays text, in upper case as objectref.js writes it.
    const system = ID_PART.test(params.get("system") || "") ? params.get("system").toUpperCase() : null;
    const wantedObject = params.get("object") || "";
    const object = OBJECT_REF.test(wantedObject)
      ? wantedObject.slice(0, wantedObject.indexOf(":") + 1) + wantedObject.slice(wantedObject.indexOf(":") + 1).toUpperCase()
      : null;
    return { stage: galaxy, sector: sector, open: params.get("open") === "1", system: system, object: object, problem: null };
  }
  if (params.has("win")) {
    const centre = parseCellKey(params.get("win"));
    if (!centre) return { stage: galaxy, sector: null, problem: "There is no cell " + params.get("win") + "." };
    const g = params.has("g") ? Number(params.get("g")) : DEFAULT_WINDOW;
    if (WINDOW_WIDTHS.indexOf(g) < 0) return { stage: galaxy, sector: null, problem: "A window cannot be " + params.get("g") + " cells across." };
    return { stage: windowStage(centre, g), sector: null, problem: null };
  }
  return { stage: galaxy, sector: null, problem: null };
}

// --- Labels -------------------------------------------------------------------

// A label's whole number, never scientific (ring, wedge and slot numbers
// are addresses, not amounts).
export function formatInt(n) {
  return Math.round(n).toLocaleString("en-US");
}

// A count of sectors: scientific past 4 digits (UX.20).
export function formatCount(n) {
  return formatNumber(Math.round(n));
}

// "7·14": a block's ring and wedge; a sector's "1,705·-20·3,225"
// (ring·layer·slot).
export function blockLabel(block) {
  if (block.m === 1) return [formatInt(block.ring), block.slab, formatInt(block.wedge)].join("·");
  return block.ring + "·" + block.wedge;
}

// The breadcrumb for `stage`: one crumb per step from the galaxy down,
// {stage, label, last}.
export function crumbs(stage, outline, edgePc) {
  const steps = stageSteps(stage, outline, edgePc);
  return steps.map(function (r, n) {
    const s = r.stage;
    return { stage: s, label: s.win ? windowLabel(s.win) : "Galaxy", last: n === steps.length - 1 };
  });
}

// Degrees counterclockwise from +X, as the wedge labels read.
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

export { drillBlockBounds, drillParent, drillWedgeCount, formatDrillKey };
