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
// A stage is {at, picks}: see "Stages" below.

const VERSION_QUERY = new URL(import.meta.url).search;
const { formatNumber } = await import(`./numberformat.js${VERSION_QUERY}`);

const {
  drillBlockBounds, drillChainOf, drillChildren, drillParent, drillSlabs, drillSlotRange, drillWedgeCount,
  formatDrillKey, parseDrillKey, ringSectorCount, sectorAddressAt, sectorDrawable,
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
// A stage is {at, picks}: `at` is the container (null for the galaxy, or a
// drill block {m, ring, wedge, slab} of 243, 27 or 3 sectors a side) and
// `picks` what has been picked inside it so far, in order. The picks
// alternate, as Boss laid them out (MAP.19: "arc, layer, region, layer,
// region, ..., layer, sector"):
// - {kind: "arc", n}: the galaxy's first pick (MAP.85), a piece of the
//   disk about 40 degrees of bearing by a third of its radius, the whole
//   height of the disk: arc n is radial band floor(n / ARCS_PER_TURN)
//   (inner, middle, outer) and bearing bin n % ARCS_PER_TURN. Its blocks
//   are the level-243 blocks whose middle bearing falls in the bin, so
//   its sides run along their wedge lines (MAP.52);
// - {kind: "quadrant", n}: the galaxy's first pick before MAP.85, a
//   quarter of the disk; still read from older links, never offered;
// - {kind: "layer", lo, hi}: the container's child slabs lo to hi -- the
//   slice; pickOptions offers thirds while more than three are left (the
//   keyboard and links use these), and the slab slider beside the map
//   (galaxystageview.js) takes any one slab;
// - {kind: "region", n}: an arc of the ring band in view, one of up to
//   3 x 3 (a third of its rings across, a third of its arc along).
// After a region comes a layer and after a layer a region, while each
// still has something to split; the first pick inside a block is a
// layer. A pick that leaves one block of one slab enters that block: it
// becomes the container, and its children are picked the same way.
// Inside a level-3 block the children are sectors, and the last pick is
// a sector of one layer.

// Up to this many slices or arcs per pick (so a pick is at least about a
// ninth of the view: big targets, MAP.19).
export const PICK_SPLIT = 3;
// The galaxy's arcs (MAP.85): this many bearing bins a turn, in
// ARC_BANDS bands of radius. 8 bins (45 degrees) is the count nearest
// Boss's 40 degrees whose lines are wedge lines of every level-243 block
// ring past the core (their wedge counts are multiples of 8), so an arc's
// sides are straight lines from its inner edge to its outer one.
export const ARCS_PER_TURN = 8;
export const ARC_BANDS = 3;
const ARC_BAND_NAMES = ["inner", "middle", "outer"];

const TWO_PI = 2 * Math.PI;
const childCache = new WeakMap();

export function sameBlock(a, b) {
  if (!a || !b) return a === b;
  return a.m === b.m && a.ring === b.ring && a.wedge === b.wedge && a.slab === b.slab;
}

export function samePick(a, b) {
  if (a.kind !== b.kind) return false;
  return a.kind === "layer" ? a.lo === b.lo && a.hi === b.hi : a.n === b.n;
}

export function sameStage(a, b) {
  return sameBlock(a.at, b.at) && a.picks.length === b.picks.length
    && a.picks.every(function (pick, n) { return samePick(pick, b.picks[n]); });
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

// stageChildren, worked out once per outline and container.
export function childrenOf(at, outline, edgePc) {
  let byKey = childCache.get(outline);
  if (!byKey) {
    byKey = new Map();
    childCache.set(outline, byKey);
  }
  const key = (at ? formatDrillKey(at) : "galaxy") + "@" + edgePc;
  let groups = byKey.get(key);
  if (!groups) {
    groups = stageChildren(at, outline, edgePc);
    byKey.set(key, groups);
  }
  return groups;
}

function angleFrom(angle, a0) {
  return ((angle - a0) % TWO_PI + TWO_PI) % TWO_PI;
}

function midAngle(block) {
  return (block.bounds.t0 + block.bounds.t1) / 2;
}

function distinct(values) {
  return Array.from(new Set(values)).sort(function (p, q) { return p - q; });
}

function slabsOf(blocks) {
  return distinct(blocks.map(function (b) { return b.slab; }));
}

function columnsOf(blocks) {
  return new Set(blocks.map(function (b) { return b.ring + "/" + b.wedge; }));
}

// The sectors of a level-27 block whose sectors all lie in one layer, or
// null: such a block is only one sector tall, so its view skips the
// level-3 blocks and shows the sectors themselves (Boss, 2026-10-01).
export function thinSectors(at, outline, edgePc) {
  if (!at || at.m !== 27 || outline.shapeless) return null;
  const groups = childrenOf(at, outline, edgePc);
  if (groups.length !== 1) return null;
  const sectors = [];
  const layers = new Set();
  groups[0].blocks.forEach(function (block) {
    childrenOf(stripBlock(block), outline, edgePc).forEach(function (g) {
      layers.add(g.slab);
      Array.prototype.push.apply(sectors, g.blocks);
    });
  });
  return layers.size === 1 && sectors.length ? sectors : null;
}

// The container's whole view: every child (or a thin block's sectors),
// and the bearings it spans.
function containerView(at, outline, edgePc) {
  let blocks = thinSectors(at, outline, edgePc);
  if (!blocks) {
    blocks = [];
    childrenOf(at, outline, edgePc).forEach(function (g) { Array.prototype.push.apply(blocks, g.blocks); });
  }
  let a0 = 0;
  let a1 = TWO_PI;
  if (at) {
    const b = drillBlockBounds(at, edgePc);
    a0 = b.t0;
    a1 = b.t1;
  }
  return { at: at, blocks: blocks, a0: a0, a1: a1, hadFirstPick: false };
}

// What can be picked next in `view` after `picks`: "arc", "layer",
// "region", or null when one block (of one slab) is left.
export function nextPickKind(view, picks) {
  const slabs = slabsOf(view.blocks).length;
  const columns = columnsOf(view.blocks).size;
  if (!view.at && !view.hadFirstPick && columns > 1) return "arc";
  const last = picks.length ? picks[picks.length - 1].kind : null;
  if (last !== "layer" && slabs > 1) return "layer";
  if (columns > 1) return "region";
  if (slabs > 1) return "layer";
  return null;
}

// The choices for a pick of `kind` in `view`: [{pick, blocks, a0, a1}],
// each holding at least one block. Layers come lowest first.
export function pickOptions(view, kind) {
  const out = [];
  if (kind === "layer") {
    const slabs = slabsOf(view.blocks);
    const parts = Math.min(PICK_SPLIT, slabs.length);
    for (let k = 0; k < parts; k++) {
      const some = slabs.slice(Math.floor((k * slabs.length) / parts), Math.floor(((k + 1) * slabs.length) / parts));
      const lo = some[0];
      const hi = some[some.length - 1];
      out.push({
        pick: { kind: "layer", lo: lo, hi: hi },
        blocks: view.blocks.filter(function (b) { return b.slab >= lo && b.slab <= hi; }),
        a0: view.a0, a1: view.a1,
      });
    }
    return out;
  }
  const rings = distinct(view.blocks.map(function (b) { return b.ring; }));
  const perRing = new Map();
  columnsOf(view.blocks).forEach(function (key) {
    const ring = Number(key.split("/")[0]);
    perRing.set(ring, (perRing.get(ring) || 0) + 1);
  });
  let widest = 0;
  perRing.forEach(function (n) { widest = Math.max(widest, n); });
  let bands = Math.min(PICK_SPLIT, rings.length);
  let arcs = Math.min(PICK_SPLIT, widest);
  if (kind === "quadrant") {
    bands = 1;
    arcs = 4;
  } else if (kind === "arc") {
    bands = Math.min(ARC_BANDS, rings.length);
    arcs = ARCS_PER_TURN;
  }
  const span = view.a1 - view.a0;
  const byOption = new Map();
  view.blocks.forEach(function (block) {
    const band = Math.floor((rings.indexOf(block.ring) * bands) / rings.length);
    const arc = Math.min(arcs - 1, Math.floor((angleFrom(midAngle(block), view.a0) / span) * arcs));
    const n = band * arcs + arc;
    if (!byOption.has(n)) {
      byOption.set(n, { pick: { kind: kind, n: n }, arc: arc, blocks: [] });
    }
    byOption.get(n).blocks.push(block);
  });
  Array.from(byOption.keys()).sort(function (p, q) { return p - q; }).forEach(function (n) {
    const option = byOption.get(n);
    // A choice reaches as far as its blocks do, not just its even share
    // of the view: a block is sorted by its middle, and near the core a
    // wedge (between meridians that run all the way in) is wider than a
    // share, so the choice zooms into the wedge it highlights (Boss,
    // 2026-10-01).
    const share = blockSpan(option.blocks, view.a0 + (span * (option.arc + 0.5)) / arcs);
    option.a0 = share.a0;
    option.a1 = share.a1;
    out.push(option);
  });
  return out;
}

// The bearings `blocks` cover, a0 to a1 (a1 - a0 under a turn), each
// measured within half a turn of `mid`.
function blockSpan(blocks, mid) {
  let a0 = Infinity;
  let a1 = -Infinity;
  blocks.forEach(function (block) {
    const b = block.bounds;
    const t0 = mid + ((((b.t0 - mid + Math.PI) % TWO_PI) + TWO_PI) % TWO_PI) - Math.PI;
    a0 = Math.min(a0, t0);
    a1 = Math.max(a1, t0 + (b.t1 - b.t0));
  });
  return { a0: a0, a1: Math.min(a1, a0 + TWO_PI) };
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

// `view` after `pick`, or null when the pick isn't one it could take.
// A layer pick may name any range of the view's slabs that narrows it
// (an older ?slab= link names one slab of nine).
function applyPick(view, pick) {
  if (pick.kind === "layer") {
    const slabs = slabsOf(view.blocks);
    const blocks = view.blocks.filter(function (b) { return b.slab >= pick.lo && b.slab <= pick.hi; });
    if (!(pick.lo <= pick.hi) || slabs.length < 2 || !blocks.length) return null;
    return { at: view.at, blocks: blocks, a0: view.a0, a1: view.a1, hadFirstPick: view.hadFirstPick };
  }
  if ((pick.kind === "arc" || pick.kind === "quadrant") && (view.at || view.hadFirstPick)) return null;
  if (pick.kind === "region" && columnsOf(view.blocks).size < 2) return null;
  const option = pickOptions(view, pick.kind).find(function (o) { return o.pick.n === pick.n; });
  if (!option) return null;
  return { at: view.at, blocks: option.blocks, a0: option.a0, a1: option.a1, hadFirstPick: true };
}

function pickText(pick) {
  if (pick.kind === "layer") return "layer " + (pick.lo === pick.hi ? pick.lo : pick.lo + " to " + pick.hi);
  if (pick.kind === "arc") return "arc " + pickToken(pick).slice(1);
  return pick.kind + " " + pick.n;
}

// Everything about `stage`: {stage (as given), view: {at, blocks, a0,
// a1}, kind (the next pick's, or null), options (pickOptions for it),
// sector (the one sector left, if that's all there is), problem (when a
// pick or the container doesn't fit the galaxy)}.
export function resolveStage(stage, outline, edgePc) {
  if (stage.at && !outline.shapeless && drillBlockTotal(stage.at, outline) === 0) {
    return { stage: stage, problem: "Block " + blockLabel(stage.at) + " is outside the galaxy." };
  }
  let view = containerView(stage.at, outline, edgePc);
  for (let n = 0; n < stage.picks.length; n++) {
    const next = applyPick(view, stage.picks[n]);
    if (!next) return { stage: stage, problem: "There is no " + pickText(stage.picks[n]) + " here." };
    view = next;
  }
  const kind = view.blocks.length ? nextPickKind(view, stage.picks) : null;
  const options = kind ? pickOptions(view, kind) : [];
  const sector = !kind && view.blocks.length === 1 && view.blocks[0].m === 1 ? view.blocks[0] : null;
  return { stage: stage, view: view, kind: kind, options: options, sector: sector, problem: null };
}

// `stage` carried on as far as it goes without a choice: a pick with only
// one option is taken, and one block of one slab is entered. Returns the
// resolveStage of where it ends.
export function settleStage(stage, outline, edgePc) {
  let resolved = resolveStage(stage, outline, edgePc);
  for (let guard = 0; guard < 200 && !resolved.problem; guard++) {
    const view = resolved.view;
    if (resolved.kind && resolved.options.length === 1) {
      const picks = resolved.stage.picks.concat([resolved.options[0].pick]);
      resolved = resolveStage({ at: resolved.stage.at, picks: picks }, outline, edgePc);
    } else if (!resolved.kind && view.blocks.length === 1 && view.blocks[0].m > 1) {
      resolved = resolveStage({ at: stripBlock(view.blocks[0]), picks: [] }, outline, edgePc);
    } else {
      break;
    }
  }
  return resolved;
}

// `block` at size `m`, or its ancestor there ({m: 1, ring, wedge: slot,
// slab: layer} for a sector).
function ancestorAt(block, m) {
  if (block.m === m) return block;
  if (block.m === 1) {
    return drillChainOf(block.ring, block.slab, block.wedge).find(function (b) { return b.m === m; }) || null;
  }
  let b = block;
  while (b && b.m < m) b = drillParent(b);
  return b && b.m === m ? b : null;
}

// The way from the galaxy toward `target` (a drill block, or a sector as
// {m: 1, ring, wedge: slot, slab: layer}): one settled step per pick,
// until `done(resolved)` says stop. Returns the steps (resolveStage
// results), or null when `target` isn't in the galaxy.
function walkToward(target, outline, edgePc, done) {
  let resolved = settleStage({ at: null, picks: [] }, outline, edgePc);
  const steps = [resolved];
  for (let guard = 0; guard < 200; guard++) {
    if (resolved.problem) return null;
    if (done(resolved)) return steps;
    if (!resolved.kind) return null;
    const blocks = resolved.view.blocks;
    const goal = blocks.length ? ancestorAt(target, blocks[0].m) : null;
    if (!goal) return null;
    const option = resolved.options.find(function (o) {
      if (o.pick.kind === "layer") return goal.slab >= o.pick.lo && goal.slab <= o.pick.hi;
      return o.blocks.some(function (b) { return b.ring === goal.ring && b.wedge === goal.wedge; });
    });
    if (!option) return null;
    const picks = resolved.stage.picks.concat([option.pick]);
    resolved = settleStage({ at: resolved.stage.at, picks: picks }, outline, edgePc);
    steps.push(resolved);
  }
  return null;
}

// Whether `resolved` shows sectors of one layer, `layer`.
function sectorLayerShown(resolved, layer) {
  const blocks = resolved.view ? resolved.view.blocks : [];
  return blocks.length > 0 && blocks[0].m === 1 && blocks.every(function (b) { return b.slab === layer; });
}

// The stage that shows a sector to pick (the sector level of the slice
// holding it, MAP.26): its level-3 block with that one layer picked. Null
// when the sector is outside the galaxy.
export function sectorStage(ring, layer, slot, outline, edgePc) {
  const target = { m: 1, ring: ring, wedge: slot, slab: layer };
  const steps = walkToward(target, outline, edgePc, function (r) {
    return sectorLayerShown(r, layer) && r.view.blocks.some(function (b) { return b.ring === ring && b.wedge === slot; });
  });
  return steps ? steps[steps.length - 1].stage : null;
}

// The steps a visitor takes, galaxy first, to reach `stage` ([settled
// resolveStage results]): into each container the way its own picks
// would go, then `stage`'s own picks.
export function stageSteps(stage, outline, edgePc) {
  let steps = [];
  if (stage.at) {
    steps = walkToward(stage.at, outline, edgePc, function (r) { return sameBlock(r.stage.at, stage.at); }) || [];
    steps.pop();
  }
  for (let n = 0; n <= stage.picks.length; n++) {
    const next = settleStage({ at: stage.at, picks: stage.picks.slice(0, n) }, outline, edgePc);
    if (next.problem) break;
    if (!steps.length || !sameStage(next.stage, steps[steps.length - 1].stage)) steps.push(next);
  }
  return steps;
}

// The stage one step back from `stage` (the last different one on the
// way to it), or null at the galaxy.
export function parentStage(stage, outline, edgePc) {
  const steps = stageSteps(stage, outline, edgePc);
  for (let n = steps.length - 1; n >= 0; n--) {
    if (!sameStage(steps[n].stage, stage)) return steps[n].stage;
  }
  return null;
}

// The smallest stage showing every one of `sectors` ({ring, layer,
// slot}), for the NAV course overlay (section 9.4): as far along the way
// to each as they all go together. An empty list, or sectors whose ways
// part at once, give the galaxy.
export function courseStage(sectors, outline, edgePc) {
  const galaxy = settleStage({ at: null, picks: [] }, outline, edgePc).stage;
  if (!sectors || !sectors.length) return galaxy;
  const ways = sectors.map(function (s) {
    const target = { m: 1, ring: s.ring, wedge: s.slot, slab: s.layer };
    return walkToward(target, outline, edgePc, function (r) { return sectorLayerShown(r, s.layer); });
  });
  if (ways.some(function (w) { return !w; })) return galaxy;
  let shared = galaxy;
  for (let n = 0; n < ways[0].length; n++) {
    const here = ways[0][n].stage;
    if (!ways.every(function (w) { return w[n] && sameStage(w[n].stage, here); })) break;
    shared = here;
  }
  return shared;
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

// A pick as it reads in a URL: "a1.90" (the arc of band 1, the middle
// third, starting at bearing 90 degrees), "r4" (region), "L-1" (one slab
// or layer) or "L-4~-2" (a range of them); "q1" (quadrant) from older
// links.
export function pickToken(pick) {
  if (pick.kind === "layer") return "L" + pick.lo + (pick.hi === pick.lo ? "" : "~" + pick.hi);
  if (pick.kind === "arc") {
    const band = Math.floor(pick.n / ARCS_PER_TURN);
    return "a" + band + "." + ((pick.n % ARCS_PER_TURN) * 360) / ARCS_PER_TURN;
  }
  return (pick.kind === "quadrant" ? "q" : "r") + pick.n;
}

export function parsePickToken(token) {
  let match = /^a(\d+)\.(\d+)$/.exec(token);
  if (match) {
    const bin = (Number(match[2]) * ARCS_PER_TURN) / 360;
    if (!Number.isInteger(bin) || bin >= ARCS_PER_TURN) return null;
    return { kind: "arc", n: Number(match[1]) * ARCS_PER_TURN + bin };
  }
  match = /^([qr])(\d+)$/.exec(token);
  if (match) return { kind: match[1] === "q" ? "quadrant" : "region", n: Number(match[2]) };
  match = /^L(-?\d+)(?:~(-?\d+))?$/.exec(token);
  if (match) return { kind: "layer", lo: Number(match[1]), hi: Number(match[2] != null ? match[2] : match[1]) };
  return null;
}

// A stage URL's query (section 8.1): "" for the galaxy, else
// "?at=m.ring.wedge.slab" and "p=" with its picks, comma-separated.
export function stageQuery(stage) {
  const parts = [];
  if (stage.at) parts.push("at=" + formatDrillKey(stage.at));
  if (stage.picks.length) parts.push("p=" + stage.picks.map(pickToken).join(","));
  return parts.length ? "?" + parts.join("&") : "";
}

// The stage a URL's query asks for: {stage, sector, problem}. `sector`
// ({ring, layer, slot}) is set by ?sector=<designation>, whose stage the
// caller works out (sectorStage needs the outline). An older link's
// ?slab=s reads as a pick of that one slab, before any others. Anything
// malformed opens the galaxy, with `problem` saying why. Whether the
// picks fit is the caller's to check against the outline (resolveStage).
export function parseStageQuery(search) {
  const params = new URLSearchParams(search || "");
  const galaxy = { at: null, picks: [] };
  const designation = params.get("sector");
  if (designation != null) {
    const sector = parseSectorDesignation(designation);
    if (!sector) return { stage: galaxy, sector: null, problem: "There is no sector " + designation + "." };
    return { stage: galaxy, sector: sector, problem: null };
  }
  let at = null;
  if (params.has("at")) {
    at = parseDrillKey(params.get("at"));
    if (!at) return { stage: galaxy, sector: null, problem: "There is no block " + params.get("at") + "." };
  }
  const picks = [];
  if (params.has("slab")) {
    const raw = params.get("slab");
    if (!/^-?\d+$/.test(raw)) return { stage: galaxy, sector: null, problem: "There is no slab " + raw + "." };
    picks.push({ kind: "layer", lo: Number(raw), hi: Number(raw) });
  }
  const tokens = (params.get("p") || "").split(",").filter(Boolean);
  for (let n = 0; n < tokens.length; n++) {
    const pick = parsePickToken(tokens[n]);
    if (!pick) return { stage: galaxy, sector: null, problem: "There is no pick " + tokens[n] + "." };
    picks.push(pick);
  }
  return { stage: { at: at, picks: picks }, sector: null, problem: null };
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

// Bearings a0 to a1 as "92°–97°" ("92.5°–93.1°" for a narrow arc).
export function arcLabel(a0, a1) {
  const narrow = ((a1 - a0) * 180) / Math.PI < 10;
  // A bearing a hair under 0 (a block edge after floating-point sums)
  // reads 0, not 360.
  const deg = function (rad) {
    const d = (((rad * 180) / Math.PI) % 360 + 360) % 360;
    const text = narrow ? d.toFixed(1) : String(Math.round(d) % 360);
    return Number(text) >= 360 ? (narrow ? "0.0" : "0") : text;
  };
  const end = a1 - a0 >= TWO_PI - 1e-9 || Number(deg(a1)) === 0 ? "360" : deg(a1);
  return deg(a0) + "°–" + end + "°";
}

// What a pick is called in the breadcrumb, given the view it led to.
export function pickLabel(pick, at, view) {
  if (pick.kind === "layer") {
    const noun = slabNoun(at);
    return pick.lo === pick.hi ? noun + " " + pick.lo : noun + "s " + pick.lo + " to " + pick.hi;
  }
  if (pick.kind === "arc") {
    const band = ARC_BAND_NAMES[Math.floor(pick.n / ARCS_PER_TURN)];
    return "Arc " + arcLabel(view.a0, view.a1) + (band ? " (" + band + ")" : "");
  }
  return (pick.kind === "quadrant" ? "Quarter " : "Arc ") + arcLabel(view.a0, view.a1);
}

// The breadcrumb for `stage`: one crumb per step from the galaxy down,
// {stage, label, last}. A crumb's label names the pick that led there, or
// the block entered.
export function crumbs(stage, outline, edgePc) {
  const steps = stageSteps(stage, outline, edgePc);
  return steps.map(function (r, n) {
    const s = r.stage;
    let label;
    if (!s.picks.length) label = s.at ? "Block " + blockLabel(s.at) : "Galaxy";
    else label = pickLabel(s.picks[s.picks.length - 1], s.at, r.view);
    return { stage: s, label: label, last: n === steps.length - 1 };
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

export { drillBlockBounds, drillSlabs, drillWedgeCount, formatDrillKey };
