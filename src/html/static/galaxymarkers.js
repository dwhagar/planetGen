// html/static/galaxymarkers.js
//
// The Galaxy Map's marker layer, the logic half (MAP.219): a registry of marker kinds, each a small glyph the map
// draws over a notable thing so it can be found from far away (black holes first; habitable worlds and chosen
// classes are MAP.220, stations and vessels MAP.221). The drawing half is galaxymarkerlayer.js; this file has no
// DOM at import and no three.js, so it runs under plain node for tests (tests/js/galaxymarkers.test.mjs).
//
// A kind is
//   id          the layer id in the map's layer registry (maplayers.js), also what the URL lists ("marker-black-hole");
//   label       the Menu button's text and the legend line;
//   glyph       the shape drawn (a key of GLYPHS), or
//   icon        the point-of-interest symbol drawn (an id of poiicons.js, MAP.320, its file under static/poi/), a kind has
//               one of the two; the symbol takes the colour of its token (`item.tier`, the habitability score's tier,
//               picks one for the habitable and the inhabited world);
//   color       the glyph's color (CSS), also the legend swatch; an icon kind has none, its token colours it;
//   maxBlock    the largest block size (shells) whose node lists this kind completely (a node of a bigger block does not
//               list it, so the kind is not drawn from it; default: every node does);
//   defaultOn   shown until the visitor says otherwise (default true);
//   source      `function (node) -> [item]`: the kind's items in one node. An item has `key` (unique within the
//               kind), `x`, `y`, `z` (galaxy frame, parsecs) and anything the kind needs to describe it;
//   rank        `function (item) -> number`: higher is kept first when the budget cuts the list;
//   extentPc    optional `function (item) -> parsecs`: the item's own size, for a kind whose items are extended
//               (a nebula): it is hidden while that is under MIN_EXTENT_PX on screen;
//   describe    `function (item) -> string`: the hover text.
//
// What shows is chosen by `selectMarkers`: only kinds that are on, only items in view, and no more than the
// budget (MARKER_BUDGET) in all, the kinds sharing it evenly and handing on what one cannot use.

const VERSION_QUERY = new URL(import.meta.url).search;
const { poiIcon, tokenOf, TIER_NAMES } = await import(`./poiicons.js${VERSION_QUERY}`);

// The most markers drawn at once (MAP.116's rule for the stars: a crowded view must stay readable). Each kind gets
// an even share of it; a kind that has less hands the rest on.
export const MARKER_BUDGET = 240;

// A marker of a kind with its own size (`extentPc`) is hidden while it would be under this many pixels across.
export const MIN_EXTENT_PX = 2;

// Pixels across a drawn glyph.
export const MARKER_PX = 18;

// The shapes a kind draws, each a function `(ctx, size, color, outline)` that draws into a canvas context `size`
// pixels square. The black hole is the one drawn here; the other kinds are point-of-interest symbols (`icon`, MAP.320),
// vector files, not code.
export const GLYPHS = {
  // A black hole: a dark disc in a ring.
  "black-hole": function (ctx, size, color, outline) {
    const c = size / 2;
    ctx.lineJoin = "round";
    ctx.beginPath();
    ctx.arc(c, c, c * 0.78, 0, Math.PI * 2);
    ctx.lineWidth = size * 0.2;
    ctx.strokeStyle = outline;
    ctx.stroke();
    ctx.lineWidth = size * 0.1;
    ctx.strokeStyle = color;
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(c, c, c * 0.34, 0, Math.PI * 2);
    ctx.fillStyle = outline;
    ctx.fill();
    ctx.lineWidth = size * 0.06;
    ctx.strokeStyle = color;
    ctx.stroke();
  },
};

// The markers of black holes (not neutron stars or quasars): the nodes' `points` (placed ones and the scattered ones
// not yet built into a sector), heaviest first. Nodes of every block size list the biggest, so every node counts.
export const BLACK_HOLE_KIND = {
  id: "marker-black-hole",
  label: "Black holes",
  glyph: "black-hole",
  color: "#c58cff",
  maxBlock: Infinity,
  defaultOn: true,
  source: function (node) {
    return (node.points || []).filter(function (point) { return point.type === "black_hole"; }).map(function (point) {
      return Object.assign({ key: String(point.id) }, point);
    });
  },
  rank: function (item) { return item.map_luminosity_sol || item.mass_solar || item.size || 0; },
  describe: function (item) { return (item.name ? item.name + ", " : "") + "Black Hole" + (item.descriptor ? " (" + item.descriptor + ")" : ""); },
};

// The reasons a star system is worth a mark (MAP.305, MAP.303): the client side of `planetgen.db.interests.REASONS`. A
// system with at least one reason is listed by the nodes' `interests` (blocks of 81 shells and smaller; the server leaves them
// out of bigger ones), each with its `reasons` ({id: {n, score, tier}}). A reason is
//   id         what the server and the layer call it ("habitable", "gas-giant");
//   label      the Menu button's text and the legend line, also the infobox's name for the reason;
//   icon       the symbol the system's mark takes when this is its first reason (poiicons.js);
//   defaultOn  shown until the visitor says otherwise;
//   describe   `function (reason) -> string`: the infobox line for a system's `{n, score, tier}` of it.
// The reasons are in the order the map prefers them: a system with several is drawn once, with the first of them that is on
// (selectMarkers), and its infobox lists them all. A reason is added here (and in the server's registry) and the Menu, the
// legend, the budget and the infobox take it up with no other change.
export const INTEREST_REASONS = [
  {
    id: "habitable",
    label: "Habitable worlds",
    icon: "habitable-world",
    defaultOn: false,
    describe: function (reason) {
      return reason.n + (reason.n === 1 ? " habitable world" : " habitable worlds") + (reason.tier == null ? "" : ", best " + TIER_NAMES[reason.tier]);
    },
  },
  {
    id: "gas-giant",
    label: "Gas giants",
    icon: "gas-giant",
    defaultOn: false,
    describe: function (reason) { return reason.n + (reason.n === 1 ? " gas giant" : " gas giants"); },
  },
];

// A node's system as an item of one reason: `key`, the system's `id`, `name` and place, how many bodies (`n`), the best
// score's tier as an index (0 Blue to 3 Red, `tier`; absent for a reason with no score) and every reason the system has
// (`reasons`, for the infobox).
function interestItems(node, reason) {
  return (node.interests || []).filter(function (system) { return system.reasons && system.reasons[reason.id]; }).map(function (system) {
    const mine = system.reasons[reason.id];
    const item = Object.assign({ key: String(system.id), interest: reason.id }, system, { n: mine.n });
    if (mine.tier != null) item.tier = TIER_NAMES.indexOf(mine.tier);
    if (mine.score != null) item.score = mine.score;
    return item;
  });
}

// The marker kind of an interest reason: its id is "marker-" and the reason's, the Menu layer and the address's `show=`.
export function interestKind(reason) {
  return {
    id: "marker-" + reason.id,
    label: reason.label,
    icon: reason.icon,
    maxBlock: 81,
    defaultOn: reason.defaultOn,
    interest: reason.id,
    source: function (node) { return interestItems(node, reason); },
    rank: function (item) { return item.n || 0; },
    describe: function (item) { return (item.name ? item.name + ", " : "") + reason.describe(item); },
  };
}

export const INTEREST_KINDS = INTEREST_REASONS.map(interestKind);

// The star systems with a habitable world (MAP.220, MAP.303): an overall PHI-6 Yellow or better planet or moon.
export const HABITABLE_KIND = INTEREST_KINDS.find(function (kind) { return kind.interest === "habitable"; });

// The star systems with a gas giant (MAP.305).
export const GAS_GIANT_KIND = INTEREST_KINDS.find(function (kind) { return kind.interest === "gas-giant"; });

// The infobox lines of a system's reasons, in the registry's order: `[{id, label, text}]`.
export function interestLines(system) {
  const lines = [];
  INTEREST_REASONS.forEach(function (reason) {
    const mine = system && system.reasons && system.reasons[reason.id];
    if (!mine) return;
    const shaped = { n: mine.n, score: mine.score, tier: mine.tier == null ? null : TIER_NAMES.indexOf(mine.tier) };
    lines.push({ id: reason.id, label: reason.label, text: reason.describe(shaped) });
  });
  return lines;
}

// The block size (shells) of a node key "m.shell.band.wedge".
export function blockOfKey(key) {
  const block = parseInt(String(key).split(".")[0], 10);
  return isFinite(block) ? block : Infinity;
}

// A registry of marker kinds, in the order they were registered (the order the legend and the budget use).
export function createMarkerKinds() {
  const kinds = [];

  function register(kind) {
    if (!kind || !kind.id || !kind.label || !(GLYPHS[kind.glyph] || (!kind.glyph && poiIcon(kind.icon)))) {
      throw new Error("a marker kind needs an id, a label and a known glyph or icon");
    }
    if (kinds.some(function (known) { return known.id === kind.id; })) throw new Error("marker kind " + kind.id + " is already registered");
    kinds.push(Object.assign({ maxBlock: Infinity, defaultOn: true, extentPc: null, rank: function () { return 0; },
      describe: function (item) { return item.name || kind.label; } }, kind));
    return kinds[kinds.length - 1];
  }

  return {
    register: register,
    kinds: function () { return kinds.slice(); },
    get: function (id) { return kinds.find(function (kind) { return kind.id === id; }) || null; },
  };
}

// Picks what to draw. `options`:
//   kinds      the registered kinds;
//   isOn       `function (kindId) -> boolean`;
//   nodes      `[{key, node}]`, the nodes of the view that are loaded;
//   inView     optional `function (item) -> boolean`: false leaves an item out (off screen, behind the camera);
//   pxOf       optional `function (item, kind) -> pixels` for an item with an `extentPc`: its on-screen size;
//   budget     optional, MARKER_BUDGET.
// Returns `[{kind, item}]`: per kind, highest rank first, at most the kind's share of the budget (an even share,
// what a kind does not use going to the others), a node of a block bigger than a kind's `maxBlock` not counted for it,
// an item listed by several nodes once, and a system with several interest reasons once (under the first reason that is on).
export function selectMarkers(options) {
  const budget = options.budget == null ? MARKER_BUDGET : options.budget;
  const inView = options.inView || function () { return true; };
  const lists = [];
  // A system with several interest reasons is drawn once, by the first reason that is on (MAP.305).
  const claimed = new Set();
  options.kinds.forEach(function (kind) {
    if (!options.isOn(kind.id)) return;
    const seen = new Map();
    options.nodes.forEach(function (entry) {
      if (!entry.node || blockOfKey(entry.key) > kind.maxBlock) return;
      kind.source(entry.node).forEach(function (item) {
        if (!seen.has(item.key)) seen.set(item.key, item);
      });
    });
    const items = Array.from(seen.values()).filter(function (item) {
      if (!inView(item)) return false;
      if (kind.extentPc && options.pxOf && !(options.pxOf(item, kind) >= MIN_EXTENT_PX)) return false;
      if (kind.interest) {
        if (claimed.has(item.id)) return false;
        claimed.add(item.id);
      }
      return true;
    });
    items.sort(function (p, q) { return kind.rank(q) - kind.rank(p) || (p.key < q.key ? -1 : p.key > q.key ? 1 : 0); });
    lists.push({ kind: kind, items: items, taken: 0 });
  });
  // Even shares, then the leftovers handed round, until the budget is spent or every list is.
  let left = budget;
  let open = lists.filter(function (list) { return list.items.length > 0; });
  while (left > 0 && open.length) {
    const share = Math.max(1, Math.floor(left / open.length));
    open.forEach(function (list) {
      const take = Math.min(share, list.items.length - list.taken, left);
      list.taken += take;
      left -= take;
    });
    open = open.filter(function (list) { return list.taken < list.items.length && left > 0; });
  }
  const chosen = [];
  lists.forEach(function (list) {
    list.items.slice(0, list.taken).forEach(function (item) { chosen.push({ kind: list.kind, item: item }); });
  });
  return chosen;
}

// The legend line of a kind for the layer registry (maplayers.js): the kind's symbol where it has one (`icon`, drawn
// by the registry in the colour of `token`), else a swatch of the glyph's color.
export function legendOf(kind) {
  if (kind.icon) return { text: kind.label, icon: kind.icon, token: tokenOf(kind.icon, null) };
  return { text: kind.label, swatch: kind.color, swatchClass: "map-marker-swatch" };
}

// The CSS custom property that colours `item` of `kind`: the symbol's own token, or the habitability tier's for a
// symbol the score colours. `null` for a glyph kind (its `color` is a color, not a token).
export function tokenOfItem(kind, item) {
  return kind.icon ? tokenOf(kind.icon, item && item.tier != null ? item.tier : null) : null;
}
