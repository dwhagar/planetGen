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
//   glyph       the shape drawn (a key of GLYPHS);
//   color       the glyph's color (CSS), also the legend swatch;
//   minLevel    the finest-than-or-equal tile level that lists this kind completely (tiles coarser than this do not
//               list it, so the kind is not drawn from them);
//   defaultOn   shown until the visitor says otherwise (default true);
//   source      `function (tile) -> [item]`: the kind's items in one tile. An item has `key` (unique within the
//               kind), `x`, `y`, `z` (galaxy frame, parsecs) and anything the kind needs to describe it;
//   rank        `function (item) -> number`: higher is kept first when the budget cuts the list;
//   extentPc    optional `function (item) -> parsecs`: the item's own size, for a kind whose items are extended
//               (a nebula): it is hidden while that is under MIN_EXTENT_PX on screen;
//   describe    `function (item) -> string`: the hover text.
//
// What shows is chosen by `selectMarkers`: only kinds that are on, only items in view, and no more than the
// budget (MARKER_BUDGET) in all, the kinds sharing it evenly and handing on what one cannot use.

// The most markers drawn at once (MAP.116's rule for the stars: a crowded view must stay readable). Each kind gets
// an even share of it; a kind that has less hands the rest on.
export const MARKER_BUDGET = 240;

// A marker of a kind with its own size (`extentPc`) is hidden while it would be under this many pixels across.
export const MIN_EXTENT_PX = 2;

// Pixels across a drawn glyph.
export const MARKER_PX = 18;

// The shapes a kind draws, each a function `(ctx, size, color, outline)` that draws into a canvas context `size`
// pixels square: the one place a shape is defined, so a symbol library (UX.101) can replace it by name.
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
  // Overall habitability (UX.101): a solid circle inside an upright triangle inside a ring.
  "habitable": function (ctx, size, color, outline) {
    const c = size / 2;
    const ring = c * 0.88;
    ctx.lineJoin = "round";
    const triangle = function () {
      ctx.beginPath();
      for (let i = 0; i < 3; i++) {
        const angle = -Math.PI / 2 + (i * 2 * Math.PI) / 3;
        const x = c + ring * 0.78 * Math.cos(angle);
        const y = c + ring * 0.78 * Math.sin(angle) + size * 0.04;
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.closePath();
    };
    [[outline, 0.2], [color, 0.09]].forEach(function (pass) {
      ctx.strokeStyle = pass[0];
      ctx.lineWidth = size * pass[1];
      ctx.beginPath();
      ctx.arc(c, c, ring, 0, Math.PI * 2);
      ctx.stroke();
      triangle();
      ctx.stroke();
    });
    ctx.beginPath();
    ctx.arc(c, c + size * 0.08, size * 0.15, 0, Math.PI * 2);
    ctx.fillStyle = outline;
    ctx.fill();
    ctx.lineWidth = size * 0.04;
    ctx.strokeStyle = outline;
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(c, c + size * 0.08, size * 0.11, 0, Math.PI * 2);
    ctx.fillStyle = color;
    ctx.fill();
  },
};

// The markers of black holes (not neutron stars or quasars): the tiles' `points` (placed ones and the scattered ones
// not yet built into a sector), heaviest first. Tiles of every level list the biggest, so it is complete from level 0.
export const BLACK_HOLE_KIND = {
  id: "marker-black-hole",
  label: "Black holes",
  glyph: "black-hole",
  color: "#c58cff",
  minLevel: 0,
  defaultOn: true,
  source: function (tile) {
    return (tile.points || []).filter(function (point) { return point.type === "black_hole"; }).map(function (point) {
      return Object.assign({ key: String(point.id) }, point);
    });
  },
  rank: function (item) { return item.map_luminosity_sol || item.mass_solar || item.size || 0; },
  describe: function (item) { return (item.name ? item.name + ", " : "") + "Black Hole" + (item.descriptor ? " (" + item.descriptor + ")" : ""); },
};

// The star systems with a habitable world (MAP.220): an overall PHI-6 Blue or Green planet or moon. The tiles' `worlds`
// list them from the finest levels (the server leaves it out of coarser tiles), most worlds first.
export const HABITABLE_KIND = {
  id: "marker-habitable",
  label: "Habitable worlds",
  glyph: "habitable",
  color: "#4fd18b",
  minLevel: 0,
  defaultOn: false,
  source: function (tile) {
    return (tile.worlds || []).map(function (system) { return Object.assign({ key: String(system.id) }, system); });
  },
  rank: function (item) { return item.worlds || 0; },
  describe: function (item) {
    return (item.name ? item.name + ", " : "") + item.worlds + (item.worlds === 1 ? " habitable world" : " habitable worlds");
  },
};

// The level of a tile key "level/ix/iy/iz".
export function levelOfKey(key) {
  const level = parseInt(String(key).split("/")[0], 10);
  return isFinite(level) ? level : 0;
}

// A registry of marker kinds, in the order they were registered (the order the legend and the budget use).
export function createMarkerKinds() {
  const kinds = [];

  function register(kind) {
    if (!kind || !kind.id || !kind.label || !GLYPHS[kind.glyph]) throw new Error("a marker kind needs an id, a label and a known glyph");
    if (kinds.some(function (known) { return known.id === kind.id; })) throw new Error("marker kind " + kind.id + " is already registered");
    kinds.push(Object.assign({ minLevel: 0, defaultOn: true, extentPc: null, rank: function () { return 0; },
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
//   tiles      `[{key, tile}]`, the tiles of the view that are loaded;
//   inView     optional `function (item) -> boolean`: false leaves an item out (off screen, behind the camera);
//   pxOf       optional `function (item, kind) -> pixels` for an item with an `extentPc`: its on-screen size;
//   budget     optional, MARKER_BUDGET.
// Returns `[{kind, item}]`: per kind, highest rank first, at most the kind's share of the budget (an even share,
// what a kind does not use going to the others), a tile coarser than a kind's `minLevel` not counted for it, an item
// listed by several tiles once.
export function selectMarkers(options) {
  const budget = options.budget == null ? MARKER_BUDGET : options.budget;
  const inView = options.inView || function () { return true; };
  const lists = [];
  options.kinds.forEach(function (kind) {
    if (!options.isOn(kind.id)) return;
    const seen = new Map();
    options.tiles.forEach(function (entry) {
      if (!entry.tile || levelOfKey(entry.key) < kind.minLevel) return;
      kind.source(entry.tile).forEach(function (item) {
        if (!seen.has(item.key)) seen.set(item.key, item);
      });
    });
    const items = Array.from(seen.values()).filter(function (item) {
      if (!inView(item)) return false;
      if (kind.extentPc && options.pxOf && !(options.pxOf(item, kind) >= MIN_EXTENT_PX)) return false;
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

// The legend line of a kind for the layer registry (maplayers.js): its swatch is the glyph's color.
export function legendOf(kind) {
  return { text: kind.label, swatch: kind.color, swatchClass: "map-marker-swatch" };
}
