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
//   defaultOn   shown until the visitor says otherwise (default true; the inhabited, habitable and gas giant reasons are on, MAP.303);
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

// The marks for worlds with life (MAP.311): a leaf for any life, a head and shoulders for a sentient species, a tower for a civilization. Each is
// a disc ringed in the outline colour with a dark sign on it, so the shape and not only the colour tells them apart.
function discWith(draw) {
  return function (ctx, size, color, outline) {
    const c = size / 2;
    ctx.beginPath();
    ctx.arc(c, c, c * 0.92, 0, Math.PI * 2);
    ctx.fillStyle = outline;
    ctx.fill();
    ctx.beginPath();
    ctx.arc(c, c, c * 0.76, 0, Math.PI * 2);
    ctx.fillStyle = color;
    ctx.fill();
    ctx.fillStyle = "#14161f";
    ctx.strokeStyle = "#14161f";
    ctx.lineJoin = "round";
    draw(ctx, size);
  };
}

GLYPHS.life = discWith(function (ctx, size) {
  // A leaf: two arcs meeting at the tip, with a stalk.
  const u = size / 24;
  ctx.beginPath();
  ctx.moveTo(6.5 * u, 17.5 * u);
  ctx.quadraticCurveTo(5 * u, 7 * u, 17.5 * u, 6 * u);
  ctx.quadraticCurveTo(18.5 * u, 16 * u, 6.5 * u, 17.5 * u);
  ctx.fill();
  ctx.lineWidth = 1.4 * u;
  ctx.beginPath();
  ctx.moveTo(5 * u, 19 * u);
  ctx.lineTo(11 * u, 13 * u);
  ctx.stroke();
});

GLYPHS.sentient = discWith(function (ctx, size) {
  // A head and shoulders.
  const u = size / 24;
  ctx.beginPath();
  ctx.arc(12 * u, 9.2 * u, 3.1 * u, 0, Math.PI * 2);
  ctx.fill();
  ctx.beginPath();
  ctx.moveTo(5.8 * u, 18.5 * u);
  ctx.quadraticCurveTo(12 * u, 7.5 * u, 18.2 * u, 18.5 * u);
  ctx.closePath();
  ctx.fill();
});

GLYPHS.civilization = discWith(function (ctx, size) {
  // A skyline of three towers.
  const u = size / 24;
  ctx.fillRect(5.2 * u, 12 * u, 3.8 * u, 7 * u);
  ctx.fillRect(10.1 * u, 6.2 * u, 3.8 * u, 12.8 * u);
  ctx.fillRect(15 * u, 10 * u, 3.8 * u, 9 * u);
});

// The planet classes the map can mark (MAP.220): `planetgen.tuning.PLANET_CLASSES`' letters; tests/test_system_interests.py holds the
// two lists to each other. A class's reason is "class-" and its letter, the same id the server's `system_interests` rows use.
export const PLANET_CLASS_LETTERS = "ABCDEFGHIJKLMNOPSTV".split("");

// The life classes the map can mark (MAP.351): `planetgen.physics.life_classes.CLASSES`' ids; tests/test_system_interests.py holds the two lists to
// each other. A life class's reason is "lifeclass-" and its id, the same id the server's `system_interests` rows use.
export const LIFE_CLASS_IDS = ["1a", "1b", "1c", "2", "3a", "3b", "4a", "4b", "4c", "5a", "5b", "6", "7", "8", "9"];

// The biggest block (shells) whose node lists the class reasons (`planetgen.db.nodes.CLASS_ABOVE`): a bigger block has too many systems.
export const CLASS_MAX_BLOCK = 3;

// The colour of every class glyph (the letter tells them apart) and of the letter on it.
const CLASS_COLOR = "#e3a63b";
const LIFE_CLASS_COLOR = "#5fc77a";
const CLASS_LETTER_COLOR = "#14161f";

// A class's glyph: a disc ringed in the outline colour with the class letter on it.
function classGlyph(letter) {   // `letter` is the text on the disc: a planet class letter or a life class id
  return function (ctx, size, color, outline) {
    const c = size / 2;
    ctx.beginPath();
    ctx.arc(c, c, c * 0.92, 0, Math.PI * 2);
    ctx.fillStyle = outline;
    ctx.fill();
    ctx.beginPath();
    ctx.arc(c, c, c * 0.76, 0, Math.PI * 2);
    ctx.fillStyle = color;
    ctx.fill();
    ctx.fillStyle = CLASS_LETTER_COLOR;
    ctx.font = "bold " + Math.round(size * (letter.length > 1 ? 0.5 : 0.64)) + "px sans-serif";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(letter, c, c + size * 0.04);
  };
}

PLANET_CLASS_LETTERS.forEach(function (letter) { GLYPHS["class-" + letter] = classGlyph(letter); });
LIFE_CLASS_IDS.forEach(function (id) { GLYPHS["lifeclass-" + id] = classGlyph(id); });

// MAP.339: black holes are chosen for the view, not for the whole list: every node lists its own, whatever its distance, so the
// heaviest 240 were mostly far from what the camera looks at (237 of 240 at a 10 pc zoom). The rule, with `d` the distance from the
// camera to the point it looks at (the view centre) and `r` a hole's distance from that point (both in parsecs, so it is the same at
// every zoom):
//   reach   a hole within NEAR_FIELD_REACH view distances (r <= 3 d) is a candidate; one beyond only if its mass is at least
//           BLACK_HOLE_MIN_MASS * (r / (3 d))^2, so a rare heavy hole stays and a crowd of light ones does not;
//   rank    mass over r squared (the nearest to what is looked at, the heaviest among those);
//   cap     as many as fit the screen at one mark per `densityPx` by `densityPx` pixels (a user preference, `map.blackHoleDensity`,
//           BLACK_HOLE_DENSITY_PX by default, mappreferences.js), a heavy hole following the same rule as a light one.
export const NEAR_FIELD_REACH = 3;
// The lightest stellar black hole, in solar masses: the floor of the reach rule.
export const BLACK_HOLE_MIN_MASS = 5;
// The side in pixels of the square one black hole mark is allowed on a screen, by default.
export const BLACK_HOLE_DENSITY_PX = 120;

// A black hole's mass in solar masses, or null when the item does not say (a hole that does not give its mass passes the reach
// rule only inside the reach).
export function blackHoleMass(item) {
  const mass = [item.mass_solar, item.system_mass_solar].find(function (value) { return Number.isFinite(value) && value > 0; });
  return mass === undefined ? null : mass;
}

// The rule for a view: `view` is `{centre: [x, y, z], distance, width, height, densityPx}` (the point the camera looks at, its distance
// from the camera, the canvas size in CSS pixels and the density preference's side). Returns `{accepts(item), rank(item), cap}`
// for `selectMarkers`' `near`.
export function nearFieldRule(view) {
  const reach = NEAR_FIELD_REACH * Math.max(view.distance, 1e-9);
  const side = view.densityPx > 0 ? view.densityPx : BLACK_HOLE_DENSITY_PX;
  const distanceOf = function (item) {
    return Math.hypot(item.x - view.centre[0], item.y - view.centre[1], item.z - view.centre[2]);
  };
  return {
    accepts: function (item) {
      const r = distanceOf(item);
      if (r <= reach) return true;
      const mass = blackHoleMass(item);
      return mass !== null && mass >= BLACK_HOLE_MIN_MASS * (r / reach) * (r / reach);
    },
    rank: function (item) {
      const mass = blackHoleMass(item) || BLACK_HOLE_MIN_MASS;
      const r = Math.max(distanceOf(item), 0.01 * view.distance, 1e-9);
      return mass / (r * r);
    },
    cap: Math.max(1, Math.floor((Math.max(view.width, 1) * Math.max(view.height, 1)) / (side * side))),
  };
}

// The markers of black holes (not neutron stars or quasars): the nodes' `points` (placed ones and the scattered ones
// not yet built into a sector), heaviest first. Nodes of every block size list the biggest, so every node counts.
export const BLACK_HOLE_KIND = {
  id: "marker-black-hole",
  label: "Black holes",
  glyph: "black-hole",
  color: "#c58cff",
  maxBlock: Infinity,
  defaultOn: true,
  crowds: true,       // MAP.340: marks too close on the screen merge into one with a count (`mergeMarks`), among black holes only
  nearField: true,    // MAP.339: chosen by `nearFieldRule` (the view's near field, then by mass over distance squared, then by screen area)
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
//   body       what one body of this reason is called (a world on the map's system level, MAP.321: "habitable world");
//   describe   `function (reason) -> string`: the infobox line for a system's `{n, score, tier}` of it.
// The reasons are in the order the map prefers them (MAP.321: the best world's symbol, an inhabited one before a habitable one before a barren one before a
// gas giant before a belt): a system with several is drawn once, with the first of them that is on (selectMarkers), and its infobox
// lists them all. A reason is added here (and in the server's registry) and the Menu, the
// legend, the budget and the infobox take it up with no other change.
export const INTEREST_REASONS = [
  {
    id: "inhabited",
    label: "Inhabited worlds",
    body: "inhabited world",
    icon: "inhabited-world",
    defaultOn: true,
    describe: function (reason) {
      return reason.n + (reason.n === 1 ? " inhabited world" : " inhabited worlds") + (reason.tier == null ? "" : ", best " + TIER_NAMES[reason.tier]);
    },
  },
  {
    id: "civilization",
    label: "Civilizations",
    body: "civilization",
    glyph: "civilization",
    color: "#f08a4b",
    defaultOn: false,
    describe: function (reason) {
      return reason.n + (reason.n === 1 ? " world with a civilization" : " worlds with a civilization") +
        (reason.band ? ", most advanced " + reason.band : "");
    },
  },
  {
    id: "sentient",
    label: "Sentient life",
    body: "world with sentient life",
    glyph: "sentient",
    color: "#e0c24b",
    defaultOn: false,
    describe: function (reason) { return reason.n + (reason.n === 1 ? " world with sentient life" : " worlds with sentient life"); },
  },
  {
    id: "life",
    label: "Life",
    body: "world with life",
    glyph: "life",
    color: "#7bd16f",
    defaultOn: false,
    describe: function (reason) { return reason.n + (reason.n === 1 ? " world with life" : " worlds with life"); },
  },
  {
    id: "habitable",
    label: "Habitable worlds",
    body: "habitable world",
    icon: "habitable-world",
    defaultOn: true,
    describe: function (reason) {
      return reason.n + (reason.n === 1 ? " habitable world" : " habitable worlds") + (reason.tier == null ? "" : ", best " + TIER_NAMES[reason.tier]);
    },
  },
  {
    id: "uninhabitable",
    label: "Uninhabitable worlds",
    body: "uninhabitable world",
    icon: "uninhabitable-world",
    defaultOn: false,
    describe: function (reason) { return reason.n + (reason.n === 1 ? " uninhabitable world" : " uninhabitable worlds"); },
  },
  {
    id: "gas-giant",
    label: "Gas giants",
    body: "gas giant",
    icon: "gas-giant",
    defaultOn: true,
    describe: function (reason) { return reason.n + (reason.n === 1 ? " gas giant" : " gas giants"); },
  },
  {
    id: "asteroid-field",
    label: "Asteroid fields",
    body: "asteroid belt",
    icon: "asteroid-field",
    defaultOn: false,
    describe: function (reason) { return reason.n + (reason.n === 1 ? " asteroid belt" : " asteroid belts"); },
  },
].concat(PLANET_CLASS_LETTERS.map(function (letter) {
  return {
    id: "class-" + letter,
    label: "Class " + letter + " worlds",
    body: "class " + letter + " world",
    glyph: "class-" + letter,
    color: CLASS_COLOR,
    maxBlock: CLASS_MAX_BLOCK,
    family: "class",
    defaultOn: false,
    describe: function (reason) { return reason.n + (reason.n === 1 ? " class " + letter + " world" : " class " + letter + " worlds"); },
  };
})).concat(LIFE_CLASS_IDS.map(function (id) {
  return {
    id: "lifeclass-" + id,
    label: "Class " + id + " life",
    body: "world with class " + id + " life",
    glyph: "lifeclass-" + id,
    color: LIFE_CLASS_COLOR,
    maxBlock: CLASS_MAX_BLOCK,
    family: "lifeclass",
    defaultOn: false,
    describe: function (reason) {
      return reason.n + (reason.n === 1 ? " world with class " + id + " life" : " worlds with class " + id + " life");
    },
  };
}));

// A node's system as an item of one reason: `key`, the system's `id`, `name` and place, how many bodies (`n`), the best
// score's tier as an index (0 Blue to 3 Red, `tier`; absent for a reason with no score) and every reason the system has
// (`reasons`, for the infobox).
function interestItems(node, reason) {
  return (node.interests || []).filter(function (system) { return system.reasons && system.reasons[reason.id]; }).map(function (system) {
    const mine = system.reasons[reason.id];
    const item = Object.assign({ key: String(system.id), interest: reason.id }, system, { n: mine.n });
    if (mine.tier != null) item.tier = TIER_NAMES.indexOf(mine.tier);
    if (mine.score != null) item.score = mine.score;
    if (mine.band != null) {
      item.tech = mine.tech;
      item.band = mine.band;
    }
    return item;
  });
}

// The marker kind of an interest reason: its id is "marker-" and the reason's, the Menu layer and the address's `show=`.
export function interestKind(reason) {
  const kind = {
    id: "marker-" + reason.id,
    label: reason.label,
    maxBlock: reason.maxBlock || 81,
    defaultOn: reason.defaultOn,
    interest: reason.id,
    family: reason.family || null,
    bodyName: reason.body,
    source: function (node) { return interestItems(node, reason); },
    rank: function (item) { return item.n || 0; },
    describe: function (item) { return (item.name ? item.name + ", " : "") + reason.describe(item); },
  };
  if (reason.glyph) {
    kind.glyph = reason.glyph;
    kind.color = reason.color;
  } else {
    kind.icon = reason.icon;
  }
  return kind;
}

export const INTEREST_KINDS = INTEREST_REASONS.map(interestKind);

// The technology bands of a civilization (MAP.344): `planetgen.tuning.TECH_BANDS`' names, lowest first; tests/test_system_interests.py holds the two
// lists to each other. A system's band is that of its most advanced species (the civilization reason's `band`, from the species' `tech_level`).
export const TECH_BAND_NAMES = ["Primal", "Agrarian", "Industrial", "Atomic", "Fusion", "Gravitic", "Antimatter", "Transcendent"];

// The colour of each band's civilization glyph, dull for the lowest rising to bright, so a band can be told at a glance as well as chosen.
const TECH_BAND_COLORS = ["#8f7a5c", "#b5903f", "#d9822b", "#e8602c", "#e04a6b", "#b45bd6", "#6f7be8", "#4fc3e8"];

// One optional marker for each technology band (MAP.344): the systems whose most advanced species is in that band, drawn as the civilization
// tower in the band's colour. They share the civilization reason (`interest`), so a system is drawn once (by its band's kind when that is on,
// else by the Civilizations kind) and the card lists the reason once; a system whose species is not drawn yet has no band and is in the
// Civilizations kind alone.
const CIVILIZATION_REASON = INTEREST_REASONS.find(function (reason) { return reason.id === "civilization"; });
export const TECH_KINDS = TECH_BAND_NAMES.map(function (band, index) {
  const kind = interestKind(CIVILIZATION_REASON);
  kind.id = "marker-tech-" + band.toLowerCase();
  kind.label = band + " civilizations";
  kind.bodyName = band.toLowerCase() + " civilization";
  kind.color = TECH_BAND_COLORS[index];
  kind.family = "tech";
  kind.defaultOn = false;
  kind.source = function (node) { return interestItems(node, CIVILIZATION_REASON).filter(function (item) { return item.band === band; }); };
  return kind;
});

// The kinds the Galaxy Map registers: the interest kinds with the technology bands ahead of the Civilizations kind, so a system with a band is
// drawn by its band's colour when that is on. (The open system's bodies use `INTEREST_KINDS`, which have no bands.)
const CIVILIZATION_AT = INTEREST_KINDS.findIndex(function (kind) { return kind.interest === "civilization"; });
export const MAP_INTEREST_KINDS = INTEREST_KINDS.slice(0, CIVILIZATION_AT).concat(TECH_KINDS, INTEREST_KINDS.slice(CIVILIZATION_AT));

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
    const shaped = { n: mine.n, score: mine.score, tier: mine.tier == null ? null : TIER_NAMES.indexOf(mine.tier), tech: mine.tech, band: mine.band };
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
//   isPulled   optional `function (kindId) -> boolean`: a kind pulled forward (MAP.306) has first call on the budget and on a
//              system's symbol, and its entries carry `pulled: true`;
//   near       optional, a view's `nearFieldRule` (MAP.339): a kind with `nearField` keeps only the items it accepts, ranks them by
//              its rank and takes at most its cap, before the kinds share the budget;
//   budget     optional, MARKER_BUDGET.
// Returns `[{kind, item}]`: per kind, highest rank first, at most the kind's share of the budget (an even share,
// what a kind does not use going to the others), a node of a block bigger than a kind's `maxBlock` not counted for it,
// an item listed by several nodes once, and a system with several interest reasons once (under the first reason that is on).
export function selectMarkers(options) {
  const budget = options.budget == null ? MARKER_BUDGET : options.budget;
  const inView = options.inView || function () { return true; };
  const lists = [];
  // A system with several interest reasons is drawn once, by the first reason that is on (MAP.305), a reason pulled forward before
  // the others.
  const claimed = new Set();
  const isPulled = options.isPulled || function () { return false; };
  const ordered = options.kinds.filter(function (kind) { return isPulled(kind.id); })
    .concat(options.kinds.filter(function (kind) { return !isPulled(kind.id); }));
  ordered.forEach(function (kind) {
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
    const near = kind.nearField ? options.near : null;
    const wanted = near ? items.filter(near.accepts) : items;
    const rankOf = near ? near.rank : kind.rank;
    wanted.sort(function (p, q) { return rankOf(q) - rankOf(p) || (p.key < q.key ? -1 : p.key > q.key ? 1 : 0); });
    lists.push({ kind: kind, items: near ? wanted.slice(0, near.cap) : wanted, taken: 0 });
  });
  // A kind pulled forward is served first, whole, before the rest share what is left.
  let left = budget;
  lists.forEach(function (list) {
    if (!isPulled(list.kind.id)) return;
    list.taken = Math.min(list.items.length, left);
    left -= list.taken;
  });
  // Even shares, then the leftovers handed round, until the budget is spent or every list is.
  let open = lists.filter(function (list) { return !isPulled(list.kind.id) && list.items.length > 0; });
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
    const pulled = isPulled(list.kind.id);
    list.items.slice(0, list.taken).forEach(function (item) {
      const entry = { kind: list.kind, item: item };
      if (pulled) entry.pulled = true;
      chosen.push(entry);
    });
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

// Marks closer together on the screen than this, in CSS pixels between their centres, are drawn as one (MAP.305): the symbols are
// MARKER_PX across, so these would lie over each other.
export const MERGE_PX = 22;

// Black holes merge only with ones about as far from the camera (MAP.340): within this ratio of the crowd's seed (and of its inverse), so a
// crowd is a crowd in space and not a near hole in front of a far one on the same line of sight, and the Zoom in on them button has
// somewhere to go. `screenOf` gives the camera distance as `depth`; a mark without one merges at any depth.
export const MERGE_DEPTH_RATIO = 1.5;

// Merges the crowded star-level interest marks of `chosen` (`selectMarkers`' `[{kind, item}]`, with any other entries after them) into one
// mark each. `screenOf(item)` is the item's place on the screen in CSS pixels, `{x, y}`, or null for one that has none. A system's
// own interest mark merges (a kind with `interest`, no `bodyRef` or `regionKey`), and so do black holes (a kind with `crowds`, MAP.340),
// each only with its own sort (a black hole also only with ones at about its depth, `MERGE_DEPTH_RATIO`): a crowd of black holes is one mark with their count, the heaviest first among them, `{key, id, x, y, z (the middle
// of its holes), count, holes: true, members, merged: true}`; a body's symbol or a region's token stays as it is. The marks are taken pulled-forward ones first (MAP.306), then in the registry's order of reasons, most bodies first, each one gathering the unmerged marks within
// `radiusPx` of it, so the same crowd always merges the same way and a mark's symbol is its crowd's best reason (an inhabited world
// before a habitable one before the rest). A crowd of two or more becomes one entry whose item is
//   `{key, id, x, y, z (the middle of its systems), n (their bodies), count (their number), tier (the best), members (their items),
//     merged: true, interest}`
// and whose kind is the best member's with the hover text of a crowd. Zooming in spreads the marks and the crowd splits.
export function mergeMarks(chosen, screenOf, radiusPx) {
  const radius = radiusPx === undefined ? MERGE_PX : radiusPx;
  const order = INTEREST_REASONS.map(function (reason) { return reason.id; });
  const placed = [];
  const rest = [];
  chosen.forEach(function (entry) {
    const item = entry.item;
    const mergeable = entry.kind.crowds || (entry.kind.interest && !item.bodyRef && !item.regionKey && !item.follow);
    const spot = mergeable ? screenOf(item) : null;
    if (spot) placed.push({ entry: entry, spot: spot });
    else rest.push(entry);
  });
  placed.sort(function (p, q) {
    const holes = (p.entry.kind.crowds ? 1 : 0) - (q.entry.kind.crowds ? 1 : 0);
    if (holes) return holes;
    if (p.entry.kind.crowds) {
      return (blackHoleMass(q.entry.item) || 0) - (blackHoleMass(p.entry.item) || 0)
        || (p.entry.item.key < q.entry.item.key ? -1 : p.entry.item.key > q.entry.item.key ? 1 : 0);
    }
    return (q.entry.pulled ? 1 : 0) - (p.entry.pulled ? 1 : 0)
      || order.indexOf(p.entry.kind.interest) - order.indexOf(q.entry.kind.interest) || (q.entry.item.n || 0) - (p.entry.item.n || 0)
      || (p.entry.item.key < q.entry.item.key ? -1 : p.entry.item.key > q.entry.item.key ? 1 : 0);
  });
  const taken = new Set();
  const out = [];
  placed.forEach(function (seed) {
    if (taken.has(seed)) return;
    taken.add(seed);
    const crowd = [seed];
    placed.forEach(function (other) {
      if (taken.has(other) || !!other.entry.kind.crowds !== !!seed.entry.kind.crowds) return;
      if (seed.entry.kind.crowds && seed.spot.depth > 0 && other.spot.depth > 0) {
        const ratio = other.spot.depth / seed.spot.depth;
        if (ratio > MERGE_DEPTH_RATIO || ratio < 1 / MERGE_DEPTH_RATIO) return;
      }
      if (Math.hypot(other.spot.x - seed.spot.x, other.spot.y - seed.spot.y) <= radius) {
        taken.add(other);
        crowd.push(other);
      }
    });
    out.push(crowd.length === 1 ? seed.entry : mergedEntry(crowd.map(function (member) { return member.entry; })));
  });
  return out.concat(rest);
}

// The one entry a crowd of marks is drawn as; `entries` are in the order of preference, the first one the best.
function mergedEntry(entries) {
  const best = entries[0];
  if (best.kind.crowds) {
    const mean = function (axis) { return entries.reduce(function (sum, entry) { return sum + entry.item[axis]; }, 0) / entries.length; };
    return {
      kind: crowdKind(best.kind),
      item: { key: "merge:" + best.item.key, id: "merge:" + best.item.id, name: null, x: mean("x"), y: mean("y"), z: mean("z"),
        count: entries.length, holes: true, merged: true, members: entries.map(function (entry) { return entry.item; }) },
    };
  }
  const total = function (key) { return entries.reduce(function (sum, entry) { return sum + (entry.item[key] || 0); }, 0); };
  const item = {
    key: "merge:" + best.item.key, id: "merge:" + best.item.id, name: null,
    x: total("x") / entries.length, y: total("y") / entries.length, z: total("z") / entries.length,
    n: total("n"), count: entries.length, interest: best.kind.interest, merged: true,
    members: entries.map(function (entry) { return entry.item; }),
  };
  if (best.item.tier != null) item.tier = best.item.tier;
  const merged = { kind: crowdKind(best.kind), item: item };
  if (entries.some(function (entry) { return entry.pulled; })) merged.pulled = true;    // a crowd takes the strongest state of its members
  return merged;
}

// The kind of a crowd's mark: its best reason's kind with the hover text of a crowd ("12 points of interest").
const crowdKinds = new Map();
function crowdKind(kind) {
  if (!crowdKinds.has(kind.id)) {
    crowdKinds.set(kind.id, Object.assign({}, kind, {
      describe: function (item) { return item.count + (item.holes ? " black holes" : " points of interest"); },
    }));
  }
  return crowdKinds.get(kind.id);
}
