// html/static/labelbudget.js
//
// The map's label budget (MAP.317): which of the things in view are offered a name, so a sparse view names everything
// it draws and a dense one names a slice of every domain instead of only the brightest stars. Pure; no DOM, no
// three.js (tests/js/labelbudget.test.mjs). The engine (maplabels.js) then places the offered names by rank, none
// overlapping.
//
// The budget N is the screen area over the area one name takes (`labelBudget`). The things offered are the ones the
// map draws (`drawn` not false: a name is never offered for a marker or star that is not on the screen as a mark), in
// five domains that share N by weight (SLICE_WEIGHTS), a domain with less than its share handing the rest to the
// others:
//
//   inhabited    worlds with a population, ranked by population and kept per type (`category`: planet, moon, station...,
//                so a station does not compete with a planet); the habitable worlds with none ranked by habitability
//                score, which is all of the slice until population is stored (MAP.318). A world that is inhabited
//                (`inhabited: true`, whatever its population is yet) is named before one that is not, and a type with an
//                inhabited world in it is served first (MAP.321);
//   uninhabited  uninhabitable worlds, gas giants, asteroid belts (fields included) and the non-luminous star-like
//                objects, an equal share per `category`, each ranked by `mass` (solar masses; a belt by `volume`, the
//                volume of space it occupies in cubic parsecs, standing in for its mass; a body whose mass the map
//                does not carry by its `size`, parsecs across);
//   star         stars, a third of the slice each by size (`radiusSol`), brightness (`luminositySol`) and mass
//                (`massSolar`, else the main-sequence estimate of the brightness); a star in two lists counts once and
//                the next in line fills its place;
//   nebula       nebulae and dust clouds by `volume` (cubic parsecs);
//   region       regions by `shareCharted` (0 to 1), ties by `population`.
//
// An item with `first: true` (a region at the farthest zoom, a sector) goes before everything else; among them a sector goes by the
// total mass of its generated contents, heaviest first (`massKg`, MAP.322).

// Percent of the budget each domain gets while it has things to name; they need not add to 100 (they are shares of
// their sum), and what a domain cannot use goes to the others in proportion.
export const SLICE_WEIGHTS = { inhabited: 30, star: 30, uninhabited: 25, nebula: 10, region: 5 };

// The screen one name takes, in CSS pixels: a short name (about 14 characters at the label font), its padding and
// the line it stands on.
export const LABEL_CELL_PX = { width: 120, height: 24 };

// The main-sequence mass-luminosity exponent (L ~ M^3.5 over the sun's range), used only to rank stars whose mass the
// map does not carry; it is a ranking, never shown.
export const MASS_LUMINOSITY_EXPONENT = 3.5;

// How many names a viewport of `viewport.width` by `viewport.height` CSS pixels can hold.
export function labelBudget(viewport, cell) {
  const size = cell || LABEL_CELL_PX;
  return Math.max(1, Math.floor((viewport.width * viewport.height) / (size.width * size.height)));
}

// A star's mass in solar masses: the one the item carries, else the estimate from its brightness.
export function starMass(item) {
  if (item.massSolar > 0) return item.massSolar;
  return item.luminositySol > 0 ? Math.pow(item.luminositySol, 1 / MASS_LUMINOSITY_EXPONENT) : 0;
}

function byKey(item) {
  return String(item.id);
}

// Sorts a copy of `items` by the `keys` (functions to numbers, higher first), then by id so the order is stable.
function ranked(items, keys) {
  return items.slice().sort(function (a, b) {
    for (let i = 0; i < keys.length; i++) {
      const d = keys[i](b) - keys[i](a);
      if (d) return d;
    }
    const ka = byKey(a);
    const kb = byKey(b);
    return ka < kb ? -1 : ka > kb ? 1 : 0;
  });
}

// Takes one from each list in turn (skipping an empty one, so what one cannot give the others do), once per item.
function roundRobin(lists) {
  const seen = new Set();
  const out = [];
  const at = lists.map(function () { return 0; });
  let open = true;
  while (open) {
    open = false;
    for (let i = 0; i < lists.length; i++) {
      while (at[i] < lists[i].length && seen.has(lists[i][at[i]])) at[i]++;
      if (at[i] >= lists[i].length) continue;
      open = true;
      const item = lists[i][at[i]++];
      seen.add(item);
      out.push(item);
    }
  }
  return out;
}

function groupBy(items, key) {
  const groups = new Map();
  items.forEach(function (item) {
    const k = key(item);
    if (!groups.has(k)) groups.set(k, []);
    groups.get(k).push(item);
  });
  return Array.from(groups.keys()).sort().map(function (k) { return groups.get(k); });
}

const num = function (value) { return value > 0 ? value : 0; };
const isInhabited = function (item) { return item.inhabited === true; };

// Each domain's own order: the thing named first first.
const DOMAINS = {
  inhabited: function (items) {
    // Per type, an equal share; within a type those with a population by population, then by habitability.
    const groups = groupBy(items, function (item) { return item.category || ""; });
    groups.sort(function (a, b) { return (b.some(isInhabited) ? 1 : 0) - (a.some(isInhabited) ? 1 : 0); });
    const lists = groups.map(function (group) {
      return ranked(group, [
        function (item) { return isInhabited(item) ? 1 : 0; },
        function (item) { return num(item.population) > 0 ? 1 : 0; },
        function (item) { return num(item.population); },
        function (item) { return num(item.habitability); },
      ]);
    });
    return roundRobin(lists);
  },
  uninhabited: function (items) {
    const lists = groupBy(items, function (item) { return item.category || ""; }).map(function (group) {
      return ranked(group, [function (item) { return num(item.mass) || num(item.volume) || num(item.size); }]);
    });
    return roundRobin(lists);
  },
  star: function (items) {
    return roundRobin([
      ranked(items, [function (item) { return num(item.radiusSol); }]),
      ranked(items, [function (item) { return num(item.luminositySol); }]),
      ranked(items, [starMass]),
    ]);
  },
  nebula: function (items) {
    return ranked(items, [function (item) { return num(item.volume); }]);
  },
  region: function (items) {
    return ranked(items, [function (item) { return num(item.shareCharted); }, function (item) { return num(item.population); }]);
  },
};

// Chooses the names to offer. `items` is [{id, domain, category?, drawn?, first?, ...the domain's numbers}]; options:
// `budget` (the count to offer, default labelBudget is the caller's to work out), `weights` (default SLICE_WEIGHTS).
// Returns the chosen items, the first to name first, each as a copy with `rank` (0 is first).
export function chooseLabels(items, options) {
  const weights = (options && options.weights) || SLICE_WEIGHTS;
  const budget = options && options.budget !== undefined ? options.budget : Infinity;
  const drawn = items.filter(function (item) { return item.drawn !== false && weights[item.domain] > 0 && DOMAINS[item.domain]; });
  // The items that go first keep their order, except that those with a `massKg` (a sector's total mass, kilograms, MAP.322) go heaviest
  // first, after any that have none: a crowded region names its heaviest sectors and a short budget drops the lightest.
  const chosen = drawn.filter(function (item) { return item.first; }).map(function (item, n) { return { item: item, n: n }; })
    .sort(function (a, b) {
      const ka = a.item.massKg === undefined ? Infinity : num(a.item.massKg);
      const kb = b.item.massKg === undefined ? Infinity : num(b.item.massKg);
      return ka !== kb ? (kb > ka ? 1 : -1) : a.n - b.n;
    }).map(function (entry) { return entry.item; });
  const rest = drawn.filter(function (item) { return !item.first; });
  const lists = Object.keys(DOMAINS).map(function (domain) {
    return { domain: domain, weight: weights[domain] || 0, items: DOMAINS[domain](rest.filter(function (item) { return item.domain === domain; })), taken: 0 };
  }).filter(function (list) { return list.items.length > 0 && list.weight > 0; });
  // Highest averages (D'Hondt): the domain with the most weight per name already taken goes next, so the names come
  // out in step with the weights and a domain that runs dry stops claiming.
  const out = chosen.slice();
  while (out.length < budget) {
    let best = null;
    lists.forEach(function (list) {
      if (list.taken >= list.items.length) return;
      if (!best || list.weight / (list.taken + 1) > best.weight / (best.taken + 1)) best = list;
    });
    if (!best) break;
    out.push(best.items[best.taken++]);
  }
  return out.slice(0, budget).map(function (item, rank) { return Object.assign({}, item, { rank: rank }); });
}
