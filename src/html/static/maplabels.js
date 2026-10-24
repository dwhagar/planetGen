// html/static/maplabels.js
//
// The map label engine (MAP.229, built as part of MAP.268): names placed by priority, with no overlaps, shown and
// hidden by zoom. It is the one label code the Galaxy Map uses (and, later, the Sector Map and the System Map).
//
// A label candidate is a plain object
//   id        unique and stable (the same thing keeps the same id from frame to frame);
//   kind      the layer it belongs to ("sector", "system", "nebula"): the caller only passes the kinds that are on;
//   text      what to write;
//   x, y, z   the anchor in the map's frame;
//   priority  a number, higher placed first (apparent brightness, size, nearness: the caller's rank);
//   pixels    optional: how large the thing is on the screen now; with `minPixels` and `maxPixels` it limits the
//             zoom range the label shows at (a nebula's name only where the nebula is large enough to read).
//
// `placeLabels` is the pure part: it projects the anchors (`project(x, y, z)` gives {x, y, depth} in CSS pixels
// of the view, or null when the point is behind the camera), drops those outside the view or outside their
// size range, and places the rest by priority on a screen grid with a padding so no two overlap, up to a budget.
// A label that was shown last time keeps its place against an equal one (hysteresis: the size range also widens
// by HOLD_SHARE for a label already shown, so it does not flicker at a range edge). Tests:
// tests/js/maplabels.test.mjs.
//
// `createLabelLayer` is the DOM part: one `<span class="map-label">` per placed label in a container laid over
// the map, pooled, faded in and out (no fade under reduced motion) and never picking (the container has
// pointer-events none: a label does not hide the thing's own pick and hover). The text is set with textContent
// only, never markup: names are database content.

export const DEFAULT_BUDGET = 60;
export const PADDING_PX = 4;
export const HOLD_SHARE = 0.15;
const GRID_PX = 96;
const LABEL_HEIGHT_PX = 16;
const GAP_PX = 6;
const FADE_MS = 200;

// The width a text takes at the label font: measured once per text with a canvas where there is one, else by
// the number of letters.
export function createMeasure(font) {
  const cache = new Map();
  let context = null;
  if (typeof document !== "undefined" && document.createElement) {
    try {
      const canvas = document.createElement("canvas");
      context = canvas.getContext && canvas.getContext("2d");
      if (context) context.font = font || "12px sans-serif";
    } catch (e) {
      context = null;
    }
  }
  return function measure(text) {
    let width = cache.get(text);
    if (width === undefined) {
      width = context && context.measureText ? context.measureText(text).width : text.length * 6.5;
      if (!(width > 0)) width = text.length * 6.5;
      cache.set(text, width);
      if (cache.size > 2000) cache.delete(cache.keys().next().value);
    }
    return width;
  };
}

// Whether the candidate's size range lets it show; a label shown last time keeps showing a little beyond the range.
function inRange(candidate, shownBefore) {
  if (candidate.pixels === undefined || candidate.pixels === null) return true;
  const grow = shownBefore ? 1 + HOLD_SHARE : 1;
  const min = candidate.minPixels;
  const max = candidate.maxPixels;
  if (min !== undefined && min !== null && candidate.pixels < min / grow) return false;
  if (max !== undefined && max !== null && candidate.pixels > max * grow) return false;
  return true;
}

// Places `candidates` and returns the shown ones, [{candidate, left, top, width, height}], most important first.
// options: project, measure, viewport {width, height}, budget, padding, previous (a Set of the ids shown last).
export function placeLabels(candidates, options) {
  const budget = options.budget !== undefined ? options.budget : DEFAULT_BUDGET;
  const padding = options.padding !== undefined ? options.padding : PADDING_PX;
  const previous = options.previous || new Set();
  const width = options.viewport.width;
  const height = options.viewport.height;
  const ranked = [];
  candidates.forEach(function (candidate) {
    const held = previous.has(candidate.id);
    if (!inRange(candidate, held)) return;
    const spot = options.project(candidate.x, candidate.y, candidate.z);
    if (!spot || !(spot.x >= 0 && spot.x <= width && spot.y >= 0 && spot.y <= height)) return;
    ranked.push({ candidate: candidate, spot: spot, held: held });
  });
  ranked.sort(function (a, b) {
    if (b.candidate.priority !== a.candidate.priority) return b.candidate.priority - a.candidate.priority;
    if (a.held !== b.held) return a.held ? -1 : 1;
    const da = a.spot.depth === undefined ? 0 : a.spot.depth;
    const db = b.spot.depth === undefined ? 0 : b.spot.depth;
    if (da !== db) return da - db;
    return a.candidate.id < b.candidate.id ? -1 : a.candidate.id > b.candidate.id ? 1 : 0;
  });
  const grid = new Map();
  const placed = [];
  function cells(box) {
    const keys = [];
    const x0 = Math.floor(box.left / GRID_PX);
    const x1 = Math.floor((box.left + box.width) / GRID_PX);
    const y0 = Math.floor(box.top / GRID_PX);
    const y1 = Math.floor((box.top + box.height) / GRID_PX);
    for (let gx = x0; gx <= x1; gx++) for (let gy = y0; gy <= y1; gy++) keys.push(gx + "," + gy);
    return keys;
  }
  function collides(box, keys) {
    for (let k = 0; k < keys.length; k++) {
      const list = grid.get(keys[k]);
      if (!list) continue;
      for (let n = 0; n < list.length; n++) {
        const other = list[n];
        if (box.left < other.left + other.width + padding && other.left < box.left + box.width + padding
          && box.top < other.top + other.height + padding && other.top < box.top + box.height + padding) return true;
      }
    }
    return false;
  }
  for (let i = 0; i < ranked.length && placed.length < budget; i++) {
    const item = ranked[i];
    const w = options.measure(item.candidate.text) + 2 * GAP_PX;
    const box = { left: item.spot.x - w / 2, top: item.spot.y - LABEL_HEIGHT_PX - GAP_PX, width: w, height: LABEL_HEIGHT_PX };
    if (box.left < 0) box.left = 0;
    if (box.left + box.width > width) box.left = Math.max(0, width - box.width);
    if (box.top < 0) box.top = item.spot.y + GAP_PX;
    const keys = cells(box);
    if (collides(box, keys)) continue;
    keys.forEach(function (key) {
      if (!grid.has(key)) grid.set(key, []);
      grid.get(key).push(box);
    });
    placed.push({ candidate: item.candidate, left: box.left, top: box.top, width: box.width, height: box.height });
  }
  return placed;
}

// The labels drawn over a map. `container` is an element laid over the map (position absolute, no pointer
// events); `project`, `viewport` (a function giving {width, height}) and `reducedMotion` come from the map.
export function createLabelLayer(options) {
  const container = options.container;
  const measure = options.measure || createMeasure("12px sans-serif");
  const pool = [];
  const live = new Map(); // id -> {element, kind}
  let shownIds = new Set();
  let lastPlaced = [];

  function elementFor(id) {
    let entry = live.get(id);
    if (entry) return entry;
    const element = pool.pop() || document.createElement("span");
    element.className = "map-label";
    element.style.opacity = "0";
    container.appendChild(element);
    entry = { element: element, shown: false };
    live.set(id, entry);
    return entry;
  }

  // A label that is no longer placed fades out (at once under reduced motion), then its element goes back to the pool.
  function retire(id, entry) {
    function remove() {
      if (live.get(id) !== entry) return;
      live.delete(id);
      entry.element.remove();
      entry.element.textContent = "";
      if (pool.length < 200) pool.push(entry.element);
    }
    entry.shown = false;
    entry.element.style.opacity = "0";
    if (options.reducedMotion && options.reducedMotion()) {
      remove();
      return;
    }
    entry.dropping = setTimeout(function () { entry.dropping = null; remove(); }, FADE_MS);
  }

  // Places `candidates` and brings the DOM to match; returns the placed list.
  function update(candidates) {
    const view = options.viewport();
    const placed = placeLabels(candidates, {
      project: options.project, measure: measure, viewport: view, budget: options.budget, previous: shownIds,
    });
    const next = new Set();
    placed.forEach(function (item) {
      const id = item.candidate.id;
      next.add(id);
      const entry = elementFor(id);
      if (entry.dropping) {
        clearTimeout(entry.dropping);
        entry.dropping = null;
      }
      const element = entry.element;
      if (element.textContent !== item.candidate.text) element.textContent = item.candidate.text;
      element.dataset.kind = item.candidate.kind;
      element.style.transform = "translate(" + Math.round(item.left) + "px, " + Math.round(item.top) + "px)";
      if (!entry.shown) {
        entry.shown = true;
        if (options.reducedMotion && options.reducedMotion()) element.style.opacity = "1";
        else requestAnimationFrame(function () { element.style.opacity = "1"; });
      }
    });
    live.forEach(function (entry, id) {
      if (next.has(id) || entry.dropping) return;
      retire(id, entry);
    });
    shownIds = next;
    lastPlaced = placed;
    return placed;
  }

  // For the browser tests: what is placed now, [{id, kind, text}].
  function shown() {
    return lastPlaced.map(function (item) { return { id: item.candidate.id, kind: item.candidate.kind, text: item.candidate.text }; });
  }

  function clear() {
    update([]);
  }

  return { update: update, shown: shown, clear: clear };
}
