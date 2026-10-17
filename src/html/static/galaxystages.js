// html/static/galaxystages.js
//
// The Galaxy Map's stages (docs/design/galaxy-drilldown-navigation.md, as
// MAP.150 and MAP.151 changed it; docs/design/spherical-block-ladder.md for
// the shell grid), without any drawing: which blocks a stage holds and how
// many sectors each can hold, stage URLs, the breadcrumb, labels, and the
// camera path of a flight into a block. galaxymap3d.js draws and drives them.
// Like galaxyprisms.js it has no three.js import, so it runs under plain node
// for tests (src/tests/test_galaxystages.py).
//
// A stage is {win}: see "Stages" below.

const VERSION_QUERY = new URL(import.meta.url).search;
const { formatNumber } = await import(`./numberformat.js${VERSION_QUERY}`);

const {
  DEFAULT_WINDOW, DRILL_LEVELS, WINDOW_WIDTHS, bandRange, bandColatitudeBounds, boundsContain, cellAt, cellCenter, cellExists,
  cellPositionPc, drillBandCount, drillBlockBounds, drillBlockHolds, drillBlockTotal, drillChainOf, drillParent,
  drillWedgeCount, formatDrillKey, galaxyOutline, isBlock, outlineHigh, parseSectorDesignation, sameBlock,
  sectorAddressAt, sectorDesignation, sectorInOutline, wedgeCount, windowCells,
} = await import(`./galaxyprisms.js${VERSION_QUERY}`);

// The flight's zoom/pan trade-off (van Wijk and Nuij's rho), and its
// duration per unit of path length, clamped (section 5.3).
export const FLIGHT_RHO = 1.4;
export const FLIGHT_SPEED = 1.1;
export const FLIGHT_SECONDS = [0.5, 1.6];
// A view fits its blocks' extent with this much margin.
export const FIT_MARGIN = 1.1;

const TWO_PI = 2 * Math.PI;

// --- The galaxy's outline -----------------------------------------------------

// The outline the browser can work out for itself: for each sector shell,
// the highest |band| allowed (-1 for none), from the analytic density bound
// (galaxyprisms.galaxyOutline).
export { galaxyOutline };

// The block shells `outline` reaches at size m, and for each the rows that can hold a sector: a row is
// worth looking at when it meets the colatitudes the outline allows in any of its sector shells.
function blockRows(shell, m, outline) {
  let half = 0;
  for (let k = shell * m; k < shell * m + m; k++) {
    const high = outlineHigh(outline, k);
    if (high < 0) continue;
    if (k === 0) {
      half = Math.max(half, Math.PI);
      continue;
    }
    const north = bandColatitudeBounds(k, high)[0];
    half = Math.max(half, Math.PI / 2 - north);
  }
  if (half <= 0) return [];
  const rows = drillBandCount(m, shell);
  const low = Math.PI / 2 - half;
  const high = Math.PI / 2 + half;
  const out = [];
  for (let row = 0; row < rows; row++) {
    if ((row + 1) * Math.PI / rows > low && row * Math.PI / rows < high) out.push((rows - 1) / 2 - row);
  }
  return out;
}

// --- Stages -------------------------------------------------------------------
//
// A stage is {win}: `win` is null for the galaxy (every level-243 block of it,
// a click on one opens the window one size down) or a window {centre, g}: `g`
// cells across around `centre`, a cell of any level (a block of 243, 27 or 3
// shells deep, or a sector as {m: 1, shell, band, wedge}).
// There is nothing to pick but a cell: a click on one takes the window one
// level down, centred on the cell clicked (a sector is opened instead), and
// Up takes it one level up, centred on the cell holding the centre (the
// galaxy above level 243). A window's data is /galaxy/window
// (planetgen/galaxy/window.py); the galaxy's is /galaxy/stage.

const galaxyCache = new WeakMap();
const totalCache = new WeakMap();

export function galaxyStage() {
  return { win: null };
}

export { sameBlock };

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

// A block without the bounds a view adds.
export function stripBlock(block) {
  return { m: block.m, shell: block.shell, band: block.band, wedge: block.wedge };
}

// How many sectors `block` can hold under `outline`, or null without a shape. Worked out when asked (a big
// block counts every band of every shell it holds) and kept.
export function blockTotal(block, outline) {
  if (outline.shapeless) return null;
  let cache = totalCache.get(outline);
  if (!cache) {
    cache = new Map();
    totalCache.set(outline, cache);
  }
  const key = formatDrillKey(block);
  let total = cache.get(key);
  if (total === undefined) {
    total = drillBlockTotal(block, outline);
    cache.set(key, total);
  }
  return total;
}

// The galaxy's own blocks (the level-243 ones that can hold a sector), each
// with its `bounds`; the sectors one holds are `blockTotal`'s to count. Worked
// out once per outline.
export function galaxyBlocks(outline, edgePc) {
  let blocks = galaxyCache.get(outline);
  if (!blocks) {
    blocks = [];
    const top = DRILL_LEVELS[0];
    for (let shell = 0; shell <= Math.floor(outline.outerShell / top); shell++) {
      blockRows(shell, top, outline).forEach(function (band) {
        const wedges = drillWedgeCount(top, shell, band);
        for (let wedge = 0; wedge < wedges; wedge++) {
          const block = { m: top, shell: shell, band: band, wedge: wedge };
          if (!drillBlockHolds(block, outline)) continue;
          block.bounds = drillBlockBounds(block, edgePc);
          blocks.push(block);
        }
      });
    }
    galaxyCache.set(outline, blocks);
  }
  return blocks;
}

// The window's blocks with their bounds (those that can hold no sector are left out), as galaxyBlocks makes the
// galaxy's.
export function windowView(win, outline, edgePc) {
  const blocks = [];
  windowCells(win.centre, edgePc, win.g).forEach(function (block) {
    if (!drillBlockHolds(block, outline)) return;
    block.bounds = drillBlockBounds(block, edgePc);
    blocks.push(block);
  });
  return { win: win, blocks: blocks };
}

// The choices of a view: each block is one, a click on it (windowInto) opens
// the window one size down.
function cellOptions(blocks) {
  return blocks.map(function (block) {
    return { pick: { kind: "cell", block: stripBlock(block) }, blocks: [block] };
  });
}

// Everything about `stage`: {stage (as given), view: {win, blocks}, kind
// ("cell"), options (one per block), problem (when the window holds nothing
// of the galaxy)}.
export function resolveStage(stage, outline, edgePc) {
  if (!stage.win) {
    const blocks = galaxyBlocks(outline, edgePc);
    return { stage: stage, view: { win: null, blocks: blocks }, kind: "cell", options: cellOptions(blocks), problem: null };
  }
  const view = windowView(stage.win, outline, edgePc);
  if (!view.blocks.length) {
    return { stage: stage, problem: "Nothing of the galaxy lies around " + windowLabel(stage.win) + "." };
  }
  return { stage: stage, view: view, kind: "cell", options: cellOptions(view.blocks), problem: null };
}

// The cell one size down that a click on `block` opens a window round: the one under `at` (a point
// {x, y, z} in parsecs, where the click landed), or holding the block's middle without one. A block
// holds the cell it was clicked in unless a rounding put the point on the block's edge, when the
// block's middle serves. Null for a sector.
export function windowChildCentre(block, at, edgePc) {
  if (block.m === 1) return null;
  const childM = block.m / 3;
  const b = block.bounds || drillBlockBounds(block, edgePc);
  const middle = cellCenter(b);
  const inside = function (point) {
    const child = cellAt(childM, point[0], point[1], point[2], edgePc);
    // A sector joins the level-3 block holding its centre, so the sector at a point need not be in the block.
    if (childM === 1 || sameBlock(drillParent(child), block)) return child;
    return null;
  };
  const found = (at && inside([at.x, at.y, at.z])) || inside(middle);
  if (found) return found;
  return cellAt(childM, middle[0], middle[1], middle[2], edgePc);
}

// Where a click went through a block, so a window opens round the middle of what was pointed at and not
// round the face the ray first met (the top of a block that reaches high above the plane): the middle of the
// stretch of the ray, from `eye` through `at` (both {x, y, z}), that lies inside `bounds`. `at` itself when the ray
// misses the block (a cell picked by its column from above) or has no direction.
export function chordMiddle(bounds, at, eye) {
  if (!at || !eye) return at;
  const d = [at.x - eye.x, at.y - eye.y, at.z - eye.z];
  const norm = Math.hypot(d[0], d[1], d[2]);
  if (!(norm > 0)) return at;
  const centre = cellCenter(bounds);
  const size = Math.max(bounds.r1 - bounds.r0, bounds.r1 * (bounds.c1 - bounds.c0),
    bounds.r1 * Math.sin((bounds.c0 + bounds.c1) / 2) * (bounds.t1 - bounds.t0));
  const reach = 2 * size + Math.hypot(at.x - centre[0], at.y - centre[1], at.z - centre[2]);
  const steps = 512;
  let first = null;
  let last = null;
  for (let k = -steps; k <= steps; k++) {
    const t = (reach * k) / steps / norm;
    const x = at.x + d[0] * t;
    const y = at.y + d[1] * t;
    const z = at.z + d[2] * t;
    if (boundsContain(bounds, x, y, z)) {
      if (!first) first = { x: x, y: y, z: z };
      last = { x: x, y: y, z: z };
    }
  }
  if (!first) return at;
  return { x: (first.x + last.x) / 2, y: (first.y + last.y) / 2, z: (first.z + last.z) / 2 };
}

export function windowStage(centre, g) {
  return { win: { centre: stripBlock(centre), g: g || DEFAULT_WINDOW } };
}

// The window stage a click on `block` (at the point `at`, see windowChildCentre) opens, or null on a
// sector.
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

// Whether sector (shell, band, wedge) can hold anything under `outline`.
function sectorInGalaxy(shell, band, wedge, outline) {
  return sectorInOutline(shell, band, wedge, outline);
}

// The stage that shows a sector among its neighbours (MAP.26): the window of
// sectors centred on it. Null when the sector is outside the galaxy.
export function sectorStage(shell, band, wedge, outline) {
  if (!sectorInGalaxy(shell, band, wedge, outline)) return null;
  return windowStage({ m: 1, shell: shell, band: band, wedge: wedge }, DEFAULT_WINDOW);
}

// The smallest stage showing every one of `sectors` ({shell, band, wedge}),
// for the NAV course overlay (section 9.4): the window of the finest level
// whose DEFAULT_WINDOW cells around the course's middle hold the cell of
// every stop. An empty list, or a course too wide for a window of 243-shell
// blocks, gives the galaxy.
export function courseStage(sectors, outline, edgePc) {
  const galaxy = galaxyStage();
  if (!sectors || !sectors.length) return galaxy;
  const stops = sectors.filter(function (s) { return sectorInGalaxy(s.shell, s.band, s.wedge, outline); });
  if (stops.length !== sectors.length) return galaxy;
  const middle = [0, 0, 0];
  stops.forEach(function (s) {
    const p = cellPositionPc({ m: 1, shell: s.shell, band: s.band, wedge: s.wedge }, edgePc);
    for (let k = 0; k < 3; k++) middle[k] += p[k] / stops.length;
  });
  const at = sectorAddressAt(middle[0], middle[1], middle[2], edgePc);
  const around = drillChainOf(at.shell, at.band, at.wedge).reverse();
  const chains = stops.map(function (s) { return drillChainOf(s.shell, s.band, s.wedge).reverse(); });
  for (let level = 0; level < around.length; level++) {
    const wanted = new Set(windowCells(around[level], edgePc, DEFAULT_WINDOW).map(formatDrillKey));
    if (chains.every(function (chain) { return wanted.has(formatDrillKey(chain[level])); })) {
      return windowStage(around[level], DEFAULT_WINDOW);
    }
  }
  return galaxy;
}

// "m.shell.band.wedge" of any level, a sector (m = 1) included, or null.
export function parseCellKey(key) {
  const parts = String(key == null ? "" : key).split(".");
  if (parts.length !== 4 || !parts.every(function (p) { return /^-?\d+$/.test(p); })) return null;
  const v = parts.map(Number);
  const cell = { m: v[0], shell: v[1], band: v[2], wedge: v[3] };
  return isBlock(cell) ? cell : null;
}

export function windowLabel(win) {
  return (win.centre.m === 1 ? "Sector " : "Block ") + blockLabel(win.centre) + " (" + win.centre.m + ")";
}

// The 12 edges of a cell given as bounds as line segments, a flat [x, y, z, x, y, z, ...] list of pairs: the
// four arcs of longitude (at each radius and colatitude), the four of colatitude (at each radius and
// longitude) and the four radial lines. Edges that collapse onto the axis or the centre are left out.
export function cellEdgePoints(b) {
  const points = [];
  const at = function (r, c, t) {
    const plane = r * Math.sin(c);
    return [plane * Math.cos(t), plane * Math.sin(t), r * Math.cos(c)];
  };
  const arc = function (r, c0, c1, t0, t1) {
    const turn = Math.max(Math.abs(t1 - t0) * Math.sin((c0 + c1) / 2), Math.abs(c1 - c0));
    const steps = Math.max(1, Math.ceil(turn / (Math.PI / 90)));
    for (let k = 0; k < steps; k++) {
      const f0 = k / steps;
      const f1 = (k + 1) / steps;
      Array.prototype.push.apply(points, at(r, c0 + (c1 - c0) * f0, t0 + (t1 - t0) * f0));
      Array.prototype.push.apply(points, at(r, c0 + (c1 - c0) * f1, t0 + (t1 - t0) * f1));
    }
  };
  const fullTurn = b.t1 - b.t0 >= TWO_PI - 1e-9;
  [b.r0, b.r1].forEach(function (r) {
    if (r <= 0) return;
    [b.c0, b.c1].forEach(function (c) {
      if (Math.sin(c) > 1e-12) arc(r, c, c, b.t0, b.t1);
    });
    if (!fullTurn) {
      [b.t0, b.t1].forEach(function (t) { arc(r, b.c0, b.c1, t, t); });
    }
  });
  if (!fullTurn) {
    [b.c0, b.c1].forEach(function (c) {
      [b.t0, b.t1].forEach(function (t) {
        if (Math.sin(c) <= 1e-12) return;
        Array.prototype.push.apply(points, at(b.r0, c, t));
        Array.prototype.push.apply(points, at(b.r1, c, t));
      });
    });
  }
  return points;
}

// Whether the vertical line through (x, y) meets a cell given as bounds: whether the whole galaxy, seen from
// above, shows the cell where a ray comes down at that point. A cell's points lie at a distance from the axis
// of r sin(colatitude), which runs from the inner radius at its less steep colatitude to the outer radius at the
// steepest (the plane's, if it holds that).
export function columnContains(b, x, y) {
  const reach = Math.hypot(x, y);
  const s0 = Math.sin(b.c0);
  const s1 = Math.sin(b.c1);
  const low = b.r0 * Math.min(s0, s1);
  const high = b.r1 * (b.c0 <= Math.PI / 2 && b.c1 >= Math.PI / 2 ? 1 : Math.max(s0, s1));
  if (reach < low || reach > high) return false;
  if (b.t1 - b.t0 >= TWO_PI - 1e-9 || reach === 0) return true;
  const theta = Math.atan2(y, x);
  const t = theta < 0 ? theta + TWO_PI : theta + 0;
  return t >= b.t0 && t < b.t1;
}

// Points round a set of blocks' bounds for fitting a view to them: each block's corners and, along its outer
// and inner arcs, a point at most 22.5 degrees apart. A flat [x, y, z, ...] list.
export function boundsPoints(blocks) {
  const points = [];
  blocks.forEach(function (block) {
    const b = block.bounds;
    const lon = Math.max(1, Math.ceil((b.t1 - b.t0) / (Math.PI / 8)));
    const lat = Math.max(1, Math.ceil((b.c1 - b.c0) / (Math.PI / 8)));
    [b.r0, b.r1].forEach(function (r) {
      if (r <= 0) {
        points.push(0, 0, 0);
        return;
      }
      for (let j = 0; j <= lat; j++) {
        const c = b.c0 + ((b.c1 - b.c0) * j) / lat;
        for (let i = 0; i <= lon; i++) {
          const t = b.t0 + ((b.t1 - b.t0) * i) / lon;
          const plane = r * Math.sin(c);
          points.push(plane * Math.cos(t), plane * Math.sin(t), r * Math.cos(c));
        }
      }
    });
  });
  return points;
}

// The extent of `blocks`: its center (the middle of the bounding box of their bounds' points) and the radius
// that reaches every point from it.
export function extentOf(blocks) {
  const points = boundsPoints(blocks);
  const low = [Infinity, Infinity, Infinity];
  const high = [-Infinity, -Infinity, -Infinity];
  for (let k = 0; k + 2 < points.length; k += 3) {
    for (let a = 0; a < 3; a++) {
      low[a] = Math.min(low[a], points[k + a]);
      high[a] = Math.max(high[a], points[k + a]);
    }
  }
  const center = [(low[0] + high[0]) / 2, (low[1] + high[1]) / 2, (low[2] + high[2]) / 2];
  let radius = 0;
  for (let k = 0; k + 2 < points.length; k += 3) {
    radius = Math.max(radius, Math.hypot(points[k] - center[0], points[k + 1] - center[1], points[k + 2] - center[2]));
  }
  return { center: center, radius: radius };
}

// --- Designations -------------------------------------------------------------

export { sectorDesignation, parseSectorDesignation };

// --- The address bar (section 9.3) -------------------------------------------

// What the address bar's text asks for:
// - {sector: {shell, band, wedge}} for a designation ("1AA5FDEC0C"),
//   "312/-3/1042" or "shell 312 band -3 wedge 1042";
// - {sector, point: [x, y, z]} for "x, y, z" in pc (with or without
//   "pc"), the sector holding that point;
// - {name} for anything else (looked up by the server);
// - {problem} for an address that can't be a sector (a wedge past its
//   band's end, say), or empty text.
// Designations are 7 to 10 hex digits (the band and wedge fields alone take 7), so a short
// name made of hex letters ("Bead") stays a name.
export function parseAddress(text, edgePc) {
  const raw = String(text == null ? "" : text).trim();
  if (!raw) return { problem: "Type a designation, shell/band/wedge, x, y, z in pc, or a name." };
  const number = "(-?\\d+)";
  let match = raw.match(new RegExp("^" + number + "\\s*/\\s*" + number + "\\s*/\\s*" + number + "$"))
    || raw.match(new RegExp("^shell\\s+" + number + "[\\s,]+band\\s+" + number + "[\\s,]+wedge\\s+" + number + "$", "i"));
  if (match) {
    const sector = { shell: Number(match[1]), band: Number(match[2]), wedge: Number(match[3]) };
    if (!cellExists(sector.shell, sector.band, sector.wedge)) return { problem: addressProblem(sector) };
    return { sector: sector };
  }
  const real = "(-?\\d+(?:\\.\\d+)?)";
  match = raw.match(new RegExp("^\\(?\\s*" + real + "\\s*[,\\s]\\s*" + real + "\\s*[,\\s]\\s*" + real + "\\s*\\)?\\s*(?:pc)?$", "i"));
  if (match) {
    const point = [Number(match[1]), Number(match[2]), Number(match[3])];
    return { sector: sectorAddressAt(point[0], point[1], point[2], edgePc), point: point };
  }
  if (/^[0-9A-Fa-f]{7,14}$/.test(raw)) {
    const sector = parseSectorDesignation(raw);
    if (!sector) return { problem: "There is no sector " + raw.toUpperCase() + "." };
    return { sector: sector };
  }
  return { name: raw };
}

// Why shell/band/wedge is not a cell of the grid.
function addressProblem(sector) {
  if (!(sector.shell >= 0) || sector.shell > 83037) return "There is no shell " + sector.shell + ".";
  const range = bandRange(sector.shell);
  if (sector.band < range[0] || sector.band > range[1]) {
    return "Shell " + sector.shell + " has bands " + range[0] + " to " + range[1] + ", not " + sector.band + ".";
  }
  return "Band " + sector.band + " of shell " + sector.shell + " has no wedge " + sector.wedge + ".";
}

// --- URLs ---------------------------------------------------------------------

// A stage URL's query (section 8.1): "" for the galaxy, else "?win=" with the
// centre's "m.shell.band.wedge" and "&g=" the cells across (when not the default).
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
// ({shell, band, wedge}) is set by ?sector=<designation>, whose stage the
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
    // The window the visitor was in when the sector was picked, if the link keeps it (`?win=...&sector=...`).
    const asked = params.has("win") ? parseCellKey(params.get("win")) : null;
    const g = params.has("g") ? Number(params.get("g")) : DEFAULT_WINDOW;
    const stage = asked && WINDOW_WIDTHS.indexOf(g) >= 0 ? windowStage(asked, g) : galaxy;
    return { stage: stage, sector: sector, open: params.get("open") === "1", system: system, object: object, problem: null };
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

// A label's whole number, never scientific (shell, band and wedge numbers
// are addresses, not amounts).
export function formatInt(n) {
  return Math.round(n).toLocaleString("en-US");
}

// A count of sectors: scientific past 4 digits (UX.20).
export function formatCount(n) {
  return formatNumber(Math.round(n));
}

// "7·-2·14": a block's or a sector's shell, band and wedge.
export function blockLabel(block) {
  return [formatInt(block.shell), block.band, formatInt(block.wedge)].join("·");
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

// MAP.227: the stage the camera is in. A stage is a window round one cell; the camera's middle may have moved to
// another cell of the same size, and the trail of containers is that cell's, worked out from the position (the
// cells are closed form) and not from the stage the view was opened at. The stage itself when the middle is still
// in its cell, or the cell it moved to lies outside the galaxy.
export function stageAtPoint(stage, point, outline, edgePc) {
  if (!stage.win) return stage;
  const centre = stage.win.centre;
  const cell = cellAt(centre.m, point[0], point[1], point[2], edgePc);
  if (sameBlock(cell, centre)) return stage;
  const next = windowStage(cell, stage.win.g);
  return resolveStage(next, outline, edgePc).problem ? stage : next;
}

// The sector holding `point` as {shell, band, wedge}, or null when the galaxy has none there (outside its outline).
export function sectorAtPoint(point, outline, edgePc) {
  const at = sectorAddressAt(point[0], point[1], point[2], edgePc);
  return sectorInOutline(at.shell, at.band, at.wedge, outline) ? at : null;
}

// Degrees counterclockwise from +X, as the longitude labels read.
export function bearingRange(bounds) {
  const deg = function (rad) { return (((rad * 180) / Math.PI) % 360 + 360) % 360; };
  return [deg(bounds.t0), deg(bounds.t1) || 360];
}

// Degrees north of the galactic plane (negative south) at a cell's northern and southern edges.
export function latitudeRange(bounds) {
  const deg = function (colat) { return 90 - (colat * 180) / Math.PI; };
  return [deg(bounds.c0), deg(bounds.c1)];
}

// --- The flight ---------------------------------------------------------------

// Van Wijk and Nuij's smooth zoom-and-pan from view (c0, w0) to (c1, w1):
// c a 2D center, w the visible width. Returns {S, at(s) -> {center, w}}
// over path length s in [0, S] (section 5.3). A flight in space runs it along the
// straight line between the two centers, as c0 = [0, 0] and c1 = [distance, 0].
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

export { drillBlockBounds, drillParent, drillWedgeCount, formatDrillKey, drillBlockHolds, drillBlockTotal };
