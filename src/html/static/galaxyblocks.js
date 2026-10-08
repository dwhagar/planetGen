// html/static/galaxyblocks.js
//
// The Galaxy Map's block scene (galaxymap3d.js): for one view, which
// blocks to draw (galaxyprisms.blocksForView, plus the interior blocks
// holding generated sectors), their colors and opacities, and their
// geometry, packed as typed arrays ready for the GPU.
//
// This is the slow part of a zoom step, so the page runs it in a Web
// Worker: loaded as a worker, this module answers "init", "filled" and
// "build" messages (see the bottom of the file) and hands the arrays back
// as transferables. The page also imports it directly, to draw the first
// frame without waiting on the worker and as the fallback where workers
// can't start. Like galaxyprisms.js it has no three.js import, so it runs
// under plain node for tests (src/tests/test_galaxyprisms.py).
//
// Colors come in as linear [r, g, b] (the page reads them off
// THREE.Color, which keeps linear components), and every blend here is
// the same component-wise lerp THREE.Color.lerp does. The fade toward the
// view's rim and the near cut are left to the shader, so a built scene
// stays right while the camera zooms and turns within it.

const VERSION_QUERY = new URL(import.meta.url).search;

// In a worker, messages can arrive while the imports below are still
// loading; they wait here until the module is ready.
const IN_WORKER = typeof WorkerGlobalScope !== "undefined" && self instanceof WorkerGlobalScope;
const waiting = [];
if (IN_WORKER) {
  self.onmessage = function (event) {
    waiting.push(event.data);
  };
}

const {
  blockAddressAt, blockAt, blockSectorCount, blockSizeForScale, blocksForView, boundsDensity, buildPrismGeometry,
  cellCoordinates,
} = await import(`./galaxyprisms.js${VERSION_QUERY}`);

// Unfilled space: this opaque at the sparsest drawn density, rising to
// BLOCK_OPACITY_DENSE at the densest (90% to 70% see-through). Boss
// (MAP.37): "everything not yet filled should be more transparent by a
// lot", so generated systems stand out; the density shape still reads,
// only faintly.
export const BLOCK_OPACITY_SPARSE = 0.1;
export const BLOCK_OPACITY_DENSE = 0.3;
// Any filled sector lifts a block at least this far from its unfilled
// opacity toward solid, however tiny its share: at large m one filled
// sector among 531,441 must still show, and well clear of the unfilled
// blocks around it (MAP.37). The rest of the way follows the
// log of the filled count against the log of the block's total, so a
// block turns solid only when every sector in it is filled.
export const FILLED_MIN_STEP = 0.6;
// MAP.128/129: a drill-down block colored by what its generated sectors
// hold (each cell's `stats`, from the stage API's sector_stats sums).
// Boss (2026-10-08 01:59Z): density sets the transparency (denser is more
// solid, never fully opaque), the stars' mean age sets the hue, and the
// summed luminosity sets the brightness. The color is the color of the
// translucent fill.
// A generated sector with stars is this opaque at the sparsest and
// FILLED_OPACITY_DENSE at the densest; one with no stars lifts the
// unfilled opacity by FILLED_EMPTY_STEP and is drawn in the unfilled
// color a shade more saturated.
export const FILLED_OPACITY_SPARSE = 0.35;
export const FILLED_OPACITY_DENSE = 0.85;
export const FILLED_EMPTY_STEP = 0.05;
const FILLED_EMPTY_SATURATION = 0.15;
// Stellar age (Gy) to hue, sRGB: blue young, slate at the disk's 4.5 Gy
// average, amber old (OLD_AGE_GY and older).
export const DISK_AGE_GY = 4.5;
export const OLD_AGE_GY = 10;
const AGE_COLORS = { young: [0.3, 0.52, 1.0], disk: [0.56, 0.6, 0.68], old: [1.0, 0.68, 0.22] };
// Brightness runs from this share of the hue at the dimmest sector in
// view up to all of it at the brightest, by the log of the luminosity.
export const BRIGHTNESS_FLOOR = 0.3;
// Block shading runs over a wide density range: the drawing floor
// (galaxyprisms.js's PRISM_MIN_DENSITY) up to the core, on a log scale,
// through a dim-to-accent-to-white ramp.
const PRISM_DENSITY_LOW = 0.02;
const PRISM_DENSITY_HIGH = 100;
// Density alone puts the arms (1 +/- arm_amplitude around the ring's
// mean) on a sliver of that 0.02-100 ramp, so the spiral barely shows.
// Each block's shade mixes its density with its arm factor (density /
// its azimuthal mean), stretched over the arm model's own range:
// inter-arm troughs dim, arm crests bright.
const PRISM_DENSITY_SHARE = 0.45;
// The light the blocks' faces are shaded by: from galactic north, a
// little off to one side, so tops, walls and sides all read apart.
const PRISM_LIGHT = normalize([0.35, -0.3, 0.9]);
const PRISM_AMBIENT = 0.35;
// Densities already worked out, kept for this many block sizes.
const DENSITY_CACHE_SIZES = 4;
// blockSectorCount results kept between builds (at large m it walks every
// member ring and layer).
const BLOCK_TOTALS_MAX = 20000;

// Per block in a built part's `cells`: ring, seg, slab, density,
// meanDensity, filled, total, centroid x, y, z, and the generated sector's
// own point (its index among the filled points, at one sector per block)
// -- NaN where there is none.
export const CELL_STRIDE = 11;
// Per filled point (setFilled): x, y, z, count, 1 when the point is one
// generated sector (else 0), and that sector's own relative density
// (placedRelativeDensity; NaN when unknown).
export const POINT_STRIDE = 6;

function normalize(v) {
  const length = Math.hypot(v[0], v[1], v[2]);
  return [v[0] / length, v[1] / length, v[2] / length];
}

function lerp3(a, b, t) {
  return [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
}

// 0 for a block with nothing filled, FILLED_MIN_STEP for one filled
// sector among many, 1 once every sector in it is filled.
export function blockFillStep(cell) {
  if (!(cell.filled > 0)) {
    return 0;
  }
  if (cell.filled >= cell.total) {
    return 1;
  }
  const share = Math.log(1 + cell.filled) / Math.log(1 + cell.total);
  return FILLED_MIN_STEP + (1 - FILLED_MIN_STEP) * share;
}

export function prismIntensity(relativeDensity) {
  const t = Math.log(Math.max(relativeDensity, 1e-9) / PRISM_DENSITY_LOW) / Math.log(PRISM_DENSITY_HIGH / PRISM_DENSITY_LOW);
  return Math.max(0, Math.min(1, t));
}

export function blockOpacity(cell) {
  const base = BLOCK_OPACITY_SPARSE + (BLOCK_OPACITY_DENSE - BLOCK_OPACITY_SPARSE) * prismIntensity(cell.density || 0);
  return base + (1 - base) * blockFillStep(cell);
}

// MAP.128: how the cells in view compare: the log of the systems per
// generated sector (the density) and of the summed luminosity, each as a
// {min, max} over the cells that hold stars. `statsOpacity` and
// `statsColor` rank a cell against these.
export function statsRanges(cells) {
  const ranges = { density: { min: Infinity, max: -Infinity }, luminosity: { min: Infinity, max: -Infinity } };
  cells.forEach(function (cell) {
    const stats = cell.stats;
    if (!stats || !(stats.stars > 0) || !(cell.filled > 0)) {
      return;
    }
    const density = Math.log1p(stats.systems / cell.filled);
    const luminosity = Math.log10(Math.max(stats.luminosity_sol, 1e-9));
    ranges.density.min = Math.min(ranges.density.min, density);
    ranges.density.max = Math.max(ranges.density.max, density);
    ranges.luminosity.min = Math.min(ranges.luminosity.min, luminosity);
    ranges.luminosity.max = Math.max(ranges.luminosity.max, luminosity);
  });
  return ranges;
}

function rangeShare(value, range) {
  if (!(range.max > range.min)) {
    return 1;
  }
  return Math.max(0, Math.min(1, (value - range.min) / (range.max - range.min)));
}

// The density share, 0..1, of a cell holding stars.
function densityShare(cell, ranges) {
  return rangeShare(Math.log1p(cell.stats.systems / cell.filled), ranges.density);
}

// MAP.128/129: a block's opacity from its sectors. Unfilled space keeps
// its own (blockOpacity); a filled block moves toward the opacity its
// stars' density gives it (FILLED_OPACITY_SPARSE up to
// FILLED_OPACITY_DENSE, never solid) by blockFillStep, so one filled
// sector among many still shows. A filled block with no stars only
// lifts the unfilled opacity by FILLED_EMPTY_STEP.
export function statsOpacity(cell, ranges) {
  const base = BLOCK_OPACITY_SPARSE + (BLOCK_OPACITY_DENSE - BLOCK_OPACITY_SPARSE) * prismIntensity(cell.density || 0);
  if (!(cell.filled > 0)) {
    return base;
  }
  const stats = cell.stats;
  if (!stats || !(stats.stars > 0)) {
    return Math.min(FILLED_OPACITY_DENSE, base + FILLED_EMPTY_STEP);
  }
  const own = FILLED_OPACITY_SPARSE + (FILLED_OPACITY_DENSE - FILLED_OPACITY_SPARSE) * densityShare(cell, ranges);
  return base + (Math.max(own, base) - base) * blockFillStep(cell);
}

function srgbToLinear(v) {
  return v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4);
}

// `color` (linear) pushed `amount` of the way from its own grey toward
// full saturation.
function saturate(color, amount) {
  const grey = (color[0] + color[1] + color[2]) / 3;
  return color.map(function (v) { return Math.max(0, grey + (v - grey) * (1 + amount)); });
}

// The sRGB hue for a mean stellar age in Gy: blue at 0, slate at
// DISK_AGE_GY, amber from OLD_AGE_GY on.
export function ageColor(ageGy) {
  const age = Math.max(0, ageGy || 0);
  if (age <= DISK_AGE_GY) {
    return lerp3(AGE_COLORS.young, AGE_COLORS.disk, age / DISK_AGE_GY);
  }
  return lerp3(AGE_COLORS.disk, AGE_COLORS.old, Math.min(1, (age - DISK_AGE_GY) / (OLD_AGE_GY - DISK_AGE_GY)));
}

// MAP.128/129: a cell's color (linear): an unfilled block in `unfilled`
// (the density ramp's color); a filled one with stars in the hue of its
// stars' mean age, at the brightness of their summed luminosity against
// the other cells in view; a filled one without stars in `unfilled` a
// shade more saturated.
export function statsColor(cell, unfilled, ranges) {
  if (!(cell.filled > 0)) {
    return unfilled;
  }
  const stats = cell.stats;
  if (!stats || !(stats.stars > 0)) {
    return saturate(unfilled, FILLED_EMPTY_SATURATION);
  }
  const luminosity = rangeShare(Math.log10(Math.max(stats.luminosity_sol, 1e-9)), ranges.luminosity);
  const brightness = BRIGHTNESS_FLOOR + (1 - BRIGHTNESS_FLOOR) * luminosity;
  return ageColor(stats.mean_age_gy).map(function (v) { return srgbToLinear(v) * brightness; });
}

// The arm model's own amplitude, clamped to 0..1; 0 without arms.
export function galaxyArmAmplitude(shape) {
  if (!shape || !(shape.arm_count > 0)) {
    return 0;
  }
  return Math.min(1, Math.abs(shape.arm_amplitude || 0));
}

// 0..1 shade for a block: density alone when the galaxy has no arms.
export function prismShade(cell, armAmplitude) {
  const t = prismIntensity(cell.density);
  if (!(armAmplitude > 0) || !(cell.meanDensity > 0)) {
    return t;
  }
  const arm = (cell.density / cell.meanDensity - (1 - armAmplitude)) / (2 * armAmplitude);
  return PRISM_DENSITY_SHARE * t + (1 - PRISM_DENSITY_SHARE) * Math.max(0, Math.min(1, arm));
}

// Builds block scenes for one galaxy. config:
// - edgePc, galaxyRadius, shape (the density shape, or null);
// - minPx, budget (galaxyprisms.blocksForView's);
// - palette: linear [r, g, b] for dim, accent and hot (the density ramp)
//   and placedLow, placedHigh (a generated sector's own color).
// Returns {setFilled(points), build(view), buildCells(cells, eye, dim)};
// see each below.
export function createBlockScene(config) {
  const edgePc = config.edgePc || 1;
  const galaxyRadius = config.galaxyRadius;
  const shape = config.shape || null;
  const palette = config.palette;
  const armAmplitude = galaxyArmAmplitude(shape);
  const densityCaches = new Map();
  const blockTotals = new Map();
  let points = new Float64Array(0);

  function densityCacheFor(m) {
    let cache = densityCaches.get(m);
    if (cache) {
      densityCaches.delete(m);
    } else {
      cache = new Map();
    }
    densityCaches.set(m, cache);
    while (densityCaches.size > DENSITY_CACHE_SIZES) {
      densityCaches.delete(densityCaches.keys().next().value);
    }
    return cache;
  }

  function prismColor(t) {
    return t < 0.6 ? lerp3(palette.dim, palette.accent, t / 0.6) : lerp3(palette.accent, palette.hot, (t - 0.6) / 0.4);
  }

  // A generated sector's own block: its real density relative to the
  // local neighborhood's on a "log2(x+1)/3" dim-to-bright curve, through
  // a warm bronze-to-gold range.
  function placedColor(relative) {
    const t = isNaN(relative) ? 0.5 : Math.max(0, Math.min(1, Math.log2(relative + 1) / 3));
    return lerp3(palette.placedLow, palette.placedHigh, t);
  }

  // The generated sectors to count into blocks, POINT_STRIDE numbers each
  // (a sector's center, or a coarser cell's middle with its count).
  function setFilled(packed) {
    points = packed || new Float64Array(0);
  }

  // The blocks for a view (the solid's surface, plus every block in view
  // holding filled sectors, interior ones included). Each gets `filled`,
  // and filled blocks also `total`, `centroid` (the mean position of
  // their filled sectors) and, at one sector per block, `sectorPoint` and
  // `sectorRelative` (the sector's point index and relative density).
  function viewBlocks(view) {
    const options = {
      pcPerPixel: view.pcPerPixel, sliceZ: view.sliceZ, viewerZ: view.viewerZ,
      minPx: config.minPx, budget: view.budget || config.budget,
    };
    const listed = shape
      ? blocksForView(view.center, view.viewRadius, edgePc, galaxyRadius, shape, densityCacheFor, options)
      : { m: blockSizeForScale(view.pcPerPixel, edgePc, config.minPx, galaxyRadius), slice: null, blocks: [] };
    const m = listed.m;
    const byKey = new Map();
    listed.blocks.forEach(function (block) {
      byKey.set(block.ring + "/" + block.seg + "/" + block.slab, block);
    });
    const sums = new Map();
    for (let p = 0; p < points.length; p += POINT_STRIDE) {
      const x = points[p];
      const y = points[p + 1];
      const z = points[p + 2];
      const count = points[p + 3];
      const a = blockAddressAt(x, y, z, m, edgePc);
      if (listed.slice !== null && a.slab > listed.slice) {
        continue;
      }
      const key = a.ring + "/" + a.seg + "/" + a.slab;
      let sum = sums.get(key);
      if (!sum) {
        sum = { address: a, count: 0, x: 0, y: 0, z: 0, sector: -1 };
        sums.set(key, sum);
      }
      sum.count += count;
      sum.x += x * count;
      sum.y += y * count;
      sum.z += z * count;
      if (points[p + 4]) {
        sum.sector = p / POINT_STRIDE;
      }
    }
    const center = view.center;
    const reach = view.viewRadius + m * edgePc;
    sums.forEach(function (sum, key) {
      let block = byKey.get(key);
      if (!block) {
        const a = sum.address;
        block = blockAt(a.ring, a.seg, a.slab, m, edgePc, shape, shape ? densityCacheFor(m) : null);
        const mid = cellCoordinates(block).cartesian;
        if (Math.hypot(mid[0] - center[0], mid[1] - center[1], mid[2] - center[2]) > reach) {
          return;
        }
        listed.blocks.push(block);
      }
      block.filled = sum.count;
      const totalKey = m + "/" + key;
      let total = blockTotals.get(totalKey);
      if (total === undefined) {
        if (blockTotals.size > BLOCK_TOTALS_MAX) {
          blockTotals.clear();
        }
        total = blockSectorCount(block.ring, block.seg, block.slab, m, edgePc, shape);
        blockTotals.set(totalKey, total);
      }
      block.total = Math.max(sum.count, total);
      block.centroid = [sum.x / sum.count, sum.y / sum.count, sum.z / sum.count];
      if (m === 1 && sum.sector >= 0) {
        block.sectorPoint = sum.sector;
        block.sectorRelative = points[sum.sector * POINT_STRIDE + 5];
      }
    });
    return listed;
  }

  // One mesh's worth of blocks, packed: positions, centers (Float32),
  // colors (Uint16, normalized: dark linear colors need more than 8 bits),
  // uvs, alphas and fills (Uint8, normalized), indices, each vertex's
  // block (owners, for picking), and per-block records (CELL_STRIDE).
  // Normals only light the colors here, so they aren't kept. Opaque blocks
  // leave out the faces they share (galaxyprisms.buildPrismGeometry).
  let ranges = statsRanges([]);

  function pack(cells, opaque) {
    const built = buildPrismGeometry(cells, { skipShared: opaque });
    const count = built.owners.length;
    const colors = new Uint16Array(count * 3);
    const centers = new Float32Array(count * 3);
    const alphas = new Uint8Array(count);
    const fills = new Uint8Array(count);
    const uvs = new Uint8Array(count * 2);
    const normals = built.normals;
    const middles = [];
    const colorOf = [];
    const alphaOf = [];
    const fillOf = [];
    const records = new Float64Array(cells.length * CELL_STRIDE);
    cells.forEach(function (cell, n) {
      middles.push(cellCoordinates(cell).cartesian);
      const placed = cell.sectorPoint != null;
      const ramp = function () { return prismColor(shape ? prismShade(cell, armAmplitude) : 0.5); };
      colorOf.push("stats" in cell ? statsColor(cell, ramp(), ranges) : placed ? placedColor(cell.sectorRelative) : ramp());
      alphaOf.push(Math.round(cell.opacity * 255));
      fillOf.push(Math.round(blockFillStep(cell) * 255));
      const o = n * CELL_STRIDE;
      records[o] = cell.ring;
      records[o + 1] = cell.seg;
      records[o + 2] = cell.slab;
      records[o + 3] = cell.density == null ? NaN : cell.density;
      records[o + 4] = cell.meanDensity == null ? NaN : cell.meanDensity;
      records[o + 5] = cell.filled || 0;
      records[o + 6] = cell.total || 0;
      records[o + 7] = cell.centroid ? cell.centroid[0] : NaN;
      records[o + 8] = cell.centroid ? cell.centroid[1] : NaN;
      records[o + 9] = cell.centroid ? cell.centroid[2] : NaN;
      records[o + 10] = placed ? cell.sectorPoint : NaN;
    });
    for (let v = 0; v < count; v++) {
      const owner = built.owners[v];
      const c = colorOf[owner];
      const lit = normals[3 * v] * PRISM_LIGHT[0] + normals[3 * v + 1] * PRISM_LIGHT[1] + normals[3 * v + 2] * PRISM_LIGHT[2];
      const shade = (PRISM_AMBIENT + (1 - PRISM_AMBIENT) * Math.max(0, lit)) * 65535;
      colors[3 * v] = Math.min(65535, Math.round(c[0] * shade));
      colors[3 * v + 1] = Math.min(65535, Math.round(c[1] * shade));
      colors[3 * v + 2] = Math.min(65535, Math.round(c[2] * shade));
      const mid = middles[owner];
      centers[3 * v] = mid[0];
      centers[3 * v + 1] = mid[1];
      centers[3 * v + 2] = mid[2];
      alphas[v] = alphaOf[owner];
      fills[v] = fillOf[owner];
      uvs[2 * v] = Math.round(built.uvs[2 * v] * 255);
      uvs[2 * v + 1] = Math.round(built.uvs[2 * v + 1] * 255);
    }
    return {
      positions: built.positions, centers: centers, colors: colors, uvs: uvs, alphas: alphas, fills: fills,
      indices: built.indices, owners: built.owners, cells: records, vertexCount: count,
    };
  }

  // Builds one view. view:
  // - center ([x, y, z]), viewRadius: the ball to list;
  // - pcPerPixel, sliceZ (null for no slice), viewerZ, budget (optional,
  //   overriding config.budget): as galaxyprisms.blocksForView takes them;
  // - eye ([x, y, z]): the camera, to sort translucent blocks back to
  //   front from.
  // Returns {m, solid, glass}: blocks whose every sector is generated
  // (opaque), and the rest (translucent, sorted far to near), each packed.
  function build(view) {
    const listed = viewBlocks(view);
    const solid = [];
    const glass = [];
    listed.blocks.forEach(function (cell) {
      cell.opacity = blockOpacity(cell);
      (cell.opacity >= 1 ? solid : glass).push(cell);
    });
    const eye = view.eye || view.center;
    glass.forEach(function (cell) {
      const mid = cellCoordinates(cell).cartesian;
      cell.eyeDistance = Math.hypot(mid[0] - eye[0], mid[1] - eye[1], mid[2] - eye[2]);
    });
    glass.sort(function (p, q) { return q.eyeDistance - p.eyeDistance; });
    return { m: listed.m, solid: pack(solid, true), glass: pack(glass, false) };
  }

  // Packs an explicit list of blocks the same way (the drill-down's
  // stages, which pick their own blocks): each cell has ring, seg, slab
  // and bounds r0..z1, and optionally filled and total. A cell with a
  // `stats` key (the stage API's {systems, expected_systems, stars,
  // mean_age_gy, luminosity_sol}, or null) is colored and made
  // translucent by its sectors' stats (MAP.128/129: statsColor,
  // statsOpacity, ranked against the other cells) instead of lifted
  // toward solid by its filled share. Density comes from the shape. `dim` (optional) is a test: cells it accepts are drawn
  // at a fifth of their opacity (the "Charted only" toggle). Returns
  // {solid, glass} as build does, plus `cells`: [solid cells, glass
  // cells], in the order each part's `owners` index them.
  function buildCells(cells, eye, dim) {
    const solid = [];
    const glass = [];
    ranges = statsRanges(cells);
    cells.forEach(function (cell) {
      if (shape && cell.density === undefined) {
        const sampled = boundsDensity(cell, shape);
        cell.density = sampled.density;
        cell.meanDensity = sampled.mean;
      }
      cell.opacity = "stats" in cell ? statsOpacity(cell, ranges) : blockOpacity(cell);
      if (dim && dim(cell)) {
        cell.opacity *= 0.2;
      }
      (cell.opacity >= 1 ? solid : glass).push(cell);
    });
    eye = eye || [0, 0, 0];
    glass.forEach(function (cell) {
      const mid = cellCoordinates(cell).cartesian;
      cell.eyeDistance = Math.hypot(mid[0] - eye[0], mid[1] - eye[1], mid[2] - eye[2]);
    });
    glass.sort(function (p, q) { return q.eyeDistance - p.eyeDistance; });
    return { solid: pack(solid, true), glass: pack(glass, false), cells: [solid, glass] };
  }

  return { setFilled: setFilled, build: build, buildCells: buildCells };
}

// Every typed array in a built scene, to hand over without copying.
export function transferablesOf(result) {
  const buffers = [];
  [result.solid, result.glass].forEach(function (part) {
    ["positions", "centers", "colors", "uvs", "alphas", "fills", "indices", "owners", "cells"].forEach(function (name) {
      buffers.push(part[name].buffer);
    });
  });
  return buffers;
}

// --- Worker ------------------------------------------------------------------
//
// Messages, answered in order:
// - {type: "init", config}: createBlockScene(config);
// - {type: "filled", points}: setFilled(points) (a POINT_STRIDE
//   Float64Array);
// - {type: "build", id, view}: replies {type: "built", id, result}, or
//   {type: "failed", id, message}.

if (IN_WORKER) {
  let scene = null;
  const handle = function (message) {
    if (message.type === "init") {
      scene = createBlockScene(message.config);
    } else if (message.type === "filled") {
      scene.setFilled(message.points);
    } else if (message.type === "build") {
      let result;
      try {
        result = scene.build(message.view);
      } catch (err) {
        self.postMessage({ type: "failed", id: message.id, message: String(err && err.message || err) });
        return;
      }
      self.postMessage({ type: "built", id: message.id, result: result }, transferablesOf(result));
    }
  };
  self.onmessage = function (event) {
    handle(event.data);
  };
  waiting.splice(0).forEach(handle);
}
