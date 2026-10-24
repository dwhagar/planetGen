// html/static/galaxyblocks.js
//
// The Galaxy Map's block scene (galaxymap3d.js): for one stage's cells, their
// colors and opacities, and their geometry, packed as typed arrays ready for
// the GPU. Like galaxyprisms.js it has no three.js import, so it runs under
// plain node for tests (src/tests/test_galaxyprisms.py).
//
// Colors come in as linear [r, g, b] (the page reads them off
// THREE.Color, which keeps linear components), and every blend here is
// the same component-wise lerp THREE.Color.lerp does. The fade toward the
// view's rim and the near cut are left to the shader, so a built scene
// stays right while the camera zooms and turns within it.

const VERSION_QUERY = new URL(import.meta.url).search;

const {
  boundsDensity, buildPrismGeometry, cellCoordinates,
} = await import(`./galaxyprisms.js${VERSION_QUERY}`);
const { scientificText, threeFigures } = await import(`./numberformat.js${VERSION_QUERY}`);

// Unfilled space: this opaque at the sparsest drawn density, rising to
// BLOCK_OPACITY_DENSE at the densest (90% to 70% see-through). Boss
// (MAP.37): "everything not yet filled should be more transparent by a
// lot", so generated systems stand out; the density shape still reads,
// only faintly.
export const BLOCK_OPACITY_SPARSE = 0.1;
export const BLOCK_OPACITY_DENSE = 0.3;
// Boss (2026-10-08 23:51Z): "never under any frame be more opaque than 50%,
// I always want to see what's in a sector or block even if it's completely
// filled." No block or sector fill is drawn more opaque than this.
export const MAX_FILL_OPACITY = 0.5;
// Any filled sector lifts a block at least this far from its unfilled
// opacity toward solid, however tiny its share: at large m one filled
// sector among 531,441 must still show, and well clear of the unfilled
// blocks around it (MAP.37). The rest of the way follows the
// log of the filled count against the log of the block's total, so a
// block reaches MAX_FILL_OPACITY only when every sector in it is filled.
export const FILLED_MIN_STEP = 0.6;
// MAP.128/129: a drill-down block colored by what its generated sectors
// hold (each cell's `stats`, from the stage API's sector_stats sums).
// Boss (2026-10-08 01:59Z): density sets the transparency (denser is more
// solid, never fully opaque), the stars' mean age sets the hue, and the
// summed luminosity sets the brightness. The color is the color of the
// translucent fill.
// A generated sector with stars is this opaque with next to nothing in it
// and FILLED_OPACITY_DENSE (MAX_FILL_OPACITY) once it holds
// FILLED_FULL_SYSTEMS systems or more, on a log scale: the opacity follows
// what the sector really holds, not how it compares with the others in
// view, so a sector with one star is faint however alone it is on screen
// (Boss, 2026-10-08 23:51Z). One with no stars lifts the unfilled opacity
// by FILLED_EMPTY_STEP and is drawn in the unfilled color a shade more
// saturated.
export const FILLED_OPACITY_SPARSE = 0.12;
export const FILLED_OPACITY_DENSE = MAX_FILL_OPACITY;
export const FILLED_FULL_SYSTEMS = 600;
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
  return base + (MAX_FILL_OPACITY - base) * blockFillStep(cell);
}

// MAP.128: how the cells in view compare: the log of the systems per
// generated sector (the density) and of the summed luminosity, each as a
// {min, max, mean, count} over the cells that hold stars. `statsColor`
// ranks a cell against these (`rangeShare`).
export function statsRanges(cells) {
  const fresh = function () { return { min: Infinity, max: -Infinity, mean: 0, count: 0 }; };
  const ranges = { density: fresh(), luminosity: fresh(), age: fresh(), stars: fresh(), habitable: fresh() };
  cells.forEach(function (cell) {
    const stats = cell.stats;
    if (!stats || !(stats.stars > 0) || !(cell.filled > 0)) {
      return;
    }
    const density = Math.log1p(stats.systems / cell.filled);
    const luminosity = Math.log10(Math.max(stats.luminosity_sol, 1e-9));
    const add = function (range, value) {
      range.min = Math.min(range.min, value);
      range.max = Math.max(range.max, value);
      range.mean += value;
      range.count += 1;
    };
    add(ranges.density, density);
    add(ranges.luminosity, luminosity);
    add(ranges.age, stats.mean_age_gy || 0);
    add(ranges.stars, stats.stars / cell.filled);
    if (stats.habitable_worlds > 0) {
      add(ranges.habitable, COLOR_MODES.habitable.value(stats, cell.filled));
    }
  });
  Object.keys(ranges).forEach(function (name) {
    const range = ranges[name];
    range.mean = range.count ? range.mean / range.count : 0;
  });
  return ranges;
}

// MAP.143: the steps of the Habitable worlds colour and their hue (sRGB, one hue, light rising with the count so
// it still reads in greyscale and against the dark map).
const HABITABLE_BINS = 5;
const HABITABLE_RAMP = [[0.10, 0.30, 0.18], [0.18, 0.49, 0.29], [0.31, 0.68, 0.37], [0.56, 0.83, 0.42], [0.83, 0.96, 0.54]];

// MAP.131: what the Galaxy Map's fill can show. "default" is the age hue,
// density opacity and luminosity brightness above; the others color every
// cell by one statistic on a single ramp. `value` is the statistic of a
// cell holding stars, `range` names its entry in statsRanges, `text` writes
// a value in the legend.
export const COLOR_MODES = {
  default: { label: "Age, density and luminosity" },
  density: {
    label: "Density", range: "density", unit: "systems per charted sector",
    value: function (stats, filled) { return Math.log1p(stats.systems / filled); },
    text: function (v) { return formatLegend(Math.expm1(v)); },
  },
  age: {
    label: "Mean age", range: "age", unit: "Gyr",
    value: function (stats) { return stats.mean_age_gy || 0; },
    text: function (v) { return formatLegend(v); },
  },
  luminosity: {
    label: "Luminosity", range: "luminosity", unit: "L☉ in the cell",
    value: function (stats) { return Math.log10(Math.max(stats.luminosity_sol, 1e-9)); },
    text: function (v) { return formatLegend(Math.pow(10, v)); },
  },
  stars: {
    label: "Star count", range: "stars", unit: "stars per charted sector",
    value: function (stats, filled) { return stats.stars / filled; },
    text: function (v) { return formatLegend(v); },
  },
  // MAP.143: the habitable planets and moons (overall PHI-6 Blue or Green) a cell's sectors hold, per charted sector so
  // sectors and blocks of any size compare, on a log scale in HABITABLE_BINS steps of one hue. A cell holding none
  // is drawn as unfilled space (`empty`), and the range runs over the cells that do hold some.
  habitable: {
    label: "Habitable worlds", range: "habitable", unit: "habitable worlds per charted sector", bins: HABITABLE_BINS,
    value: function (stats, filled) { return Math.log10(stats.habitable_worlds / filled); },
    text: function (v) { return formatLegend(Math.pow(10, v)); },
    empty: function (stats) { return !(stats.habitable_worlds > 0); },
  },
};

function formatLegend(value) {
  if (!isFinite(value)) {
    return "0";
  }
  if (value !== 0 && (Math.abs(value) >= 1e5 || Math.abs(value) < 0.01)) {
    return scientificText(value);
  }
  return threeFigures(value);
}

// The ramp of a single statistic, low to high (sRGB).
const STAT_RAMP = [[0.18, 0.10, 0.45], [0.10, 0.55, 0.60], [0.95, 0.85, 0.25]];

// The step (0 to HABITABLE_BINS - 1) a share of a binned scale falls in.
export function binOf(share, bins) {
  return Math.max(0, Math.min(bins - 1, Math.floor(share * bins)));
}

export function habitableBinColor(bin) {
  return HABITABLE_RAMP[Math.max(0, Math.min(HABITABLE_BINS - 1, bin))];
}

export function statRampColor(share) {
  const t = Math.max(0, Math.min(1, share));
  return t < 0.5 ? lerp3(STAT_RAMP[0], STAT_RAMP[1], t * 2) : lerp3(STAT_RAMP[1], STAT_RAMP[2], (t - 0.5) * 2);
}

// The legend of a color mode over the cells' ranges, or null for the
// default and when no cell holds stars: {mode, label, unit, lowText,
// highText, stops (sRGB hex, low to high)}.
export function legendFor(mode, ranges) {
  const info = COLOR_MODES[mode];
  if (!info || !info.range) {
    return null;
  }
  const range = ranges[info.range];
  if (!(range.max >= range.min)) {
    return null;
  }
  // The scale's own ends (padded around the mean for a few cells); a count
  // or an age never reads below 0.
  const ends = scaleEnds(range, info.range);
  const low = info.range === "luminosity" || info.range === "habitable" ? ends[0] : Math.max(0, ends[0]);
  const hex = function (color) {
    return "#" + color.map(function (v) { return ("0" + Math.round(v * 255).toString(16)).slice(-2); }).join("");
  };
  // A binned mode's legend is its steps side by side (hard stops), a continuous one's the ramp.
  const stops = info.bins ? [].concat.apply([], HABITABLE_RAMP.map(function (color, n) {
    return [hex(color) + " " + (100 * n / info.bins) + "%", hex(color) + " " + (100 * (n + 1) / info.bins) + "%"];
  })) : [0, 0.25, 0.5, 0.75, 1].map(function (t) { return hex(statRampColor(t)); });
  return { mode: mode, label: info.label, unit: info.unit, lowText: info.text(low), highText: info.text(ends[1]), stops: stops };
}

// Boss (2026-10-08 23:52Z): with only one or two stars to compare, don't
// stretch the scale from the lowest to the highest of them. Up to
// FEW_ITEMS of them the scale instead puts their mean in the middle and
// pads both ends, at least SCALE_MIN_HALF_SPAN each way (in the
// statistic's own units: logs for density and luminosity), so each item
// sits about where it belongs around the average.
export const FEW_ITEMS = 3;
export const SCALE_MIN_HALF_SPAN = { density: 1, luminosity: 1, age: 3, stars: 2, habitable: 1 };

// The two ends of a statistic's scale, [low, high], from `range` (the
// statsRanges entry called `name`).
export function scaleEnds(range, name) {
  if (!(range.max >= range.min)) {
    return [0, 1];
  }
  if (range.count > FEW_ITEMS) {
    return [range.min, range.max];
  }
  const half = Math.max(range.max - range.mean, range.mean - range.min, SCALE_MIN_HALF_SPAN[name]);
  return [range.mean - half, range.mean + half];
}

function rangeShare(value, range, name) {
  const ends = scaleEnds(range, name);
  if (!(ends[1] > ends[0])) {
    return 0.5;
  }
  return Math.max(0, Math.min(1, (value - ends[0]) / (ends[1] - ends[0])));
}

// MAP.128/129: a block's opacity from its sectors. Unfilled space keeps
// its own (blockOpacity); a filled block moves toward the opacity its
// stars' density gives it (FILLED_OPACITY_SPARSE up to
// FILLED_OPACITY_DENSE, never more than MAX_FILL_OPACITY) by
// blockFillStep, so one filled sector among many still shows. The density
// is what the block holds (systems per charted sector against
// FILLED_FULL_SYSTEMS, on a log scale), not a rank among the cells in
// view. A filled block with no stars only lifts the unfilled opacity by
// FILLED_EMPTY_STEP.
export function unfilledOpacity(cell) {
  return BLOCK_OPACITY_SPARSE + (BLOCK_OPACITY_DENSE - BLOCK_OPACITY_SPARSE) * prismIntensity(cell.density || 0);
}

export function statsOpacity(cell) {
  const base = unfilledOpacity(cell);
  if (!(cell.filled > 0)) {
    return base;
  }
  const stats = cell.stats;
  if (!stats || !(stats.stars > 0)) {
    return Math.min(FILLED_OPACITY_DENSE, base + FILLED_EMPTY_STEP);
  }
  const held = Math.min(1, Math.log1p(stats.systems / cell.filled) / Math.log1p(FILLED_FULL_SYSTEMS));
  const own = FILLED_OPACITY_SPARSE + (FILLED_OPACITY_DENSE - FILLED_OPACITY_SPARSE) * held;
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
export function statsColor(cell, unfilled, ranges, mode) {
  if (!(cell.filled > 0)) {
    return unfilled;
  }
  const stats = cell.stats;
  if (!stats || !(stats.stars > 0)) {
    return saturate(unfilled, FILLED_EMPTY_SATURATION);
  }
  const single = mode && COLOR_MODES[mode] && COLOR_MODES[mode].range ? COLOR_MODES[mode] : null;
  if (single) {
    if (single.empty && single.empty(stats)) {
      return unfilled;
    }
    const share = rangeShare(single.value(stats, cell.filled), ranges[single.range], single.range);
    return (single.bins ? habitableBinColor(binOf(share, single.bins)) : statRampColor(share)).map(srgbToLinear);
  }
  const luminosity = rangeShare(Math.log10(Math.max(stats.luminosity_sol, 1e-9)), ranges.luminosity, "luminosity");
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
// - shape (the density shape, or null);
// - palette: linear [r, g, b] for dim, accent and hot (the density ramp)
//   and placedLow, placedHigh (a generated sector's own color).
// Returns {buildCells(cells, eye, dim, colorMode, stageRanges), rangesOf(cells), legendOf(mode, ranges)};
// see each below.
export function createBlockScene(config) {
  const shape = config.shape || null;
  const palette = config.palette;
  const armAmplitude = galaxyArmAmplitude(shape);

  function prismColor(t) {
    return t < 0.6 ? lerp3(palette.dim, palette.accent, t / 0.6) : lerp3(palette.accent, palette.hot, (t - 0.6) / 0.4);
  }

  // One mesh's worth of cells, packed: positions, centers (Float32),
  // colors (Uint16, normalized: dark linear colors need more than 8 bits),
  // uvs, alphas and fills (Uint8, normalized), indices, and each vertex's
  // cell (owners, for picking). Normals only light the colors here, so
  // they aren't kept.
  let ranges = statsRanges([]);
  let colorBy = "default";

  function pack(cells) {
    const built = buildPrismGeometry(cells);
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
    cells.forEach(function (cell) {
      middles.push(cellCoordinates(cell).cartesian);
      const ramp = function () { return prismColor(shape ? prismShade(cell, armAmplitude) : 0.5); };
      colorOf.push(statsColor(cell, ramp(), ranges, colorBy));
      alphaOf.push(Math.round(cell.opacity * 255));
      fillOf.push(cell.noFill ? 0 : Math.round(blockFillStep(cell) * 255));
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
      indices: built.indices, owners: built.owners, vertexCount: count,
    };
  }

  // Packs an explicit list of cells (the drill-down's stages, which pick their own blocks): each has
  // its bounds r0..t1 (galaxyprisms.drillBlockBounds) and optionally filled and total. `stageRanges` (optional,
  // from `rangesOf`) ranks the cells against a wider set than `cells`, so every choice of a stage is colored on one
  // scale; `colorMode` (MAP.131, a key of COLOR_MODES) colors the filled cells by one statistic instead; the
  // result's `legend` says what the colors mean. A cell's `stats` (the stage API's {systems, expected_systems,
  // stars, mean_age_gy, luminosity_sol}, or null) colors and makes it translucent by what its sectors hold
  // (MAP.128/129: statsColor, statsOpacity). Density comes from the shape. `dim` (optional) is a test: cells it
  // accepts are drawn at a fifth of their opacity (the "Charted only" toggle). Returns {solid, glass}: cells
  // opaque enough to write depth and the rest (translucent, sorted far to near from `eye`), each packed, plus
  // `cells`: [solid cells, glass cells], in the order each part's `owners` index them.
  function buildCells(cells, eye, dim, colorMode, stageRanges) {
    const solid = [];
    const glass = [];
    ranges = stageRanges || statsRanges(cells);
    colorBy = COLOR_MODES[colorMode] ? colorMode : "default";
    cells.forEach(function (cell) {
      if (shape && cell.density === undefined) {
        const sampled = boundsDensity(cell, shape);
        cell.density = sampled.density;
        cell.meanDensity = sampled.mean;
      }
      cell.opacity = statsOpacity(cell);
      const mode = COLOR_MODES[colorBy];
      cell.noFill = !!(mode.empty && cell.filled > 0 && cell.stats && mode.empty(cell.stats));
      if (cell.noFill) {
        cell.opacity = unfilledOpacity(cell);
      }
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
    return { solid: pack(solid), glass: pack(glass), cells: [solid, glass], legend: legendFor(colorBy, ranges) };
  }

  return { buildCells: buildCells, rangesOf: statsRanges, legendOf: legendFor };
}
