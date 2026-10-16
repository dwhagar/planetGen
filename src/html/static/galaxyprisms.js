// html/static/galaxyprisms.js
//
// The 3D Galaxy Map's sector grid and predicted-density shading
// (galaxymap3d.js), on the spherical shell grid the generator cuts space
// into (planetgen/galaxy/spherical.py, docs/design/spherical-sector-
// coordinates.md). A cell is named (shell, band, wedge): a shell is a
// spherical skin one sector edge deep around the galactic origin (shell 0
// the core ball), a band a latitude ring of the shell, a wedge one of the
// band's equal longitude steps, pole to pole. Each polar cap is one cell.
// This module is the page's copy of that arithmetic: the grid itself, the
// block ladder (planetgen/galaxy/drill.py), the data window
// (planetgen/galaxy/window.py), the galaxy's outline from the density bound
// (planetgen/galaxy/skeleton.py) and how many sectors a block can hold under
// it, and the geometry of a cell: a spherical box with curved inner and
// outer walls, two cones for the band's colatitudes and two flat meridian
// sides. Every function repeats the Python step for step (JavaScript's
// round-half-up included) and tests/test_galaxyprisms.py checks the two
// agree.
//
// Density comes from the galaxy's own analytic model
// (planetgen.galaxy.density.relative_density), evaluated right here from the
// shape parameters the page embeds -- no server round trip, no sampled point
// cloud. Each cell is shaded by its mean density over a few sample points.
//
// No three.js import: this module only does the arithmetic and hands back
// plain typed arrays, so it runs under plain node for tests
// (src/tests/test_galaxyprisms.py) and in the page's Web Worker
// (galaxyblocks.js).

// The drawing floor: cells thinner than this relative density are left out -- only when the
// shape doesn't carry the galaxy's own sector threshold (shape.sector_min_density, the density
// a sector needs to expect one system). With it, a sector is drawn exactly when the
// generator's skeleton would allow it, so at one sector per cell the map's outline is the
// galaxy's.
export var PRISM_MIN_DENSITY = 0.02;
// Curved walls get one segment per this many radians of arc.
var ARC_SEGMENT_RAD = 0.12;
var MAX_ARC_SEGMENTS = 24;
// Sample points per cell: radial x colatitude x longitude.
var SAMPLES_R = 2;
var SAMPLES_LON = 2;
var SAMPLES_COLAT = 3;

function sechSquared(x) {
  var ax = Math.abs(x);
  if (ax > 700) {
    return 0;
  }
  var e = Math.exp(-2 * ax);
  return (4 * e) / ((1 + e) * (1 + e));
}

// Port of galaxyDensity.relative_density, same field names as the
// GalaxyShape namedtuple plus galaxyDensity.model_terms (the thick disk
// and the bar), which the page embeds alongside them.
export function relativeDensity(x, y, z, shape) {
  return densityParts(x, y, z, shape)[0];
}

// A disk's vertical profile: sech^2(z / 2h), exponential scale height h
// (galaxyDensity._vertical).
function vertical(z, h) {
  return sechSquared(z / (2 * h));
}

// The galactic core (galaxyDensity._core): a Gaussian at the centre, rSq the squared distance from it.
function core(rSq, shape) {
  if (!shape.core_amplitude) {
    return 0;
  }
  return shape.core_amplitude * Math.exp(-0.5 * rSq / (shape.core_scale_radius_pc * shape.core_scale_radius_pc));
}

// The boxy bar bulge (galaxyDensity._bulge), plus the core.
function bulge(x, y, z, shape) {
  var along = (x * shape.bar_cos + y * shape.bar_sin) / shape.bulge_scale_radius_pc;
  var across = (y * shape.bar_cos - x * shape.bar_sin) / shape.bulge_scale_y_pc;
  var up = z / shape.bulge_scale_z_pc;
  var inPlane = along * along + across * across;
  return shape.bulge_amplitude * Math.exp(-0.5 * Math.sqrt(inPlane * inPlane + up * up * up * up))
    + core(x * x + y * y + z * z, shape);
}

// The bulge's maximum over azimuth: along the bar (galaxyDensity.bulge_bound).
function bulgeBound(r, z, shape) {
  var along = r / shape.bulge_scale_radius_pc;
  var up = z / shape.bulge_scale_z_pc;
  return shape.bulge_amplitude * Math.exp(-0.5 * Math.sqrt(along * along * along * along + up * up * up * up))
    + core(r * r + z * z, shape);
}

// The stellar halo (galaxyDensity._halo): a broken power law on a spheroid of axis ratio q, flat inside the core.
function halo(rCyl, z, shape) {
  var m = Math.hypot(rCyl, z / shape.halo_axis_ratio);
  var ratio = Math.max(m, shape.halo_core_pc) / shape.halo_break_pc;
  return shape.halo_amplitude * Math.pow(ratio, -(ratio < 1 ? shape.halo_inner_slope : shape.halo_outer_slope));
}

function thickDisk(rCyl, z, shape) {
  return shape.thick_disk_amplitude * Math.exp(-rCyl / shape.thick_disk_scale_length_pc)
    * vertical(z, shape.thick_disk_scale_height_pc);
}

// [density, the same with the arm factor held at 1]: the second is the
// density's azimuthal mean around the ring as far as the arms go (the
// arm factor averages to 1 around any circle), so density / mean is the
// arm factor's pull.
function densityParts(x, y, z, shape) {
  var rCyl = Math.hypot(x, y);
  var rest = bulge(x, y, z, shape) + thickDisk(rCyl, z, shape) + halo(rCyl, z, shape);
  var disk = shape.arm_level * Math.exp(-rCyl / shape.disk_scale_length_pc) * vertical(z, shape.disk_scale_height_pc);
  var armFactor = 1;
  if (rCyl > 1e-9) {
    var theta = Math.atan2(y, x);
    var thetaArm = shape.spiral_reference_angle_rad
      + Math.log(rCyl / shape.spiral_reference_radius_pc) / Math.tan(shape.pitch_angle_rad);
    armFactor = 1 + shape.arm_amplitude * Math.cos(shape.arm_count * (theta - thetaArm));
  }
  return [shape.k_norm * (rest + disk * armFactor), shape.k_norm * (rest + disk)];
}

// The most any point in the prism could have -- cheap enough to run on
// every candidate, so empty space is skipped before it's sampled. Same
// bound as galaxySkeleton.bound_relative_density_at.
function densityUpperBound(r0, zMinAbs, shape) {
  var disk = shape.arm_level * Math.exp(-r0 / shape.disk_scale_length_pc)
    * vertical(zMinAbs, shape.disk_scale_height_pc)
    * (1 + Math.abs(shape.arm_amplitude));
  return shape.k_norm * (bulgeBound(r0, zMinAbs, shape) + disk + thickDisk(r0, zMinAbs, shape));
}

// --- The sector grid (mirrors planetgen/galaxy/spherical.py) ---------------
//
// Shell k: radius [(k - 1/2) e, (k + 1/2) e), centreline k * e, shell 0 the
// core ball of radius e / 2. A shell from 2 out has two polar caps (colatitude
// below alpha = asin(rho / R), rho the six-edge rim's chord radius) and an
// odd number of bands of equal colatitude step between them; band 0 is the
// equatorial band, positive bands lie north. A band holds 6 * m wedges, m
// with no prime factor above 7, nearest one edge of arc; a cap and the core
// are one cell. Shell 1 is cut by hand. Every point of space falls in
// exactly one cell.

var TWO_PI = 2 * Math.PI;
var CAP_RIM_WEDGES = 6;
var BAND_ARC_EDGES = 1.0;
var CAP_CHORD_RADIUS_EDGES = CAP_RIM_WEDGES / TWO_PI;
// Shells the rim rule cannot cut: shell: [cap colatitude (degrees), wedges per band north to south].
var HAND_CUT_SHELLS = { 1: [37.5, [12]] };
// The designation packs shell, biased band and wedge into 56 bits (object_uid.LAYOUT: 19, 18 and 19 bits).
var DESIGNATION_WEDGE_BITS = 19n;
var DESIGNATION_BAND_BITS = 18n;
var DESIGNATION_BAND_BIAS = 131072;
// The outermost shell the designation's band field can name (spherical.check_reach).
var DESIGNATION_MAX_SHELL = 83037;
var SMOOTH_LIMIT = 1e7;

// Every number up to SMOOTH_LIMIT whose prime factors are only 2, 3, 5 and 7, ascending.
var SMOOTH_MULTIPLIERS = (function () {
  var found = [];
  for (var a = 1; a <= SMOOTH_LIMIT; a *= 2) {
    for (var b = a; b <= SMOOTH_LIMIT; b *= 3) {
      for (var c = b; c <= SMOOTH_LIMIT; c *= 5) {
        for (var d = c; d <= SMOOTH_LIMIT; d *= 7) found.push(d);
      }
    }
  }
  return found.sort(function (p, q) { return p - q; });
}());

// bisect_left: the first index whose value is not below x.
function bisectLeft(values, x) {
  var lo = 0;
  var hi = values.length;
  while (lo < hi) {
    var mid = (lo + hi) >> 1;
    if (values[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

// bisect_right: the first index whose value is above x.
function bisectRight(values, x) {
  var lo = 0;
  var hi = values.length;
  while (lo < hi) {
    var mid = (lo + hi) >> 1;
    if (x < values[mid]) hi = mid;
    else lo = mid + 1;
  }
  return lo;
}

// The 7-smooth number nearest x in ratio (a tie goes to the lower), at least 1.
function nearestSmooth(x) {
  if (x <= 1.0) return 1;
  var high = bisectLeft(SMOOTH_MULTIPLIERS, x);
  if (high >= SMOOTH_MULTIPLIERS.length) throw new Error("wedge multiplier " + x + " is past the table");
  var upper = SMOOTH_MULTIPLIERS[high];
  if (upper === x) return upper;
  var lower = SMOOTH_MULTIPLIERS[high - 1];
  return x * x <= lower * upper ? lower : upper;
}

function makeLayout(shell, cap, bandCount, table) {
  return {
    shell: shell, cap: cap, bandCount: bandCount,
    bandStep: bandCount ? (Math.PI - 2 * cap) / bandCount : 0.0,
    highBand: Math.floor((bandCount - 1) / 2),
    table: table,
  };
}

var LAYOUT_CACHE_MAX = 20000;
var layouts = new Map();

// How one shell is cut: {shell, cap (alpha), bandCount (odd; 0 for the core),
// bandStep, highBand, table (wedges per band for a hand-cut shell)}.
export function shellLayout(shell) {
  var layout = layouts.get(shell);
  if (layout) return layout;
  if (!(shell >= 0) || Math.floor(shell) !== shell) throw new Error("shell must be a whole number >= 0, got " + shell);
  if (shell === 0) {
    layout = makeLayout(0, 0.0, 0, null);
  } else if (HAND_CUT_SHELLS[shell]) {
    var cut = HAND_CUT_SHELLS[shell];
    layout = makeLayout(shell, cut[0] * (Math.PI / 180), cut[1].length, cut[1]);
  } else {
    var alpha = Math.asin(CAP_CHORD_RADIUS_EDGES / shell);
    var n = 2 * Math.floor(shell * (Math.PI - 2 * alpha) / BAND_ARC_EDGES / 2) + 1;
    layout = makeLayout(shell, alpha, n, null);
  }
  if (layouts.size >= LAYOUT_CACHE_MAX) layouts.clear();
  layouts.set(shell, layout);
  return layout;
}

// Colatitude of the i-th band edge (0 .. bandCount), north to south, symmetric about the equator.
function layoutEdge(layout, i) {
  var n = layout.bandCount;
  if (i <= Math.floor(n / 2)) return layout.cap + i * layout.bandStep;
  return Math.PI - (layout.cap + (n - i) * layout.bandStep);
}

// [lowest, highest] band of a shell, caps included; [0, 0] for the core.
export function bandRange(shell) {
  if (shell === 0) return [0, 0];
  var layout = shellLayout(shell);
  return [-(layout.highBand + 1), layout.highBand + 1];
}

export function isCap(shell, band) {
  return shell > 0 && Math.abs(band) === shellLayout(shell).highBand + 1;
}

function checkBand(shell, band) {
  var range = bandRange(shell);
  if (!(range[0] <= band && band <= range[1])) {
    throw new Error("band " + band + " out of range for shell " + shell + " (holds " + range[0] + ".." + range[1] + ")");
  }
}

// Wedges of the ordinary band `i` places south of the shell's top: the rule is worked out for the
// northern one of a mirrored pair, so the two always agree.
function wedgesAt(layout, i) {
  if (layout.table) return layout.table[i];
  i = Math.min(i, layout.bandCount - 1 - i);
  var centre = (layoutEdge(layout, i) + layoutEdge(layout, i + 1)) / 2;
  var arc = TWO_PI * layout.shell * Math.sin(centre) / BAND_ARC_EDGES;
  return CAP_RIM_WEDGES * nearestSmooth(arc / CAP_RIM_WEDGES);
}

// How many wedges a band holds; a cap and the core are one cell.
export function wedgeCount(shell, band) {
  checkBand(shell, band);
  if (shell === 0 || isCap(shell, band)) return 1;
  var layout = shellLayout(shell);
  return wedgesAt(layout, layout.highBand - band);
}

// [inner, outer] radius of shell `shell`, parsecs: the core is [0, e / 2].
export function shellBoundsPc(shell, edgePc) {
  if (shell === 0) return [0.0, edgePc / 2];
  return [(shell - 0.5) * edgePc, (shell + 0.5) * edgePc];
}

// [north, south] colatitude of a band, radians from galactic north: a cap is [0, alpha] or
// [pi - alpha, pi], the core [0, pi].
export function bandColatitudeBounds(shell, band) {
  checkBand(shell, band);
  if (shell === 0) return [0.0, Math.PI];
  var layout = shellLayout(shell);
  if (isCap(shell, band)) {
    return band > 0 ? [0.0, layout.cap] : [Math.PI - layout.cap, Math.PI];
  }
  var i = layout.highBand - band;
  return [layoutEdge(layout, i), layoutEdge(layout, i + 1)];
}

// [start, end] longitude of one wedge, radians from +X counterclockwise.
export function wedgeLongitudeBounds(shell, band, wedge) {
  var n = checkWedge(shell, band, wedge);
  var step = TWO_PI / n;
  return [wedge * step, (wedge + 1) * step];
}

function checkWedge(shell, band, wedge) {
  var n = wedgeCount(shell, band);
  if (!(wedge >= 0 && wedge < n) || Math.floor(wedge) !== wedge) {
    throw new Error("wedge " + wedge + " out of range for shell " + shell + " band " + band + " (holds " + n + " wedges)");
  }
  return n;
}

// [cylindrical radius, z] of a band's centre line, parsecs: where every wedge of the band has its centre.
export function bandCentrePc(shell, band, edgePc) {
  checkBand(shell, band);
  if (shell === 0) return [0.0, 0.0];
  var r = shell * edgePc;
  if (isCap(shell, band)) return [0.0, band > 0 ? r : -r];
  if (band === 0) return [r, 0.0];
  var bounds = bandColatitudeBounds(shell, band);
  var colat = (bounds[0] + bounds[1]) / 2;
  return [r * Math.sin(colat), r * Math.cos(colat)];
}

// A sector's centre [x, y, z], parsecs: on its shell's centreline at the middle of its band's
// colatitude and its wedge's longitude; a cap's is on the axis, the core's the origin.
export function sectorPositionPc(shell, band, wedge, edgePc) {
  checkWedge(shell, band, wedge);
  var centre = bandCentrePc(shell, band, edgePc);
  if (centre[0] === 0.0) return [0.0, 0.0, centre[1]];
  var lon = wedgeLongitudeBounds(shell, band, wedge);
  var mid = (lon[0] + lon[1]) / 2;
  return [centre[0] * Math.cos(mid), centre[0] * Math.sin(mid), centre[1]];
}

// The shell whose [inner, outer) range holds a radius.
export function shellIndexAt(radiusPc, edgePc) {
  if (radiusPc < edgePc / 2) return 0;
  var guess = Math.floor(radiusPc / edgePc + 0.5);
  var bounds = shellBoundsPc(guess, edgePc);
  if (radiusPc < bounds[0]) return guess - 1;
  if (radiusPc >= bounds[1]) return guess + 1;
  return guess;
}

function layoutBandAt(layout, colatitude) {
  var alpha = layout.cap;
  if (colatitude < alpha) return layout.highBand + 1;
  if (colatitude >= Math.PI - alpha) return -(layout.highBand + 1);
  var n = layout.bandCount;
  var i = Math.min(n - 1, Math.max(0, Math.trunc((colatitude - alpha) / layout.bandStep)));
  if (colatitude < layoutEdge(layout, i) && i > 0) i -= 1;
  else if (colatitude >= layoutEdge(layout, i + 1) && i < n - 1) i += 1;
  return layout.highBand - i;
}

// The band of shell `shell` holding a colatitude (radians from galactic north), a cap included.
export function bandAt(shell, colatitude) {
  if (shell === 0) return 0;
  return layoutBandAt(shellLayout(shell), colatitude);
}

// Python's atan2(y, x) % (2 * pi), to the bit: a tiny negative angle becomes 2 * pi itself and
// lands in the last wedge.
function bearing(y, x) {
  var theta = Math.atan2(y, x);
  return theta < 0 ? theta + TWO_PI : theta + 0;
}

// The wedge of a band of `n` holding bearing `theta` in [0, 2 * pi]: each division-based guess is
// checked against the wedge's own edges and moved by one if it rounded across one.
function wedgeOfBearing(theta, n) {
  if (n === 1) return 0;
  var wedge = Math.min(n - 1, Math.floor((theta * n) / TWO_PI));
  var step = TWO_PI / n;
  if (wedge > 0 && theta < wedge * step) wedge -= 1;
  else if (wedge < n - 1 && theta >= (wedge + 1) * step) wedge += 1;
  return wedge;
}

// The (shell, band, wedge) cell holding galaxy-frame point (x, y, z).
export function sectorAddressAt(x, y, z, edgePc) {
  var radius = Math.sqrt(x * x + y * y + z * z);
  var shell = shellIndexAt(radius, edgePc);
  if (shell === 0) return { shell: 0, band: 0, wedge: 0 };
  var layout = shellLayout(shell);
  var band = layoutBandAt(layout, Math.atan2(Math.sqrt(x * x + y * y), z));
  var n = wedgeCount(shell, band);
  return { shell: shell, band: band, wedge: wedgeOfBearing(bearing(y, x), n) };
}

// One sector's cell as bounds: {r0, r1, c0, c1, t0, t1} -- radius, colatitude and longitude, the
// colatitudes and longitudes in radians.
export function sectorBounds(shell, band, wedge, edgePc) {
  var radii = shellBoundsPc(shell, edgePc);
  var colat = bandColatitudeBounds(shell, band);
  var lon = wedgeLongitudeBounds(shell, band, wedge);
  return { r0: radii[0], r1: radii[1], c0: colat[0], c1: colat[1], t0: lon[0], t1: lon[1] };
}

// The 8 corners of a cell given as bounds, galaxy frame. Index is 4*r_bit + 2*lat_bit + lon_bit, as
// spherical.sector_cell_vertices_pc: corners i and i ^ 1, i ^ 2, i ^ 4 share an edge.
export function cellVertices(b) {
  var out = [];
  [b.r0, b.r1].forEach(function (r) {
    [b.c0, b.c1].forEach(function (c) {
      [b.t0, b.t1].forEach(function (t) {
        var s = r * Math.sin(c);
        out.push([s * Math.cos(t), s * Math.sin(t), r * Math.cos(c)]);
      });
    });
  });
  return out;
}

// Whether a cell's longitude runs the whole way round (a cap or the core): it has no meridian sides.
export function isFullTurn(b) {
  return b.t1 - b.t0 >= TWO_PI - 1e-9;
}

// A cell's middle: halfway between its radii, colatitudes and longitudes; a cell all the way round
// (a cap, the core) sits on the axis.
export function cellCenter(b) {
  var r = (b.r0 + b.r1) / 2;
  var c = (b.c0 + b.c1) / 2;
  if (isFullTurn(b)) return [0, 0, r * Math.cos(c)];
  var t = (b.t0 + b.t1) / 2;
  var s = r * Math.sin(c);
  return [s * Math.cos(t), s * Math.sin(t), r * Math.cos(c)];
}

// A cell's center in the three coordinate systems the map shows: Cartesian (x, y, z), cylindrical
// (R, theta, z) and spherical (r, theta, polar angle from galactic north).
export function cellCoordinates(b) {
  var p = cellCenter(b);
  var t = isFullTurn(b) ? 0 : (b.t0 + b.t1) / 2;
  var rho = Math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
  return {
    cartesian: p,
    cylindrical: [Math.hypot(p[0], p[1]), t, p[2]],
    spherical: [rho, t, rho > 0 ? Math.acos(p[2] / rho) : 0],
  };
}

// Whether galaxy-frame point (x, y, z) lies in a cell given as bounds (lower faces inclusive, upper
// exclusive, as the grid holds a point).
export function boundsContain(b, x, y, z) {
  var r = Math.sqrt(x * x + y * y + z * z);
  if (r < b.r0 || r >= b.r1) return false;
  if (r === 0) return true;
  var colat = Math.atan2(Math.sqrt(x * x + y * y), z);
  if (colat < b.c0 || (colat >= b.c1 && b.c1 < Math.PI)) return false;
  if (isFullTurn(b) || (x === 0 && y === 0)) return true;
  var t = bearing(y, x);
  return t >= b.t0 && t < b.t1;
}

// provisional_sector_designation's hex code: shell << 28 | (band + 2048) << 16 | wedge (BigInt, since it runs
// past 32 bits).
export function sectorDesignation(shell, band, wedge) {
  var packed = (((BigInt(shell) << DESIGNATION_BAND_BITS) | BigInt(band + DESIGNATION_BAND_BIAS)) << DESIGNATION_WEDGE_BITS)
    | BigInt(wedge);
  return packed.toString(16).toUpperCase();
}

// sectorDesignation's inverse, or null for anything that isn't one (not hex, past 56 bits, or not a cell of the grid).
export function parseSectorDesignation(text) {
  var hex = String(text == null ? "" : text).trim();
  if (!/^[0-9A-Fa-f]{1,14}$/.test(hex)) return null;
  var packed = BigInt("0x" + hex);
  var wedge = Number(packed & ((1n << DESIGNATION_WEDGE_BITS) - 1n));
  var band = Number((packed >> DESIGNATION_WEDGE_BITS) & ((1n << DESIGNATION_BAND_BITS) - 1n)) - DESIGNATION_BAND_BIAS;
  var shell = Number(packed >> (DESIGNATION_WEDGE_BITS + DESIGNATION_BAND_BITS));
  if (!cellExists(shell, band, wedge)) return null;
  return { shell: shell, band: band, wedge: wedge };
}

// Whether (shell, band, wedge) is a cell of the grid.
export function cellExists(shell, band, wedge) {
  if (!Number.isSafeInteger(shell) || !Number.isSafeInteger(band) || !Number.isSafeInteger(wedge)) return false;
  if (shell < 0 || shell > DESIGNATION_MAX_SHELL || wedge < 0) return false;
  var range = bandRange(shell);
  if (band < range[0] || band > range[1]) return false;
  return wedge < wedgeCount(shell, band);
}

// --- Drill-down ladder ------------------------------------------------------
//
// The drill-down's nested blocks (docs/design/spherical-block-ladder.md):
// 243, 81, 27, 9 and 3 shells deep and a sector (1), each block wholly inside
// one block of the next size up by counts. A block is {m, shell, band, wedge}:
// block shell `shell` holds the sector shells shell*m .. shell*m + m - 1 (its
// centreline radius is shell*m + (m - 1)/2 edges), `band` is a signed row of
// blocks (0 the equatorial row) and `wedge` one of the row's equal longitude
// steps. At m = 1 it is a sector as spherical.py names it. A sector joins the
// level-3 block holding its centre. planetgen/galaxy/drill.py is the same
// rules for the server, and tests/test_galaxyprisms.py checks the two agree.

export var DRILL_LEVELS = [243, 81, 27, 9, 3, 1];

// How many child block shells deep one block is, at every level.
var DRILL_STEP = 3;
var DRILL_TOP = 243;

function jsRound(x) {
  return Math.floor(x + 0.5);
}

// The odd integer nearest x, a tie going up, and at least 1.
function oddRound(x) {
  return Math.max(1, 2 * Math.floor(x / 2) + 1);
}

// Centreline radius of block shell `shell` at size m, in edges.
export function blockRadius(m, shell) {
  return shell * m + (m - 1) / 2;
}

var LEVEL_INDEX = { 243: 0, 81: 1, 27: 2, 9: 3, 3: 4 };

function ladderKey(m, shell, band) {
  return ((LEVEL_INDEX[m] * 8192 + shell) * 16384) + band + 8192;
}

var bandCounts = new Map();
var wedgeCounts = new Map();
var LADDER_CACHE_MAX = 600000;

function checkLevel(m) {
  if (LEVEL_INDEX[m] === undefined) throw new Error("block size must be one of 243, 81, 27, 9, 3, got " + m);
}

// Rows of blocks in block shell `shell` at size m (always odd): 243 takes the odd integer nearest
// pi * R / m; the others the parent block shell's count times the odd integer nearest that over the parent's.
export function drillBandCount(m, shell) {
  checkLevel(m);
  var key = LEVEL_INDEX[m] * 8192 + shell;
  var cached = bandCounts.get(key);
  if (cached !== undefined) return cached;
  var target = Math.PI * blockRadius(m, shell) / m;
  var count;
  if (m === DRILL_TOP) {
    count = oddRound(target);
  } else {
    var parent = drillBandCount(m * DRILL_STEP, Math.floor(shell / DRILL_STEP));
    count = parent * oddRound(target / parent);
  }
  bandCounts.set(key, count);
  return count;
}

// The signed row of the parent block shell holding row `band`.
function parentBandOf(m, shell, band) {
  var count = drillBandCount(m, shell);
  var parentCount = drillBandCount(m * DRILL_STEP, Math.floor(shell / DRILL_STEP));
  var row = (count - 1) / 2 - band;
  var parentRow = Math.floor(row / (count / parentCount));
  return (parentCount - 1) / 2 - parentRow;
}

// The north-to-south index of a signed row.
function rowOf(m, shell, band) {
  var count = drillBandCount(m, shell);
  var row = (count - 1) / 2 - band;
  if (!(row >= 0 && row < count)) {
    throw new Error("band " + band + " is not a row of block shell " + shell + " at size " + m);
  }
  return row;
}

// Wedges in row `band` of block shell `shell` at size m: 243 takes round(2 pi R sin(theta) / m), theta the
// row's centre colatitude; the others the parent row's count times max(1, round(that over the parent's)).
export function drillWedgeCount(m, shell, band) {
  if (m === 1) return wedgeCount(shell, band);
  checkLevel(m);
  if (wedgeCounts.size > LADDER_CACHE_MAX) wedgeCounts.clear();
  var key = ladderKey(m, shell, band);
  var cached = wedgeCounts.get(key);
  if (cached !== undefined) return cached;
  var count = drillBandCount(m, shell);
  var row = rowOf(m, shell, band);
  var theta = (row + 0.5) * Math.PI / count;
  var target = 2 * Math.PI * blockRadius(m, shell) * Math.sin(theta) / m;
  var wedges;
  if (m === DRILL_TOP) {
    wedges = Math.max(1, jsRound(target));
  } else {
    var parent = drillWedgeCount(m * DRILL_STEP, Math.floor(shell / DRILL_STEP), parentBandOf(m, shell, band));
    wedges = parent * Math.max(1, jsRound(target / parent));
  }
  wedgeCounts.set(key, wedges);
  return wedges;
}

// The row (signed) of block shell `shell` at size m holding a colatitude; the poles go to the end rows.
function centreRow(m, shell, colatitude) {
  var count = drillBandCount(m, shell);
  var row = Math.min(count - 1, Math.max(0, Math.floor(colatitude * count / Math.PI)));
  return (count - 1) / 2 - row;
}

// The colatitude a sector's centre lies at: the middle of its band, the plane for band 0, the axis for a
// cap, and pi / 2 for the core.
export function cellColatitude(shell, band) {
  if (shell === 0) return Math.PI / 2;
  if (isCap(shell, band)) return band > 0 ? 0.0 : Math.PI;
  if (band === 0) return Math.PI / 2;
  var bounds = bandColatitudeBounds(shell, band);
  return (bounds[0] + bounds[1]) / 2;
}

// The block wedge holding a sector wedge's centre longitude, in integers; a centre exactly on a meridian
// goes to the higher wedge.
function cellWedgeIn(cellWedges, cellWedge, blockWedges) {
  return Math.floor(((2 * cellWedge + 1) * blockWedges) / (2 * cellWedges));
}

// The block one size up holding `block` ({m, shell, band, wedge}, or a sector at m = 1), or null at 243.
export function drillParent(block) {
  var m = block.m;
  if (m === 1) {
    var parentShell = Math.floor(block.shell / DRILL_STEP);
    var parentBand = centreRow(DRILL_STEP, parentShell, cellColatitude(block.shell, block.band));
    if (block.shell === 0 || isCap(block.shell, block.band)) {
      return { m: DRILL_STEP, shell: parentShell, band: parentBand, wedge: 0 };
    }
    var parentWedges = drillWedgeCount(DRILL_STEP, parentShell, parentBand);
    return {
      m: DRILL_STEP, shell: parentShell, band: parentBand,
      wedge: cellWedgeIn(wedgeCount(block.shell, block.band), block.wedge, parentWedges),
    };
  }
  checkLevel(m);
  if (m === DRILL_TOP) return null;
  var shell = Math.floor(block.shell / DRILL_STEP);
  var band = parentBandOf(m, block.shell, block.band);
  var q = drillWedgeCount(m, block.shell, block.band) / drillWedgeCount(m * DRILL_STEP, shell, band);
  return { m: m * DRILL_STEP, shell: shell, band: band, wedge: Math.floor(block.wedge / q) };
}

// Sector (shell, band, wedge)'s blocks from the top: [243, 81, 27, 9, 3, sector].
export function drillChainOf(shell, band, wedge) {
  var chain = [{ m: 1, shell: shell, band: band, wedge: wedge }];
  while (chain[chain.length - 1].m !== DRILL_TOP) chain.push(drillParent(chain[chain.length - 1]));
  return chain.reverse();
}

// "m.shell.band.wedge", a stage URL's `at`.
export function formatDrillKey(block) {
  return [block.m, block.shell, block.band, block.wedge].join(".");
}

export function sameBlock(a, b) {
  if (!a || !b) return a === b;
  return a.m === b.m && a.shell === b.shell && a.band === b.band && a.wedge === b.wedge;
}

// A block's bounds, {r0, r1, c0, c1, t0, t1}: the radii of its sector shells, its row's colatitudes (equal
// steps of the block shell's rows) and its wedge's longitudes. Block shell 0's inner radius is 0. At m = 1,
// the sector's own cell.
export function drillBlockBounds(block, edgePc) {
  var m = block.m;
  if (m === 1) return sectorBounds(block.shell, block.band, block.wedge, edgePc);
  checkLevel(m);
  var rows = drillBandCount(m, block.shell);
  var row = rowOf(m, block.shell, block.band);
  var step = TWO_PI / drillWedgeCount(m, block.shell, block.band);
  return {
    r0: block.shell === 0 ? 0.0 : (block.shell * m - 0.5) * edgePc,
    r1: ((block.shell + 1) * m - 0.5) * edgePc,
    c0: row * Math.PI / rows, c1: (row + 1) * Math.PI / rows,
    t0: block.wedge * step, t1: (block.wedge + 1) * step,
  };
}

// Whether {m, shell, band, wedge} is a block (or sector) of the ladder.
export function isBlock(block) {
  var m = block.m;
  if (DRILL_LEVELS.indexOf(m) < 0) return false;
  if (!Number.isSafeInteger(block.shell) || !Number.isSafeInteger(block.band) || !Number.isSafeInteger(block.wedge)) return false;
  if (block.shell < 0 || block.wedge < 0) return false;
  if (m === 1) return cellExists(block.shell, block.band, block.wedge);
  if (block.shell >= 8192) return false;
  var rows = drillBandCount(m, block.shell);
  if (Math.abs(block.band) > (rows - 1) / 2) return false;
  return block.wedge < drillWedgeCount(m, block.shell, block.band);
}

// The colatitude of a block row's centre (rows are equal steps from the north pole).
function blockColatitude(m, shell, band) {
  var rows = drillBandCount(m, shell);
  return ((rows - 1) / 2 - band + 0.5) * Math.PI / rows;
}

// A cell's centre in parsecs, [x, y, z]: a sector's own, a block's at its block shell's centreline radius, its
// row's centre colatitude and its wedge's centre longitude (window.cell_position_pc).
export function cellPositionPc(cell, edgePc) {
  if (cell.m === 1) return sectorPositionPc(cell.shell, cell.band, cell.wedge, edgePc);
  var radius = blockRadius(cell.m, cell.shell) * edgePc;
  var colat = blockColatitude(cell.m, cell.shell, cell.band);
  var lon = (cell.wedge + 0.5) * TWO_PI / drillWedgeCount(cell.m, cell.shell, cell.band);
  var plane = radius * Math.sin(colat);
  return [plane * Math.cos(lon), plane * Math.sin(lon), radius * Math.cos(colat)];
}

// The cell of size `m` (a block size, or 1 for a sector) holding galaxy-frame point (x, y, z). A block is found
// by its block shell, row and wedge; a sector by the grid.
export function cellAt(m, x, y, z, edgePc) {
  if (m === 1) {
    var s = sectorAddressAt(x, y, z, edgePc);
    return { m: 1, shell: s.shell, band: s.band, wedge: s.wedge };
  }
  checkLevel(m);
  var radius = Math.sqrt(x * x + y * y + z * z) / edgePc;
  var shell = Math.max(0, Math.floor((radius + 0.5) / m));
  var rows = drillBandCount(m, shell);
  var colat = Math.atan2(Math.sqrt(x * x + y * y), z);
  var row = Math.min(rows - 1, Math.max(0, Math.floor(colat * rows / Math.PI)));
  var band = (rows - 1) / 2 - row;
  var wedges = drillWedgeCount(m, shell, band);
  var wedge = wedgeOfBearing(bearing(y, x), wedges);
  return { m: m, shell: shell, band: band, wedge: wedge };
}

// --- The data window ---------------------------------------------------------
//
// planetgen/galaxy/window.py (MAP.151) for the page: `g` cells across around
// a chosen cell of one level, a ball: every cell of the centre's level whose
// centre lies within g * m / 2 edges of the centre cell's centre, found per
// shell in closed form (a ball meets a sphere in a cap, so the bands to look
// at are those within the cap's angular radius, and in each the wedges within
// a longitude half-width). tests/test_galaxyprisms.py checks the two agree.

export var WINDOW_WIDTHS = [3, 5, 7, 9, 11];
export var DEFAULT_WINDOW = 9;
var TOLERANCE_PC = 1e-9;

var bandCache = new Map();
var BAND_CACHE_MAX = 400;

// The sector bands of a shell north to south with their centre colatitudes: {bands, centres}.
function shellBands(shell) {
  var cached = bandCache.get(shell);
  if (cached) return cached;
  var range = bandRange(shell);
  var bands = [];
  var centres = [];
  for (var b = range[1]; b >= range[0]; b--) {
    bands.push(b);
    centres.push(cellColatitude(shell, b));
  }
  if (bandCache.size >= BAND_CACHE_MAX) bandCache.clear();
  cached = { bands: bands, centres: centres };
  bandCache.set(shell, cached);
  return cached;
}

// The shell ids of level m (block shells, or sector shells for 1) whose centreline radius lies in [low, high] (pc),
// with that radius.
function levelShells(m, lowPc, highPc, edgePc) {
  var out = [];
  if (m === 1) {
    var first = Math.max(0, Math.floor(lowPc / edgePc));
    var last = Math.ceil(highPc / edgePc);
    for (var k = first; k <= last; k++) {
      if (lowPc <= k * edgePc && k * edgePc <= highPc) out.push([k, k * edgePc]);
    }
    return out;
  }
  var firstBlock = Math.max(0, Math.floor((lowPc / edgePc - (m - 1) / 2) / m));
  var lastBlock = Math.ceil((highPc / edgePc - (m - 1) / 2) / m);
  for (var s = firstBlock; s <= lastBlock; s++) {
    var radius = blockRadius(m, s) * edgePc;
    if (lowPc <= radius && radius <= highPc) out.push([s, radius]);
  }
  return out;
}

// [band, centre colatitude] of the rows (blocks) or bands (sectors) of one shell whose centre colatitude lies in
// [lowColat, highColat], north to south.
function windowBands(m, shell, lowColat, highColat) {
  var bands;
  var centres;
  if (m === 1) {
    if (shell === 0) return [[0, Math.PI / 2]];
    var cached = shellBands(shell);
    bands = cached.bands;
    centres = cached.centres;
  } else {
    var rows = drillBandCount(m, shell);
    bands = [];
    centres = [];
    for (var b = (rows - 1) / 2; b >= -(rows - 1) / 2; b--) {
      bands.push(b);
      centres.push(blockColatitude(m, shell, b));
    }
  }
  var start = bisectLeft(centres, lowColat);
  var stop = bisectRight(centres, highColat);
  var out = [];
  for (var i = start; i < stop; i++) out.push([bands[i], centres[i]]);
  return out;
}

// The wedges of a band of `count` whose centre longitude may lie within `halfWidth` of `lon0` (a superset).
function wedgeCandidates(count, lon0, halfWidth) {
  var all = [];
  if (halfWidth >= Math.PI || count === 1) {
    for (var w = 0; w < count; w++) all.push(w);
    return all;
  }
  var step = TWO_PI / count;
  var first = Math.floor((lon0 - halfWidth) / step - 0.5);
  var last = Math.ceil((lon0 + halfWidth) / step - 0.5);
  var seen = new Set();
  for (var k = first; k <= last; k++) seen.add(((k % count) + count) % count);
  return Array.from(seen).sort(function (p, q) { return p - q; });
}

// The cells of level m (a block size, or 1 for sectors) whose centre lies within `reachPc` of `centrePc`, ordered
// by shell, then band north to south, then wedge.
export function cellsInBall(m, centrePc, reachPc, edgePc) {
  var x0 = centrePc[0];
  var y0 = centrePc[1];
  var z0 = centrePc[2];
  var rho0 = Math.sqrt(x0 * x0 + y0 * y0 + z0 * z0);
  var theta0 = rho0 ? Math.acos(Math.max(-1.0, Math.min(1.0, z0 / rho0))) : 0.0;
  var lon0 = Math.atan2(y0, x0);
  var cells = [];
  levelShells(m, rho0 - reachPc - TOLERANCE_PC, rho0 + reachPc + TOLERANCE_PC, edgePc).forEach(function (entry) {
    var shell = entry[0];
    var radius = entry[1];
    var psi;
    var cosNeed;
    if (rho0 === 0 || radius === 0) {
      psi = Math.PI;
      cosNeed = -1.0;
    } else {
      cosNeed = (radius * radius + rho0 * rho0 - reachPc * reachPc) / (2 * radius * rho0);
      psi = cosNeed <= -1 ? Math.PI : (cosNeed >= 1 ? 0.0 : Math.acos(cosNeed));
    }
    // Block rows are wide: look one row's width beyond the cap for rows whose centre may still be inside.
    var slack = m > 1 ? Math.PI / drillBandCount(m, shell) : Math.PI / Math.max(1, bandRange(shell)[1] * 2 + 1);
    windowBands(m, shell, theta0 - psi - slack, theta0 + psi + slack).forEach(function (row) {
      var band = row[0];
      var colat = row[1];
      var count = m > 1 ? drillWedgeCount(m, shell, band) : wedgeCount(shell, band);
      var denominator = Math.sin(colat) * Math.sin(theta0);
      var wedges;
      if (Math.abs(denominator) < 1e-9 || radius === 0) {
        wedges = [];
        for (var w = 0; w < count; w++) wedges.push(w);
      } else {
        var cosLon = (cosNeed - Math.cos(colat) * Math.cos(theta0)) / denominator;
        var half = cosLon <= -1 ? Math.PI : (cosLon >= 1 ? 0.0 : Math.acos(cosLon));
        wedges = wedgeCandidates(count, lon0, half + TWO_PI / count);
      }
      wedges.forEach(function (wedge) {
        var cell = { m: m, shell: shell, band: band, wedge: wedge };
        var p = cellPositionPc(cell, edgePc);
        var dx = p[0] - x0;
        var dy = p[1] - y0;
        var dz = p[2] - z0;
        if (Math.sqrt(dx * dx + dy * dy + dz * dz) <= reachPc + TOLERANCE_PC) cells.push(cell);
      });
    });
  });
  return cells;
}

// The window's cells, [{m, shell, band, wedge}]: every cell of `centre`'s level whose centre lies within
// g * m / 2 edges of centre's. `centre` itself is always among them.
export function windowCells(centre, edgePc, g) {
  g = g || DEFAULT_WINDOW;
  var reach = g * centre.m * edgePc / 2.0;
  var cells = cellsInBall(centre.m, cellPositionPc(centre, edgePc), reach, edgePc);
  if (!cells.some(function (c) { return sameBlock(c, centre); })) {
    cells.push({ m: centre.m, shell: centre.shell, band: centre.band, wedge: centre.wedge });
    cells.sort(function (p, q) { return p.shell - q.shell || q.band - p.band || p.wedge - q.wedge; });
  }
  return cells;
}

// --- The galaxy's outline and how many sectors a block can hold ---------------------
//
// The galaxy's outline on the shell grid (planetgen.galaxy.skeleton.build_shell_extents), worked out from the
// analytic density bound: every shell out to the last whose equatorial band can hold a sector, each with the
// highest band it reaches. A shell's qualifying bands are one run around the plane, -high .. high (both caps when
// high is the shell's top band). Sectors allowed in a shell are those wedges of those bands.

// Whether sector band (shell, band) can hold a sector: the density bound at its centre line reaches the galaxy's
// sector threshold (shape.sector_min_density, the density a sector needs to expect one system), else the drawing
// floor.
export function bandAllowed(shell, band, edgePc, shape) {
  var centre = bandCentrePc(shell, band, edgePc);
  var threshold = shape.sector_min_density > 0 ? shape.sector_min_density : PRISM_MIN_DENSITY;
  return densityUpperBound(centre[0], Math.abs(centre[1]), shape) >= threshold;
}

// The highest |band| of shell `shell` with a cell that can qualify: the pole cap is checked first (if it qualifies
// the whole shell is in), and otherwise the qualifying bands are one run from the equator, found by bisection.
// Assumes the equatorial band qualifies.
function highestQualifyingBand(shell, edgePc, shape) {
  var top = bandRange(shell)[1];
  if (bandAllowed(shell, top, edgePc, shape)) return top;
  var low = 0;
  var high = top - 1;
  while (low < high) {
    var middle = Math.floor((low + high + 1) / 2);
    if (bandAllowed(shell, middle, edgePc, shape)) low = middle;
    else high = middle - 1;
  }
  return low;
}

// The outline: {high: Int32Array (per shell, the highest band reached, -1 for none), outerShell (the last shell
// with any, -1 for none), shapeless}. Without a shape the outline is every band of every shell out to the
// galaxy's radius.
export function galaxyOutline(edgePc, shape, galaxyRadius) {
  var limit = Math.max(0, Math.ceil(galaxyRadius / edgePc));
  if (!shape) {
    var all = new Int32Array(limit + 1);
    for (var k = 0; k <= limit; k++) all[k] = bandRange(k)[1];
    return { high: all, outerShell: limit, shapeless: true };
  }
  var extents = [];
  if (bandAllowed(0, 0, edgePc, shape)) {
    extents.push(0);
    var shell = 0;
    while (shell < limit + 1 && bandAllowed(shell + 1, 0, edgePc, shape)) {
      shell += 1;
      extents.push(highestQualifyingBand(shell, edgePc, shape));
    }
  }
  return { high: Int32Array.from(extents), outerShell: extents.length - 1, shapeless: false };
}

// The highest band sector shell `shell` reaches under `outline`, or -1.
export function outlineHigh(outline, shell) {
  return shell >= 0 && shell < outline.high.length ? outline.high[shell] : -1;
}

// Whether sector (shell, band, wedge) can hold anything under `outline`.
export function sectorInOutline(shell, band, wedge, outline) {
  if (!cellExists(shell, band, wedge)) return false;
  return Math.abs(band) <= outlineHigh(outline, shell);
}

// The sectors of the sector shell `k` that join block row `row` of `rows` equal colatitude steps at block level 3's
// own rows (`rows3` of them) -- the bands whose centre lies in level-3 rows [lo, hi) -- as [firstIndex, lastIndex+1) of the
// shell's north-to-south band list, held back to the bands the outline allows.
function bandSlice(k, lo, hi, rows3, high) {
  var cached = shellBands(k);
  var centres = cached.centres;
  var start = lo === 0 ? 0 : bisectLeft(centres, lo * Math.PI / rows3);
  var stop = hi === rows3 ? centres.length : bisectLeft(centres, hi * Math.PI / rows3);
  var top = bandRange(k)[1];
  // Bands top - high .. top + high are allowed (indices count down from the north cap).
  return [Math.max(start, top - high), Math.min(stop, top + high + 1)];
}

// The sectors of wedge range [first, last] of a band of `n`, whose centres lie in level-3 wedges [a, b) of a row of
// `w3`: integer arithmetic, so wedge lines on sector walls never round the wrong way (drill._wedge_cells).
function wedgeRun(n, a, b, w3) {
  var first = Math.max(0, Math.ceil((2 * a * n - w3) / (2 * w3)));
  var last = Math.min(n - 1, Math.ceil((2 * b * n - w3) / (2 * w3)) - 1);
  return Math.max(0, last - first + 1);
}

// How many sectors block `block` (level 3 or more) can hold under `outline`: the sectors of its sector shells whose centre
// joins it, held back to the bands the outline allows (drill.block_addresses). With `stopAtOne`, only whether there
// is any: 1 or 0, quickly.
function countSectors(block, outline, stopAtOne) {
  var m = block.m;
  var rows = drillBandCount(m, block.shell);
  var row = rowOf(m, block.shell, block.band);
  var blockWedges = drillWedgeCount(m, block.shell, block.band);
  var x = block.wedge;
  var total = 0;
  for (var k = block.shell * m; k < block.shell * m + m; k++) {
    var high = outlineHigh(outline, k);
    if (high < 0) continue;
    if (k === 0) {
      if (block.band === 0 && x === 0) {
        total += 1;
        if (stopAtOne) return 1;
      }
      continue;
    }
    var shell3 = Math.floor(k / DRILL_STEP);
    var rows3 = drillBandCount(DRILL_STEP, shell3);
    var f = rows3 / rows;
    var top = bandRange(k)[1];
    for (var row3 = row * f; row3 < (row + 1) * f; row3++) {
      var slice = bandSlice(k, row3, row3 + 1, rows3, high);
      if (slice[0] >= slice[1]) continue;
      var band3 = (rows3 - 1) / 2 - row3;
      var w3 = drillWedgeCount(DRILL_STEP, shell3, band3);
      var q = w3 / blockWedges;
      for (var i = slice[0]; i < slice[1]; i++) {
        var band = top - i;
        if (isCap(k, band)) {
          if (x === 0) {
            total += 1;
            if (stopAtOne) return 1;
          }
          continue;
        }
        var run = wedgeRun(wedgeCount(k, band), x * q, (x + 1) * q, w3);
        if (run > 0) {
          total += run;
          if (stopAtOne) return 1;
        }
      }
    }
  }
  return total;
}

// How many sectors drill block `block` can hold under `outline` (a sector: 1 or 0).
export function drillBlockTotal(block, outline) {
  if (block.m === 1) return sectorInOutline(block.shell, block.band, block.wedge, outline) ? 1 : 0;
  return countSectors(block, outline, false);
}

// Whether the block can hold any sector under `outline`.
export function drillBlockHolds(block, outline) {
  if (block.m === 1) return sectorInOutline(block.shell, block.band, block.wedge, outline);
  return countSectors(block, outline, true) > 0;
}

// A cell's mean density and its azimuthal mean ({density, mean}) over a few
// sample points, each weighted by the volume it stands for (r^2 sin(colatitude)):
// the outer samples of a deep cell cover more space than the inner ones. The page
// shades by both: density / mean is the arm factor, which picks out the spiral.
export function boundsDensity(b, shape) {
  var total = 0;
  var totalMean = 0;
  var weights = 0;
  for (var i = 0; i < SAMPLES_R; i++) {
    var r = b.r0 + ((i + 0.5) / SAMPLES_R) * (b.r1 - b.r0);
    for (var j = 0; j < SAMPLES_COLAT; j++) {
      var c = b.c0 + ((j + 0.5) / SAMPLES_COLAT) * (b.c1 - b.c0);
      var weight = r * r * Math.sin(c);
      var plane = r * Math.sin(c);
      var z = r * Math.cos(c);
      for (var k = 0; k < SAMPLES_LON; k++) {
        var t = b.t0 + ((k + 0.5) / SAMPLES_LON) * (b.t1 - b.t0);
        var parts = densityParts(plane * Math.cos(t), plane * Math.sin(t), z, shape);
        total += weight * parts[0];
        totalMean += weight * parts[1];
        weights += weight;
      }
    }
  }
  return { density: total / weights, mean: totalMean / weights };
}

// --- Geometry --------------------------------------------------------------
//
// One mesh for every cell in view: positions, per-vertex normals and a
// per-vertex cell index (the caller turns that into a color). A cell is a
// spherical box: an outer and an inner wall (spheres), a north and a south
// face (cones at the band's colatitudes) and two meridian sides (flat). A
// cap has no north face and no sides, the core is one sphere. Faces are wound
// counter-clockwise seen from outside, so front-face culling shows each cell's
// outside only. Each face carries its own 0..1 place (u along the longitude or
// the radius, v along the colatitude) for the shader to find its edges by.

// Points of a cell: radius r, colatitude c, longitude t.
function spherePoint(r, c, t) {
  var plane = r * Math.sin(c);
  return [plane * Math.cos(t), plane * Math.sin(t), r * Math.cos(c)];
}

export function buildPrismGeometry(cells) {
  var vertexCount = 0;
  var indexCount = 0;
  var plans = cells.map(function (b) {
    var lon = Math.max(1, Math.min(MAX_ARC_SEGMENTS, Math.ceil((b.t1 - b.t0) / ARC_SEGMENT_RAD)));
    var lat = Math.max(1, Math.min(MAX_ARC_SEGMENTS, Math.ceil((b.c1 - b.c0) / ARC_SEGMENT_RAD)));
    var plan = {
      lon: lon, lat: lat,
      inner: b.r0 > 0, north: b.c0 > 0, south: b.c1 < Math.PI, sides: !isFullTurn(b),
    };
    var walls = 1 + (plan.inner ? 1 : 0);
    vertexCount += walls * (lat + 1) * (lon + 1);
    indexCount += walls * 6 * lat * lon;
    var cones = (plan.north ? 1 : 0) + (plan.south ? 1 : 0);
    vertexCount += cones * 2 * (lon + 1);
    indexCount += cones * 6 * lon;
    if (plan.sides) {
      vertexCount += 2 * 2 * (lat + 1);
      indexCount += 2 * 6 * lat;
    }
    return plan;
  });
  var positions = new Float32Array(vertexCount * 3);
  var normals = new Float32Array(vertexCount * 3);
  var owners = new Uint32Array(vertexCount);
  var uvs = new Float32Array(vertexCount * 2);
  var indices = vertexCount > 65535 ? new Uint32Array(indexCount) : new Uint16Array(indexCount);
  var nv = 0;
  var ni = 0;

  function vertex(p, n, owner, u, v) {
    uvs[2 * nv] = u;
    uvs[2 * nv + 1] = v;
    var o = 3 * nv;
    positions[o] = p[0];
    positions[o + 1] = p[1];
    positions[o + 2] = p[2];
    normals[o] = n[0];
    normals[o + 1] = n[1];
    normals[o + 2] = n[2];
    owners[nv] = owner;
    nv++;
  }

  // Whether quad a-b-c-d, taken in that order, winds counter-clockwise around a's normal. (c - a) x (d - b) is the
  // quad's normal even when two of its corners coincide, as they do where a wedge meets the axis.
  function facesNormal(a, b, c, d) {
    var ux = positions[3 * c] - positions[3 * a];
    var uy = positions[3 * c + 1] - positions[3 * a + 1];
    var uz = positions[3 * c + 2] - positions[3 * a + 2];
    var vx = positions[3 * d] - positions[3 * b];
    var vy = positions[3 * d + 1] - positions[3 * b + 1];
    var vz = positions[3 * d + 2] - positions[3 * b + 2];
    var nx = uy * vz - uz * vy;
    var ny = uz * vx - ux * vz;
    var nz = ux * vy - uy * vx;
    return nx * normals[3 * a] + ny * normals[3 * a + 1] + nz * normals[3 * a + 2] >= 0;
  }

  // A grid of `rows + 1` rows of `columns + 1` vertices already written from `base`, joined into quads each wound to
  // face its own normal.
  function stitch(base, columns, rows) {
    for (var j = 0; j < rows; j++) {
      for (var i = 0; i < columns; i++) {
        var a = base + j * (columns + 1) + i;
        var b = a + 1;
        var d = a + columns + 1;
        var c = d + 1;
        if (facesNormal(a, b, c, d)) {
          indices[ni++] = a; indices[ni++] = b; indices[ni++] = c;
          indices[ni++] = a; indices[ni++] = c; indices[ni++] = d;
        } else {
          indices[ni++] = a; indices[ni++] = c; indices[ni++] = b;
          indices[ni++] = a; indices[ni++] = d; indices[ni++] = c;
        }
      }
    }
  }

  cells.forEach(function (b, n) {
    var plan = plans[n];
    var lon = plan.lon;
    var lat = plan.lat;
    var base;
    var j;
    var i;
    var t;
    var c;
    // The two walls: outer, then inner (none at the very centre).
    [b.r1, b.r0].forEach(function (r, wall) {
      if (wall === 1 && !plan.inner) return;
      var sign = wall === 0 ? 1 : -1;
      base = nv;
      for (j = 0; j <= lat; j++) {
        c = b.c0 + (j / lat) * (b.c1 - b.c0);
        for (i = 0; i <= lon; i++) {
          t = b.t0 + (i / lon) * (b.t1 - b.t0);
          var p = spherePoint(r, c, t);
          var radial = spherePoint(1, c, t);
          vertex(p, [sign * radial[0], sign * radial[1], sign * radial[2]], n, i / lon, j / lat);
        }
      }
      stitch(base, lon, lat);
    });
    // The two cones: north (colatitude c0), south (c1).
    [b.c0, b.c1].forEach(function (cone, face) {
      if (face === 0 ? !plan.north : !plan.south) return;
      var sign = face === 0 ? -1 : 1;
      base = nv;
      [b.r0, b.r1].forEach(function (r, row) {
        for (i = 0; i <= lon; i++) {
          t = b.t0 + (i / lon) * (b.t1 - b.t0);
          var p = spherePoint(r, cone, t);
          // The cone's normal points away from the cell: along or against increasing colatitude.
          vertex(p, [sign * Math.cos(cone) * Math.cos(t), sign * Math.cos(cone) * Math.sin(t), -sign * Math.sin(cone)], n, i / lon, row);
        }
      });
      stitch(base, lon, 1);
    });
    // The two meridian sides: start (t0), end (t1).
    if (plan.sides) {
      [b.t0, b.t1].forEach(function (side, face) {
        var sign = face === 0 ? -1 : 1;
        var normal = [-Math.sin(side) * sign, Math.cos(side) * sign, 0];
        base = nv;
        for (j = 0; j <= lat; j++) {
          c = b.c0 + (j / lat) * (b.c1 - b.c0);
          [b.r0, b.r1].forEach(function (r, row) {
            vertex(spherePoint(r, c, side), normal, n, row, j / lat);
          });
        }
        // Rows here run along the colatitude (two vertices each: inner, outer): a strip of `lat` quads.
        stitchStrip(base, lat);
      });
    }
  });

  // A strip of `segments + 1` pairs of vertices from `base`, [inner, outer] each, joined into quads.
  function stitchStrip(base, segments) {
    for (var s = 0; s < segments; s++) {
      var a = base + 2 * s;
      var b = a + 1;
      var d = a + 2;
      var c = d + 1;
      if (facesNormal(a, b, c, d)) {
        indices[ni++] = a; indices[ni++] = b; indices[ni++] = c;
        indices[ni++] = a; indices[ni++] = c; indices[ni++] = d;
      } else {
        indices[ni++] = a; indices[ni++] = c; indices[ni++] = b;
        indices[ni++] = a; indices[ni++] = d; indices[ni++] = c;
      }
    }
  }

  return { positions: positions, normals: normals, uvs: uvs, owners: owners, indices: indices };
}
