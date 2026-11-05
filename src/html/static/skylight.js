// html/static/skylight.js
//
// VIEW.32 (S14 of the System Sky plan, reports/system-sky-from-the-system-2026-11-04.md): how one star of the
// sky behind a system is drawn. The arithmetic only, free of the page and of three.js, so the sky points
// object (VIEW.33, skysphere.js), the band layer, the tests and the 9.1 view from a planet (VIEW.7 and
// VIEW.8, docs/design/orbits/sky-view.md) share one definition. It reuses the Galaxy Map's laws instead of
// making new ones: starmagnitude.js (absolute and apparent magnitude, the fade at the limiting magnitude)
// and starlight.js (the smallest a point of light is drawn).
//
//     m        = M + 5 log10(d / 10 pc)            apparent magnitude (starmagnitude.js)
//     peak     = min(1, 10^(-0.4 (m - 4.5)))       how bright the middle of the point is
//     opacity  = peak * smoothstep((m_lim - m) / RAMP)    peak, faded out over the last 1.5 magnitudes
//     sigma    = 0.7 px + 0.45 px per magnitude brighter than 4.5      the point's Gaussian width
//     size     = SIZE_SIGMAS * sigma, never under the Galaxy Map's core and never over MAX_SIZE_PX
//     colour   = the star's blackbody (cool stars take the model's offset), blended toward white as it
//                dims: sat = clip(0.35 + 0.65 (5.5 - m) / 4, 0.35, 1), rgb' = 1 - (1 - rgb) sat
//
// The limiting magnitude is a setting: SKY_LIMIT (7.5) by default and SKY_LIMIT_PHONE (6.5) on a phone, so
// the card draws about the same count the plan names (15,000 on a phone). The browser twin of
// planetgen/physics/photometry.py's `star_colour` (tests/js/skylight.test.mjs holds the two to one table)
// and, for the GLSL, of its own functions (SKY_LIGHT_GLSL is built from the constants below).
//
// Stand-ins (registered in GEN.725): the peak knee 4.5, the sigma constants and the desaturation are
// chosen by eye, as sky-view.md says; tier Unverified, tweak level L3. The limiting magnitudes follow the
// plan's counts and are tuned by eye too.

const VERSION_QUERY = new URL(import.meta.url).search;
const { STAR_MAGNITUDE_GLSL, absoluteMagnitude, apparentMagnitude, magnitudeShare } = await import(`./starmagnitude.js${VERSION_QUERY}`);
const { LIGHT_CORE_PX } = await import(`./starlight.js${VERSION_QUERY}`);

// The limiting magnitude of the sky: all a person sees under a dark sky, and the lower one for a phone.
export const SKY_LIMIT = 7.5;
export const SKY_LIMIT_PHONE = 6.5;

// A star this bright (or brighter) is drawn at full strength in its middle.
export const PEAK_MAGNITUDE = 4.5;
// The point's Gaussian width, px: SIGMA_PX at the peak magnitude and SIGMA_PER_MAG_PX more for every magnitude brighter.
export const SIGMA_PX = 0.7;
export const SIGMA_PER_MAG_PX = 0.45;
// A point is drawn this many sigmas wide, so the Gaussian has faded to nothing at its rim.
export const SIZE_SIGMAS = 4;
// No point is drawn wider than this (a star of magnitude -1.5 would otherwise take 14 px).
export const MAX_SIZE_PX = 12;
// The smallest a point is drawn: the Galaxy Map's smallest core.
export const MIN_SIZE_PX = LIGHT_CORE_PX[0];

// A dim star shows less of its colour: the share kept is SAT_MIN at SAT_FAINT_MAG and dimmer and 1 at SAT_FULL_MAG and brighter.
export const SAT_MIN = 0.35;
export const SAT_SLOPE = 0.65;
export const SAT_REF_MAG = 5.5;
export const SAT_SPAN_MAG = 4;

function clamp(x, lo, hi) {
  return Math.min(hi, Math.max(lo, x));
}

// The limiting magnitude to use: the phone's lower one when `phone`, else the setting's `limit` (default SKY_LIMIT).
export function skyLimit(phone, limit) {
  if (Number.isFinite(limit)) return limit;
  return phone ? SKY_LIMIT_PHONE : SKY_LIMIT;
}

// The strength of the middle of a star of apparent magnitude `m`: 1 at PEAK_MAGNITUDE and brighter, 0.2 at 6.5.
export function skyPeak(m) {
  return Math.min(1, Math.pow(10, -0.4 * (m - PEAK_MAGNITUDE)));
}

// How strongly a star of apparent magnitude `m` is drawn against the limit `limit`: its peak, faded out over the last RAMP magnitudes.
export function skyOpacity(m, limit) {
  return skyPeak(m) * magnitudeShare(m, limit);
}

// The width of the point's Gaussian, px, of a star of apparent magnitude `m`.
export function skySigmaPx(m) {
  return SIGMA_PX + SIGMA_PER_MAG_PX * Math.max(0, PEAK_MAGNITUDE - m);
}

// The full width of the point of a star of apparent magnitude `m`, px.
export function skySizePx(m) {
  return clamp(SIZE_SIGMAS * skySigmaPx(m), MIN_SIZE_PX, MAX_SIZE_PX);
}

// The share of its colour a star of apparent magnitude `m` keeps (1 whole, SAT_MIN at the faintest).
export function skySaturation(m) {
  return clamp(SAT_MIN + (SAT_SLOPE * (SAT_REF_MAG - m)) / SAT_SPAN_MAG, SAT_MIN, 1);
}

// Blends a linear colour `[r, g, b]` toward white for a star of apparent magnitude `m`.
export function skyTint(rgb, m) {
  const s = skySaturation(m);
  return [1 - (1 - rgb[0]) * s, 1 - (1 - rgb[1]) * s, 1 - (1 - rgb[2]) * s];
}

// ---- Colour: photometry.py's `star_colour` in the browser ----

const SECOND_RADIATION_CONSTANT_UM_K = 14387.76877; // c2 = h c / k, micrometre kelvin
const WAVELENGTH_MIN_NM = 380;
const WAVELENGTH_MAX_NM = 780;
export const BLACKBODY_ABOVE_K = 4500;
export const COOL_TAPER_FROM_K = 4000;
// [Teff K, dx, dy]: the BT-Settl chromaticity minus the blackbody's (photometry.py COOL_OFFSETS).
export const COOL_OFFSETS = [
  [2300, -0.01, -0.0222],
  [2500, -0.0495, -0.0154],
  [2800, -0.0273, -0.0104],
  [3000, -0.0208, -0.0079],
  [3300, -0.0018, -0.0029],
  [3500, 0.0026, -0.0039],
  [4000, 0.018, -0.0023],
];
const XYZ_TO_LINEAR_SRGB = [
  [3.2404542, -1.5371385, -0.4985314],
  [-0.969266, 1.8760108, 0.041556],
  [0.0556434, -0.2040259, 1.0572252],
];

function lobe(lam, mu, sigmaLow, sigmaHigh) {
  const t = (lam - mu) / (lam < mu ? sigmaLow : sigmaHigh);
  return Math.exp(-0.5 * t * t);
}

function cieBars(lam) {
  return [
    1.056 * lobe(lam, 599.8, 37.9, 31.0) + 0.362 * lobe(lam, 442.0, 16.0, 26.7) - 0.065 * lobe(lam, 501.1, 20.4, 26.2),
    0.821 * lobe(lam, 568.8, 46.9, 40.5) + 0.286 * lobe(lam, 530.9, 16.3, 31.1),
    1.217 * lobe(lam, 437.0, 11.8, 36.0) + 0.681 * lobe(lam, 459.0, 26.0, 13.8),
  ];
}

// The tristimulus values (Y = 1) of a blackbody of `temperatureK`, summed on the 1 nm grid with the trapezoid rule.
function blackbodyXyz(temperatureK) {
  let xs = 0;
  let ys = 0;
  let zs = 0;
  let prev = null;
  for (let nm = WAVELENGTH_MIN_NM; nm <= WAVELENGTH_MAX_NM; nm++) {
    const um = nm * 1e-3;
    const x = SECOND_RADIATION_CONSTANT_UM_K / (um * temperatureK);
    const power = x > 700 ? 0 : 1 / (Math.pow(um, 5) * Math.expm1(x));
    const bars = cieBars(nm);
    const cur = [power * bars[0], power * bars[1], power * bars[2]];
    if (prev) {
      xs += 0.5 * (prev[0] + cur[0]);
      ys += 0.5 * (prev[1] + cur[1]);
      zs += 0.5 * (prev[2] + cur[2]);
    }
    prev = cur;
  }
  return [xs / ys, 1, zs / ys];
}

function linearColour(xyz) {
  const rgb = XYZ_TO_LINEAR_SRGB.map((row) => Math.max(0, row[0] * xyz[0] + row[1] * xyz[1] + row[2] * xyz[2]));
  const peak = Math.max(rgb[0], rgb[1], rgb[2]);
  return peak > 0 ? rgb.map((c) => c / peak) : [1, 1, 1];
}

// The chromaticity offset (dx, dy) of a cool star from a blackbody, and its weight (1 at 4000 K and cooler, 0 at 4500 K and hotter).
function coolOffset(temperatureK) {
  if (temperatureK >= BLACKBODY_ABOVE_K) return [0, 0, 0];
  const weight = temperatureK <= COOL_TAPER_FROM_K ? 1 : (BLACKBODY_ABOVE_K - temperatureK) / (BLACKBODY_ABOVE_K - COOL_TAPER_FROM_K);
  const rows = COOL_OFFSETS;
  if (temperatureK <= rows[0][0]) return [rows[0][1], rows[0][2], weight];
  const last = rows[rows.length - 1];
  if (temperatureK >= last[0]) return [last[1] * weight, last[2] * weight, weight];
  for (let i = 0; i + 1 < rows.length; i++) {
    const a = rows[i];
    const b = rows[i + 1];
    if (temperatureK >= a[0] && temperatureK <= b[0]) {
      const f = (temperatureK - a[0]) / (b[0] - a[0]);
      return [a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2]), weight];
    }
  }
  return [0, 0, 0];
}

// The linear sRGB `[r, g, b]` (brightest channel 1: a hue, not a brightness) of a star of `temperatureK`, worked out exactly
// (a blackbody from 4500 K up, a cool star with the model's offset and no spread). The same as photometry.py's `star_colour`.
export function starColourExact(temperatureK) {
  const xyz = blackbodyXyz(temperatureK);
  const off = coolOffset(temperatureK);
  if (off[2] === 0) return linearColour(xyz);
  const total = xyz[0] + xyz[1] + xyz[2];
  const x = xyz[0] / total + off[0];
  const y = xyz[1] / total + off[1];
  return linearColour([x / y, 1, (1 - x - y) / y]);
}

// A table of the exact colour on COLOUR_GRID_POINTS temperatures spaced evenly in log between COLOUR_GRID_K[0] and [1], so a
// sky of tens of thousands of stars reads each colour with a lookup. Built on first use.
export const COLOUR_GRID_K = [1500, 50000];
export const COLOUR_GRID_POINTS = 128;
let colourTable = null;

function table() {
  if (!colourTable) {
    colourTable = [];
    const lo = Math.log(COLOUR_GRID_K[0]);
    const hi = Math.log(COLOUR_GRID_K[1]);
    for (let i = 0; i < COLOUR_GRID_POINTS; i++) {
      colourTable.push(starColourExact(Math.exp(lo + ((hi - lo) * i) / (COLOUR_GRID_POINTS - 1))));
    }
  }
  return colourTable;
}

// The linear colour of a star of `temperatureK` from the table (the ends hold beyond the grid); white with no temperature.
export function starColour(temperatureK) {
  if (!(temperatureK > 0)) return [1, 1, 1];
  const rows = table();
  const lo = Math.log(COLOUR_GRID_K[0]);
  const hi = Math.log(COLOUR_GRID_K[1]);
  const at = clamp(((Math.log(temperatureK) - lo) / (hi - lo)) * (COLOUR_GRID_POINTS - 1), 0, COLOUR_GRID_POINTS - 1);
  const i = Math.min(COLOUR_GRID_POINTS - 2, Math.floor(at));
  const f = at - i;
  const a = rows[i];
  const b = rows[i + 1];
  return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])];
}

// ---- One star ----

// How the sky draws a star of `luminositySol` (L_sun) and `temperatureK` at `distancePc` against the limit `limit`:
// `m` (apparent magnitude), `opacity`, `sizePx` and `rgb` (linear, blended toward white as it dims). A star with no
// luminosity is too faint to show (opacity 0).
export function skyStar(luminositySol, temperatureK, distancePc, limit = SKY_LIMIT) {
  const m = apparentMagnitude(absoluteMagnitude(luminositySol), distancePc);
  return { m: m, opacity: skyOpacity(m, limit), sizePx: skySizePx(m), rgb: skyTint(starColour(temperatureK), m) };
}

// The same arithmetic in GLSL, on top of STAR_MAGNITUDE_GLSL (`apparentMagnitude`, `magnitudeShare`, `fluxGain`):
// `float skyPeak(float m)`, `float skyOpacity(float m, float limit)`, `float skySizePx(float m)`,
// `float skySaturation(float m)` and `vec3 skyTint(vec3 rgb, float m)`.
export const SKY_LIGHT_GLSL = [
  STAR_MAGNITUDE_GLSL,
  "float skyPeak(float m) {",
  "  return min(1.0, pow(10.0, -0.4 * (m - " + PEAK_MAGNITUDE.toFixed(4) + ")));",
  "}",
  "float skyOpacity(float m, float limit) {",
  "  return skyPeak(m) * magnitudeShare(m, limit);",
  "}",
  "float skySizePx(float m) {",
  "  float sigma = " + SIGMA_PX.toFixed(4) + " + " + SIGMA_PER_MAG_PX.toFixed(4) + " * max(0.0, " + PEAK_MAGNITUDE.toFixed(4) + " - m);",
  "  return clamp(" + SIZE_SIGMAS.toFixed(4) + " * sigma, " + MIN_SIZE_PX.toFixed(4) + ", " + MAX_SIZE_PX.toFixed(4) + ");",
  "}",
  "float skySaturation(float m) {",
  "  return clamp(" + SAT_MIN.toFixed(4) + " + " + SAT_SLOPE.toFixed(4) + " * (" + SAT_REF_MAG.toFixed(4) + " - m) / " + SAT_SPAN_MAG.toFixed(4) + ", " + SAT_MIN.toFixed(4) + ", 1.0);",
  "}",
  "vec3 skyTint(vec3 rgb, float m) {",
  "  return vec3(1.0) - (vec3(1.0) - rgb) * skySaturation(m);",
  "}",
].join("\n");
