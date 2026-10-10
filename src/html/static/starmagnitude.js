// html/static/starmagnitude.js
//
// MAP.148: the star visibility law of the Galaxy Map. A star is drawn at an
// opacity that follows its APPARENT MAGNITUDE from the camera, against a
// limiting magnitude chosen each time the view changes so about N stars are
// on screen. docs/design/fly-through-view-distance.md (4.1) has the rule;
// this file is its arithmetic, free of the page, so the shader (which
// repeats `magnitudeShare` and `fluxGain` in GLSL, built from the constants
// below: STAR_MAGNITUDE_GLSL), the page's histogram and the tests share one
// definition.
//
//     M      = 4.83 - 2.5 log10(L / L_sun)          absolute magnitude
//     m      = M + 5 log10(d / 10 pc)               apparent magnitude
//     share  = smoothstep((m_lim - m) / RAMP)       the opacity
//     gain   = mix(FLUX_GAIN_MIN, 1, F / (F + 1))   brightness by flux,
//                F = 10^(-0.4 (m - (m_lim - FLUX_HALF_BELOW)))
//
// m_lim is the larger of two limits (limitingMagnitude takes it):
// - the HISTOGRAM limit: the magnitude at which the stars on screen,
//   counted by the zoom opacity MAP.153 already gives them, reach the
//   target (20,000; 8,000 on a phone), and
// - the CALIBRATION limit: where the dimmest stars the tiles carry still
//   show whole at the focus distance D. So the law only dims what stands
//   much farther than the thing looked at, and never thins the same stars
//   twice with MAP.153's rank birth radius.
// A black hole, neutron star or quasar has no luminosity to rank; its
// magnitude is PHENOMENON_MAGNITUDE, bright enough to stay lit.

export const SUN_ABSOLUTE_MAGNITUDE = 4.83;
// The width of the opacity ramp (1.5 magnitudes: the factor 4 the tile floors step by).
export const RAMP = 1.5;
export const TARGET_STARS = 20000;
export const TARGET_STARS_PHONE = 8000;
// Landmarks (black holes, neutron stars, quasars) are never dimmed away.
export const PHENOMENON_MAGNITUDE = -40;
// A star with no luminosity on record is too faint to show.
export const UNKNOWN_MAGNITUDE = 60;
// The brightness of a star at the faint end of the ramp, as a share of its own.
export const FLUX_GAIN_MIN = 0.5;
// Flux 1 sits this many magnitudes brighter than the limit.
export const FLUX_HALF_BELOW = 3;
export const HISTOGRAM_BINS = 64;
// The share of the shown stars' absolute magnitudes the calibration covers.
export const CALIBRATION_QUANTILE = 0.98;

function ease01(x) {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
}

export function absoluteMagnitude(luminositySol) {
  return luminositySol > 0 ? SUN_ABSOLUTE_MAGNITUDE - 2.5 * Math.log10(luminositySol) : UNKNOWN_MAGNITUDE;
}

// The apparent magnitude at `distancePc` of a star of absolute magnitude `M`.
export function apparentMagnitude(M, distancePc) {
  return M + 5 * Math.log10(Math.max(distancePc, 1e-6) / 10);
}

// How much of a star of apparent magnitude `m` shows against the limit.
export function magnitudeShare(m, limit) {
  return ease01((limit - m) / RAMP);
}

// The brightness of a star of apparent magnitude `m` against the limit, as
// a share of its own: tone-mapped flux, so a very bright star is capped at
// 1 and a star at the faint end shows at FLUX_GAIN_MIN.
export function fluxGain(m, limit) {
  const flux = Math.pow(10, -0.4 * (m - (limit - FLUX_HALF_BELOW)));
  return FLUX_GAIN_MIN + (1 - FLUX_GAIN_MIN) * (flux / (flux + 1));
}

// The histogram limit. `magnitudes` and `weights` (the zoom opacity of each
// star on screen; 1 when omitted) as arrays: the magnitude at which the
// weighted count of stars brighter than `limit - RAMP / 2` (the middle of
// the ramp) reaches `target`, found in HISTOGRAM_BINS bins from the
// brightest to the dimmest. With no more than `target` in all, one ramp past
// the dimmest, so all show whole.
export function histogramLimit(magnitudes, weights, target = TARGET_STARS) {
  let lo = Infinity;
  let hi = -Infinity;
  let total = 0;
  for (let i = 0; i < magnitudes.length; i++) {
    const w = weights ? weights[i] : 1;
    if (!(w > 0)) continue;
    total += w;
    if (magnitudes[i] < lo) lo = magnitudes[i];
    if (magnitudes[i] > hi) hi = magnitudes[i];
  }
  if (!(total > target)) return (Number.isFinite(hi) ? hi : 0) + RAMP;
  const span = Math.max(hi - lo, 1e-9);
  const bins = new Float64Array(HISTOGRAM_BINS);
  for (let i = 0; i < magnitudes.length; i++) {
    const w = weights ? weights[i] : 1;
    if (!(w > 0)) continue;
    const bin = Math.min(HISTOGRAM_BINS - 1, Math.floor(((magnitudes[i] - lo) / span) * HISTOGRAM_BINS));
    bins[bin] += w;
  }
  let sum = 0;
  for (let b = 0; b < HISTOGRAM_BINS; b++) {
    if (sum + bins[b] >= target) {
      const within = bins[b] > 0 ? (target - sum) / bins[b] : 0;
      return lo + ((b + within) / HISTOGRAM_BINS) * span + RAMP / 2;
    }
    sum += bins[b];
  }
  return hi + RAMP;
}

// The absolute magnitude of the dimmest stars shown (the quantile of
// `absolute` among those with a weight over 0.5).
export function dimmestShown(absolute, weights, quantile = CALIBRATION_QUANTILE) {
  const shown = [];
  for (let i = 0; i < absolute.length; i++) {
    if ((weights ? weights[i] : 1) > 0.5 && absolute[i] < UNKNOWN_MAGNITUDE) shown.push(absolute[i]);
  }
  if (!shown.length) return null;
  shown.sort((p, q) => p - q);
  return shown[Math.min(shown.length - 1, Math.floor(quantile * shown.length))];
}

// The calibration limit: the star of absolute magnitude `dimmest` at the
// focus distance `focusPc` shows whole.
export function calibrationLimit(dimmest, focusPc) {
  return dimmest === null ? -Infinity : apparentMagnitude(dimmest, focusPc) + RAMP;
}

// The limit to use: whichever is larger, so nothing at or nearer the
// focus is thinned by it.
export function limitingMagnitude(histogram, calibration) {
  return Math.max(histogram, calibration);
}

// Moves `current` toward `goal` over `seconds` with time constant `tau`
// (no jump: the limit changes smoothly as the view does).
export function approachLimit(current, goal, seconds, tau = 0.3) {
  if (!Number.isFinite(current) || !Number.isFinite(goal) || !(tau > 0)) return goal;
  return goal + (current - goal) * Math.exp(-Math.max(0, seconds) / tau);
}

// The same arithmetic in GLSL: `float magnitudeShare(float m, float limit)`
// and `float fluxGain(float m, float limit)`, plus
// `float apparentMagnitude(float M, float distancePc)`.
export const STAR_MAGNITUDE_GLSL = [
  "float apparentMagnitude(float M, float distancePc) {",
  "  return M + 5.0 * 0.30103 * log2(max(distancePc, 1e-6) / 10.0);",
  "}",
  "float magnitudeShare(float m, float limit) {",
  "  float t = clamp((limit - m) / " + RAMP.toFixed(4) + ", 0.0, 1.0);",
  "  return t * t * (3.0 - 2.0 * t);",
  "}",
  "float fluxGain(float m, float limit) {",
  "  float flux = pow(10.0, -0.4 * (m - (limit - " + FLUX_HALF_BELOW.toFixed(4) + ")));",
  "  return " + FLUX_GAIN_MIN.toFixed(4) + " + " + (1 - FLUX_GAIN_MIN).toFixed(4) + " * (flux / (flux + 1.0));",
  "}",
].join("\n");
