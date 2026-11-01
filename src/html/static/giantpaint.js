// html/static/giantpaint.js
//
// The painter for gas and ice giants (R6a of the body-rendering plan, docs/design/interface/body-rendering.md section 5). A pure
// function from a look record (planetgen/web/maps/bodylook.py) and a size to equirectangular RGBA pixels: belts and zones in
// latitude, wavy edges, a few oval storms, and a palette set by what the clouds are made of.
//
// What the numbers mean on the picture (all display choices, a picture and not physics; the constants are stand-ins):
//   * the band count comes from the spin and the size: a fast spinner has narrower bands (the Rossby scale shrinks), so
//     `n = clamp(round(BAND_K * (R / R_JUPITER) / P_hours), 4, 24)`; with no stored spin the count is `BANDS_NO_SPIN`;
//   * the edge sharpness comes from the band count (many narrow bands are crisp, a few wide ones soft);
//   * the contrast between belts and zones comes from the cloud-top temperature, the way the clouds change: methane over
//     everything is nearly featureless (Neptune), ammonia clouds contrast most (Jupiter and Saturn), water clouds less, and
//     a cloudless or silicate-hot giant almost none;
//   * the palette comes from the stored gases: methane tints the air blue-green, ammonia clouds are cream, hydrogen sulfide
//     (the source of the ammonium hydrosulfide browns) warms the belts, water whitens the zones;
//   * 1 to 4 storms, more for a fast spinner, tinted by the same minor chromophore.
//
// Same look, same size: same pixels. Do not change the noise or the constants without a changelog note.

// Versioned like every sibling import (a static import would drop `?v=`).
const VERSION_QUERY = new URL(import.meta.url).search;
const { fbm, hashPixels, makeRandom, parseHex } = await import(`./bodytexture.js${VERSION_QUERY}`);

export var JUPITER_RADIUS_KM = 71492;
export var BAND_K = 120;
export var BANDS_MIN = 4;
export var BANDS_MAX = 24;
export var BANDS_NO_SPIN = 8;

// Belt and zone contrast (0 to 1, the lightness step) by cloud-top temperature, K: linear between the rows.
var CONTRAST_BY_TEMP = [
  [60, 0.06],     // methane air over everything: nearly featureless
  [130, 0.08],
  [160, 0.3],     // ammonia clouds
  [240, 0.3],
  [300, 0.18],    // water clouds
  [600, 0.12],
  [900, 0.06],    // cloudless or alkali-dark
  [1400, 0.04],
  [2400, 0.03],
];
var CREAM = [240, 232, 214];
var WARM_BELT = [150, 95, 60];
var COLD_BELT = [60, 120, 150];
var STORM_RED = [176, 84, 58];
var HOT_ONSET_K = 130;
var AMMONIA_END_K = 250;
var WATER_END_K = 800;
var H2S_FULL = 3e-4;
var CH4_FULL = 0.02;
var H2O_FULL = 1e-3;

function clamp01(x) {
  return x < 0 ? 0 : x > 1 ? 1 : x;
}

function mix(a, b, t) {
  return [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
}

function scaled(rgb, k) {
  return [Math.min(255, rgb[0] * k), Math.min(255, rgb[1] * k), Math.min(255, rgb[2] * k)];
}

function smoothstep(lo, hi, x) {
  var t = clamp01((x - lo) / (hi - lo));
  return t * t * (3 - 2 * t);
}

function interpolate(table, x) {
  if (x <= table[0][0]) return table[0][1];
  for (var i = 1; i < table.length; i++) {
    if (x <= table[i][0]) {
      var t = (x - table[i - 1][0]) / (table[i][0] - table[i - 1][0]);
      return table[i - 1][1] + (table[i][1] - table[i - 1][1]) * t;
    }
  }
  return table[table.length - 1][1];
}

// How many belt and zone bands a giant of this spin and size gets.
export function bandCount(spinH, radiusKm) {
  if (!(spinH > 0)) return BANDS_NO_SPIN;
  var size = radiusKm > 0 ? radiusKm / JUPITER_RADIUS_KM : 1;
  return Math.max(BANDS_MIN, Math.min(BANDS_MAX, Math.round((BAND_K * size) / spinH)));
}

// The recipe for a giant's picture: every number the painter uses, from the look record alone.
export function giantRecipe(look) {
  var rand = makeRandom((look.seed >>> 0) ^ 0x6a09e667);
  var temp = look.tempK > 0 ? look.tempK : 160;
  var gas = look.gas || {};
  var known = look.gas !== undefined;
  var base = parseHex(look.color);
  var bands = bandCount(look.spinH, look.radiusKm);
  var contrast = interpolate(CONTRAST_BY_TEMP, temp) * (0.8 + 0.4 * rand());

  // With no stored mixture the temperature regime stands in: ammonia in its band, methane below it.
  var nh3 = known ? clamp01((gas.NH3 || 0) / 1e-4) : temp >= HOT_ONSET_K && temp < AMMONIA_END_K ? 1 : 0;
  var ch4 = known ? clamp01((gas.CH4 || 0) / CH4_FULL) : temp < HOT_ONSET_K ? 1 : 0;
  var h2s = known ? clamp01((gas.H2S || 0) / H2S_FULL) : 0.5;
  var h2o = known ? clamp01((gas.H2O || 0) / H2O_FULL) : temp >= AMMONIA_END_K && temp < WATER_END_K ? 1 : 0;

  var ammoniaRegime = temp >= HOT_ONSET_K && temp < AMMONIA_END_K;
  var waterRegime = temp >= AMMONIA_END_K && temp < WATER_END_K;
  var chromo = ammoniaRegime ? clamp01(0.3 + 0.7 * h2s) * (0.7 + 0.3 * nh3) : waterRegime ? 0.15 * (0.5 + 0.5 * h2s) : 0;

  var zone = mix(base, CREAM, (ammoniaRegime ? 0.4 * (0.5 + 0.5 * nh3) : waterRegime ? 0.45 * (0.4 + 0.6 * h2o) : 0.12) * clamp01(contrast * 4));
  zone = scaled(zone, 1 + contrast * 0.35);
  var belt = scaled(base, 1 - contrast * 1.4);
  belt = mix(belt, ammoniaRegime || waterRegime ? WARM_BELT : COLD_BELT, ammoniaRegime || waterRegime ? chromo * 0.55 : ch4 * 0.3);
  var storm = chromo > 0.45 ? mix(zone, STORM_RED, 0.5 + 0.4 * chromo) : mix(zone, [255, 255, 255], 0.5);

  var fast = look.spinH > 0 ? clamp01((20 - look.spinH) / 14) : 0.4;
  var stormCount = Math.max(1, Math.min(4, 1 + Math.floor(rand() * 2) + Math.round(fast * 1.5)));
  if (contrast < 0.1) stormCount = 1;

  return {
    bands: bands,
    sharp: 0.25 + 0.7 * clamp01((bands - BANDS_MIN) / 16),
    contrast: contrast,
    chromo: chromo,
    waviness: 0.25 + 0.5 * rand(),
    storms: stormCount,
    stormStrength: clamp01(contrast * 4) * 0.6 + 0.2,
    zone: zone,
    belt: belt,
    stormColor: storm,
  };
}

// Equirectangular RGBA pixels, w by h (row 0 is the north pole, column 0 the 0-degree meridian).
export function paintGiant(look, w, h) {
  var r = giantRecipe(look);
  var seed = look.seed >>> 0;
  var rand = makeRandom(seed ^ 0xbb67ae85);
  var n = r.bands;

  // Band edges, evenly spread with a seeded jitter, and a seeded strength for each band.
  var edges = [0];
  var strength = [];
  for (var i = 1; i < n; i++) edges.push(i / n + ((rand() - 0.5) * 0.6) / n);
  edges.push(1);
  for (var k = 0; k < n; k++) strength.push(0.6 + 0.4 * rand());

  // Storms: a latitude inside a band, a longitude, and an oval size.
  var storms = [];
  for (var s = 0; s < r.storms; s++) {
    var band = Math.floor(rand() * n);
    var lat = Math.max(0.12, Math.min(0.88, (edges[band] + edges[band + 1]) / 2));
    var wide = 0.025 + 0.05 * rand();
    storms.push({ u: rand(), v: lat, rx: wide, ry: wide * 0.55 });
  }

  var out = new Uint8ClampedArray(w * h * 4);
  var edgeWidth = (1 - r.sharp) * 0.45 + 0.05; // transition width as a share of a band
  for (var y = 0; y < h; y++) {
    var v = (y + 0.5) / h;
    var polar = 1 - 0.6 * Math.pow(Math.abs(v - 0.5) * 2, 3);
    for (var x = 0; x < w; x++) {
      var u = (x + 0.5) / w;
      // Wavy edges: shift the latitude with seeded noise along the longitude.
      var wob = (fbm(u, v * 3, seed + 11, 6, 2) - 0.5) * r.waviness * 2;
      var t = Math.max(0, Math.min(0.9999, v + wob / n));
      var b = 0;
      while (b < n - 1 && t >= edges[b + 1]) b++;
      var width = edges[b + 1] - edges[b];
      // 0 for a belt, 1 for a zone, softened across the edges.
      var here = b % 2;
      var near = Math.min(t - edges[b], edges[b + 1] - t) / (width * edgeWidth);
      var tone = here + (0.5 - here) * (1 - smoothstep(0, 1, near));
      var level = clamp01(0.5 + (tone - 0.5) * strength[b] * 2 * polar);
      var rgb = mix(r.belt, r.zone, level);
      // Streaks: fine noise stretched along the longitude.
      var streak = (fbm(u, v * 6, seed + 29, 8, 3) - 0.5) * 2 * r.contrast * 0.9;
      var shade = 1 + streak;
      var pr = rgb[0] * shade;
      var pg = rgb[1] * shade;
      var pb = rgb[2] * shade;
      for (var q = 0; q < storms.length; q++) {
        var st = storms[q];
        var du = Math.abs(u - st.u);
        if (du > 0.5) du = 1 - du;
        var d = Math.sqrt((du / st.rx) * (du / st.rx) + ((v - st.v) / st.ry) * ((v - st.v) / st.ry));
        if (d < 1) {
          var a = (1 - smoothstep(0.4, 1, d)) * r.stormStrength;
          pr += (r.stormColor[0] - pr) * a;
          pg += (r.stormColor[1] - pg) * a;
          pb += (r.stormColor[2] - pb) * a;
        }
      }
      var o = (y * w + x) * 4;
      out[o] = Math.round(pr);
      out[o + 1] = Math.round(pg);
      out[o + 2] = Math.round(pb);
      out[o + 3] = 255;
    }
  }
  return out;
}

// The picture's fingerprint, for the determinism tests.
export function giantFingerprint(look, w, h) {
  return hashPixels(paintGiant(look, w, h));
}
