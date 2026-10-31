// html/static/solidpaint.js
//
// The painter for solid worlds (planets and moons) -- MAP.273, R4a of the body-rendering plan
// (docs/design/interface/body-rendering.md). A pure function from a look record (planetgen/web/maps/bodylook.py)
// and a size to equirectangular RGBA pixels: ocean, land and ice, in a deliberately rough few-tone vector style,
// with soft edges and a thin lighter coast.
//
// What the numbers mean on the picture:
//   * the sea level is chosen so the area under it (weighted by cos(latitude), the area an equirectangular row
//     really covers) equals the stored `ocean` fraction, so the page's numbers and the picture agree;
//   * ice sits where the mean temperature, lowered toward the poles, is below the freezing point of water, with a
//     ragged edge from the same seeded noise;
//   * land is one colour, two tones by height: the lab colour of the stored ground types (look.quality "composition",
//     MAP.274), or a stand-in from temperature (look.quality "estimate") where no ground is stored;
//   * a world with no air shows seeded craters;
//   * the air is seen through (MAP.350): `veil` (0 to 1, the air's optical depth) tints the ground toward the sky colour `glow` and
//     lowers its contrast, `haze` (haze and dust depth) lowers the contrast and softens it further, and `cloudCover` with `cloudTau`
//     lay a seeded, whitish cloud layer over it. A look without those fields is painted exactly as it was.
//
// Same look, same size: same pixels. Do not change the noise without a changelog note.

// Versioned like every sibling import (a static import would drop `?v=`).
const VERSION_QUERY = new URL(import.meta.url).search;
const { fbm, hashPixels, makeRandom, parseHex } = await import(`./bodytexture.js${VERSION_QUERY}`);

export var ICE_COLOR = "#f2f6fa";
export var SEA_ICE_COLOR = "#dbe8f2";
var FREEZING_K = 271;
var POLE_COOLING_K = 35;
var HISTOGRAM_BINS = 1024;

function clamp01(x) {
  return x < 0 ? 0 : x > 1 ? 1 : x;
}

function scaled(rgb, k) {
  return [Math.min(255, rgb[0] * k), Math.min(255, rgb[1] * k), Math.min(255, rgb[2] * k)];
}

function mix(a, b, t) {
  return [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
}

// The elevation at which `fraction` of the sphere's area lies below, from the elevation field `elev` (w * h values
// in [0, 1)), weighting each row by cos(latitude).
export function seaLevel(elev, w, h, fraction) {
  var bins = new Float64Array(HISTOGRAM_BINS);
  var total = 0;
  for (var y = 0; y < h; y++) {
    var weight = Math.cos((0.5 - (y + 0.5) / h) * Math.PI);
    for (var x = 0; x < w; x++) {
      var bin = Math.min(HISTOGRAM_BINS - 1, Math.floor(elev[y * w + x] * HISTOGRAM_BINS));
      bins[bin] += weight;
      total += weight;
    }
  }
  var target = clamp01(fraction) * total;
  var run = 0;
  for (var b = 0; b < HISTOGRAM_BINS; b++) {
    run += bins[b];
    if (run >= target) return (b + 1) / HISTOGRAM_BINS;
  }
  return 1;
}

// How the stored water state turns into what is drawn: the sea fraction, and whether the whole surface is ice or
// hidden under a haze deck.
function waterState(look) {
  var hydro = look.hydro;
  var ocean = typeof look.ocean === "number" ? clamp01(look.ocean) : 0;
  if (hydro === "ice") return { ocean: 0, allIce: true, haze: false };
  if (hydro === "ice-covered ocean") return { ocean: 0, allIce: true, haze: false, underIce: true };
  if (hydro === "vapour") return { ocean: 0, allIce: false, haze: true };
  if (hydro === "hycean") return { ocean: 1, allIce: false, haze: false };
  if (hydro === "dry") return { ocean: 0, allIce: false, haze: false };
  return { ocean: ocean, allIce: false, haze: false };
}

export var CLOUD_COLOR = "#f4f6f8";
var VEIL_TINT_MAX = 0.55;
var VEIL_FLATTEN_MAX = 0.4;
var HAZE_FLATTEN_MAX = 0.6;
var CLOUD_OPACITY_MAX = 0.85;
var CLOUD_TAU_SCALE = 2;

function luminance(c) {
  return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2];
}

// The air over the ground (MAP.350): flattens the picture toward its mean, tints it toward the sky colour at the same brightness, and lays
// the cloud over it. Changes `rgb` in place; returns how many extra blur passes the picture takes (a thick haze softens it).
function applyAir(look, rgb, w, h, seed) {
  var veil = clamp01(typeof look.veil === "number" ? look.veil : 0);
  var haze = clamp01(typeof look.haze === "number" ? look.haze : 0);
  var cover = clamp01(typeof look.cloudCover === "number" ? look.cloudCover : 0);
  var cloudTau = typeof look.cloudTau === "number" ? look.cloudTau : 0;
  var n = w * h;
  if (veil > 0 || haze > 0) {
    var mean = [0, 0, 0];
    for (var i = 0; i < n; i++) {
      mean[0] += rgb[i * 3];
      mean[1] += rgb[i * 3 + 1];
      mean[2] += rgb[i * 3 + 2];
    }
    mean = [mean[0] / n, mean[1] / n, mean[2] / n];
    var keep = (1 - VEIL_FLATTEN_MAX * veil) * (1 - HAZE_FLATTEN_MAX * haze);
    var tint = look.glow ? parseHex(look.glow) : null;
    var tintLum = tint ? luminance(tint) : 0;
    var amount = tint && tintLum > 0 ? VEIL_TINT_MAX * veil : 0;
    for (var j = 0; j < n; j++) {
      var c = [rgb[j * 3], rgb[j * 3 + 1], rgb[j * 3 + 2]];
      c = mix(mean, c, keep);
      if (amount > 0) {
        var k = luminance(c) / tintLum;
        c = mix(c, [tint[0] * k, tint[1] * k, tint[2] * k], amount);
      }
      rgb[j * 3] = Math.min(255, c[0]);
      rgb[j * 3 + 1] = Math.min(255, c[1]);
      rgb[j * 3 + 2] = Math.min(255, c[2]);
    }
  }
  if (cover > 0 && cloudTau > 0) {
    var field = new Float32Array(n);
    for (var y = 0; y < h; y++) {
      var v = (y + 0.5) / h;
      // More cloud near the equator and the middle latitudes than at about 30 degrees and the poles (a loose Hadley pattern).
      var lat = Math.abs(0.5 - v) * 2;
      var band = 0.08 * Math.cos(lat * Math.PI * 3);
      for (var x = 0; x < w; x++) field[y * w + x] = fbm((x + 0.5) / w, v, seed ^ 0xc10d, 5, 4) + band;
    }
    var cut = seaLevel(field.map(function (f) { return clamp01(f); }), w, h, 1 - cover);
    var opacity = CLOUD_OPACITY_MAX * (1 - Math.exp(-cloudTau / CLOUD_TAU_SCALE));
    var white = parseHex(CLOUD_COLOR);
    for (var m = 0; m < n; m++) {
      var t = clamp01((field[m] - cut) / 0.12 + 0.5);
      if (cover >= 1) t = 1;
      var a = opacity * t;
      if (a <= 0) continue;
      var shade = 0.9 + 0.2 * field[m];
      rgb[m * 3] += (Math.min(255, white[0] * shade) - rgb[m * 3]) * a;
      rgb[m * 3 + 1] += (Math.min(255, white[1] * shade) - rgb[m * 3 + 1]) * a;
      rgb[m * 3 + 2] += (Math.min(255, white[2] * shade) - rgb[m * 3 + 2]) * a;
    }
  }
  return haze > 0.5 || veil > 0.8 ? 1 : 0;
}

// Paints `look` into a w x h RGBA array (Uint8ClampedArray).
export function paintSolid(look, w, h) {
  var seed = look.seed >>> 0;
  var water = waterState(look);
  var land = parseHex(look.color);
  var sea = parseHex(look.oceanColor || "#1f4f8f");
  var ice = parseHex(ICE_COLOR);
  var seaIce = parseHex(SEA_ICE_COLOR);
  var temp = typeof look.tempK === "number" ? look.tempK : null;
  var hasAir = !!look.atm;

  var elev = new Float32Array(w * h);
  var iceNoise = new Float32Array(w * h);
  for (var y = 0; y < h; y++) {
    for (var x = 0; x < w; x++) {
      var u = (x + 0.5) / w;
      var v = (y + 0.5) / h;
      elev[y * w + x] = fbm(u, v, seed, 4, 4);
      iceNoise[y * w + x] = fbm(u, v, seed ^ 0x51ed, 6, 2);
    }
  }

  var level = water.ocean > 0 ? (water.ocean >= 1 ? 1.01 : seaLevel(elev, w, h, water.ocean)) : -1;
  var base = level < 0 ? 0 : Math.min(level, 1);
  var depthT = clamp01((look.oceanKm || 2) / 4);
  var deepMul = 0.85 - 0.3 * depthT;
  var shallowSpan = 0.08 + 0.05 * (1 - depthT);

  // Tone classes: 0 deep sea, 1 shallow sea, 2 low land, 3 high land.
  var tone = new Uint8Array(w * h);
  for (var i = 0; i < w * h; i++) {
    var e = elev[i];
    if (e < level) tone[i] = e < level - shallowSpan ? 0 : 1;
    else tone[i] = e < base + (1 - base) * 0.5 ? 2 : 3;
  }

  var rgb = new Float32Array(w * h * 3);
  for (var yy = 0; yy < h; yy++) {
    var lat = (0.5 - (yy + 0.5) / h) * Math.PI;
    var s = Math.abs(Math.sin(lat));
    var cap = s * s * s * s;
    for (var xx = 0; xx < w; xx++) {
      var k = yy * w + xx;
      var t = tone[k];
      var isSea = t < 2;
      var c;
      if (water.haze) {
        c = mix(land, [255, 255, 255], 0.35 + 0.25 * iceNoise[k]);
      } else if (water.allIce) {
        c = scaled(water.underIce ? mix(ice, sea, 0.18) : ice, 0.92 + 0.12 * elev[k]);
      } else {
        if (t === 0) c = scaled(sea, deepMul);
        else if (t === 1) c = scaled(sea, 1.08);
        else if (t === 2) c = scaled(land, 0.95);
        else c = scaled(land, 1.2);
        // Thin lighter coast: a land pixel with sea beside it.
        if (!isSea) {
          var x0 = xx > 0 ? xx - 1 : w - 1;
          var x1 = xx < w - 1 ? xx + 1 : 0;
          var nearSea =
            tone[yy * w + x0] < 2 ||
            tone[yy * w + x1] < 2 ||
            (yy > 0 && tone[(yy - 1) * w + xx] < 2) ||
            (yy < h - 1 && tone[(yy + 1) * w + xx] < 2);
          if (nearSea) c = mix(c, [255, 255, 255], 0.3);
        }
        if (temp !== null) {
          var localT = temp - POLE_COOLING_K * cap + (iceNoise[k] - 0.5) * 8;
          if (localT < FREEZING_K) c = mix(c, isSea ? seaIce : ice, isSea ? 0.85 : 0.95);
        }
      }
      rgb[k * 3] = c[0];
      rgb[k * 3 + 1] = c[1];
      rgb[k * 3 + 2] = c[2];
    }
  }

  var extraBlur = look.atm ? applyAir(look, rgb, w, h, seed) : 0;

  // Craters on a world with no air (and no ice sheet over it): seeded darker discs.
  if (!hasAir && !water.allIce && !water.haze) {
    var rnd = makeRandom(seed ^ 0x5bd1e995);
    var count = Math.round((w * h) / 160);
    for (var n = 0; n < count; n++) {
      var cx = Math.floor(rnd() * w);
      var cy = Math.floor(rnd() * h);
      var rad = 1 + Math.floor(rnd() * 3);
      var dark = 0.78 + 0.12 * rnd();
      for (var dy = -rad; dy <= rad; dy++) {
        for (var dx = -rad; dx <= rad; dx++) {
          if (dx * dx + dy * dy > rad * rad) continue;
          var py = cy + dy;
          if (py < 0 || py >= h) continue;
          var px = (cx + dx + w) % w;
          var idx = (py * w + px) * 3;
          rgb[idx] *= dark;
          rgb[idx + 1] *= dark;
          rgb[idx + 2] *= dark;
        }
      }
    }
  }

  // Soft edges: one 3 x 3 weighted blur, wrapping in longitude (a second pass under a thick haze).
  var out = new Uint8ClampedArray(w * h * 4);
  var weights = [1, 2, 1, 2, 4, 2, 1, 2, 1];
  for (var pass = 0; pass <= extraBlur; pass++) {
    for (var by = 0; by < h; by++) {
      for (var bx = 0; bx < w; bx++) {
        var r = 0;
        var g = 0;
        var b = 0;
        var sum = 0;
        var wi = 0;
        for (var oy = -1; oy <= 1; oy++) {
          var ny = by + oy;
          for (var ox = -1; ox <= 1; ox++, wi++) {
            if (ny < 0 || ny >= h) continue;
            var nx = (bx + ox + w) % w;
            var j = (ny * w + nx) * 3;
            var wt = weights[wi];
            r += rgb[j] * wt;
            g += rgb[j + 1] * wt;
            b += rgb[j + 2] * wt;
            sum += wt;
          }
        }
        var o = (by * w + bx) * 4;
        out[o] = Math.round(r / sum);
        out[o + 1] = Math.round(g / sum);
        out[o + 2] = Math.round(b / sum);
        out[o + 3] = 255;
      }
    }
    if (pass < extraBlur) {
      for (var q = 0; q < w * h; q++) {
        rgb[q * 3] = out[q * 4];
        rgb[q * 3 + 1] = out[q * 4 + 1];
        rgb[q * 3 + 2] = out[q * 4 + 2];
      }
    }
  }
  return out;
}

// The picture's fingerprint, for the determinism tests.
export function solidFingerprint(look, w, h) {
  return hashPixels(paintSolid(look, w, h));
}
