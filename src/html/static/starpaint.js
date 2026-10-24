// html/static/starpaint.js
//
// The painter for stars (MAP.276, R7 of the body-rendering plan, docs/design/interface/body-rendering.md): a pure
// function from a look record (planetgen/web/maps/bodylook.py) and a size to equirectangular RGBA pixels, plus the
// glow numbers a luminosity gives. The disc colour is the look's `color` (the blackbody colour of the star's
// temperature, the one the sky colour assumes too); this file only adds the texture the class calls for:
//
//   * O, B, A: nearly featureless, tiny granulation contrast;
//   * F, G: fine granulation, a few dark spots;
//   * K, M: larger cooler cells, more and bigger spots (M most), and on a flare star an occasional bright patch;
//   * giants and supergiants: a few huge convection cells, no spots;
//   * white dwarfs: a smooth hot disc, no granulation.
//
// A star is self-luminous, so the texture never goes below 0.7 of the base colour. Same look, same size: same
// pixels. Do not change the noise without a changelog note.

// Versioned like every sibling import (a static import would drop `?v=`).
const VERSION_QUERY = new URL(import.meta.url).search;
const { fbm, hashPixels, makeRandom, parseHex } = await import(`./bodytexture.js${VERSION_QUERY}`);

var GIANT_CLASSES = { "0": 1, IA: 1, "IA+": 1, IAB: 1, IB: 1, II: 1, III: 1 };

// "dwarf", "giant" or "white dwarf" from the Yerkes class; a record without one is treated as a dwarf.
export function starKind(look) {
  var yerkes = (look.yerkes || "").toUpperCase();
  if (yerkes === "VII" || yerkes === "D") return "white dwarf";
  if (GIANT_CLASSES[yerkes]) return "giant";
  return "dwarf";
}

// The spectral letter, "" when the type is missing.
export function spectralLetter(look) {
  return ((look.spec || "")[0] || "").toUpperCase();
}

// How the surface looks, from the class: granulation cell count around the star, octaves, contrast, spot count and
// spot size (as a fraction of the texture width), and the chance of a bright flare patch.
export function surfaceRecipe(look) {
  var kind = starKind(look);
  var letter = spectralLetter(look);
  if (kind === "white dwarf") return { cells: 0, octaves: 0, contrast: 0, spots: 0, spotSize: 0, flare: 0 };
  if (kind === "giant") return { cells: 3, octaves: 2, contrast: 0.22, spots: 0, spotSize: 0, flare: 0 };
  switch (letter) {
    case "O":
    case "B":
    case "A":
      return { cells: 16, octaves: 3, contrast: 0.04, spots: 0, spotSize: 0, flare: 0 };
    case "F":
      return { cells: 24, octaves: 3, contrast: 0.1, spots: 1, spotSize: 0.02, flare: 0 };
    case "G":
      return { cells: 24, octaves: 3, contrast: 0.12, spots: 3, spotSize: 0.02, flare: 0 };
    case "K":
      return { cells: 14, octaves: 3, contrast: 0.16, spots: 5, spotSize: 0.03, flare: 0 };
    case "M":
      return { cells: 10, octaves: 3, contrast: 0.2, spots: 8, spotSize: 0.045, flare: 0.3 };
    default:
      return { cells: 16, octaves: 3, contrast: 0.1, spots: 0, spotSize: 0, flare: 0 };
  }
}

// The corona shell around a star as a multiple of its radius, and its strength, from the luminosity in solar units:
// a brighter star gets a wider, stronger glow, capped; a faint one a small one. 1 solar luminosity gives 1.4 and 0.75.
export function glowFor(lumSolar) {
  var lum = typeof lumSolar === "number" && lumSolar > 0 ? lumSolar : 1;
  var log = Math.log10(Math.max(lum, 1e-4));
  return {
    scale: Math.max(1.2, Math.min(1.9, 1.4 + 0.07 * log)),
    strength: Math.max(0.5, Math.min(1, 0.75 + 0.08 * log)),
  };
}

// Paints `look` into a w x h RGBA array (Uint8ClampedArray).
export function paintStar(look, w, h) {
  var seed = look.seed >>> 0;
  var base = parseHex(look.color);
  var recipe = surfaceRecipe(look);
  var shade = new Float32Array(w * h);
  for (var y = 0; y < h; y++) {
    var v = (y + 0.5) / h;
    for (var x = 0; x < w; x++) {
      var n = recipe.cells ? fbm((x + 0.5) / w, v, seed, recipe.cells, recipe.octaves) : 0.5;
      shade[y * w + x] = 1 + (n - 0.5) * 2 * recipe.contrast;
    }
  }

  var rnd = makeRandom(seed ^ 0x2c1b3c6d);
  function blot(multiplier, count, size) {
    for (var s = 0; s < count; s++) {
      var cx = rnd() * w;
      var lat = (rnd() * 2 - 1) * 0.7; // spots sit within about 40 degrees of the equator
      var cy = (0.5 - lat / Math.PI) * h;
      var rad = Math.max(1.5, size * w * (0.6 + 0.8 * rnd()));
      for (var dy = -Math.ceil(rad); dy <= Math.ceil(rad); dy++) {
        var py = Math.round(cy) + dy;
        if (py < 0 || py >= h) continue;
        for (var dx = -Math.ceil(rad); dx <= Math.ceil(rad); dx++) {
          var d = Math.sqrt(dx * dx + dy * dy) / rad;
          if (d >= 1) continue;
          var px = (((Math.round(cx) + dx) % w) + w) % w;
          var soft = 1 - d * d;
          shade[py * w + px] *= 1 + (multiplier - 1) * soft;
        }
      }
    }
  }
  if (recipe.spots) blot(0.55, recipe.spots, recipe.spotSize);
  if (recipe.flare && rnd() < recipe.flare) blot(1.5, 1, recipe.spotSize * 0.7);

  var out = new Uint8ClampedArray(w * h * 4);
  for (var i = 0; i < w * h; i++) {
    var k = Math.max(0.7, Math.min(1.25, shade[i]));
    out[i * 4] = Math.min(255, Math.round(base[0] * k));
    out[i * 4 + 1] = Math.min(255, Math.round(base[1] * k));
    out[i * 4 + 2] = Math.min(255, Math.round(base[2] * k));
    out[i * 4 + 3] = 255;
  }
  return out;
}

export function starFingerprint(look, w, h) {
  return hashPixels(paintStar(look, w, h));
}
