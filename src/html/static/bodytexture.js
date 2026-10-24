// html/static/bodytexture.js
//
// The shared parts of painting a body's surface from its look record (MAP.272, planetgen/web/maps/bodylook.py):
// seeded random numbers and wrapped value noise, the texture cache, the paint queue that does one paint per frame
// nearest first, the size tiers that decide how much detail a body gets, and a pixel hash for the determinism
// tests. The painters themselves (static/solidpaint.js for solid worlds) are pure functions from a look record and
// a size to RGBA pixels, so the same body gets the same picture on every machine and the tests need no browser.
//
// Everything here is display-only noise. It must not change between releases without a changelog note: a changed
// look is a visible change.

// 32-bit seeded generator (mulberry32). `seed` is the look record's `seed`.
export function makeRandom(seed) {
  var state = seed >>> 0;
  return function () {
    state = (state + 0x6d2b79f5) >>> 0;
    var t = state;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// A hash of two lattice coordinates and a seed, in [0, 1).
function lattice(ix, iy, seed) {
  var h = Math.imul(ix | 0, 0x27d4eb2d) ^ Math.imul(iy | 0, 0x165667b1) ^ Math.imul(seed | 0, 0x9e3779b1);
  h = Math.imul(h ^ (h >>> 15), 0x85ebca6b);
  h = Math.imul(h ^ (h >>> 13), 0xc2b2ae35);
  return ((h ^ (h >>> 16)) >>> 0) / 4294967296;
}

function smooth(t) {
  return t * t * (3 - 2 * t);
}

// Value noise in [0, 1) at (x, y), repeating every `wrapX` lattice cells in x so a texture's left and right edges
// meet (longitude wraps; latitude does not).
export function valueNoise(x, y, seed, wrapX) {
  var x0 = Math.floor(x);
  var y0 = Math.floor(y);
  var fx = smooth(x - x0);
  var fy = smooth(y - y0);
  function cell(ix, iy) {
    var wrapped = wrapX ? ((ix % wrapX) + wrapX) % wrapX : ix;
    return lattice(wrapped, iy, seed);
  }
  var a = cell(x0, y0);
  var b = cell(x0 + 1, y0);
  var c = cell(x0, y0 + 1);
  var d = cell(x0 + 1, y0 + 1);
  return a + (b - a) * fx + (c - a) * fy + (a - b - c + d) * fx * fy;
}

// Fractal noise in [0, 1): `octaves` layers of valueNoise, each twice the frequency and half the weight. `baseX` is
// the number of lattice cells around the whole texture at the first octave (so the wrap stays whole at every octave).
export function fbm(u, v, seed, baseX, octaves) {
  var total = 0;
  var weight = 0;
  var amp = 1;
  var cells = baseX;
  for (var o = 0; o < octaves; o++) {
    total += amp * valueNoise(u * cells, v * (cells / 2), seed + o * 101, cells);
    weight += amp;
    amp *= 0.5;
    cells *= 2;
  }
  return total / weight;
}

// A least-recently-used cache. `onEvict(value)` is called for each value pushed out (for texture.dispose()).
export function makeCache(capacity, onEvict) {
  var map = new Map();
  return {
    get: function (key) {
      if (!map.has(key)) return undefined;
      var value = map.get(key);
      map.delete(key);
      map.set(key, value);
      return value;
    },
    set: function (key, value) {
      if (map.has(key)) map.delete(key);
      map.set(key, value);
      while (map.size > capacity) {
        var oldest = map.keys().next().value;
        var evicted = map.get(oldest);
        map.delete(oldest);
        if (onEvict) onEvict(evicted, oldest);
      }
    },
    has: function (key) {
      return map.has(key);
    },
    size: function () {
      return map.size;
    },
  };
}

// The paint queue: jobs are added with a priority (smaller is more urgent, such as the body's distance from the
// camera) and `runOne()` does the most urgent, at most once per frame so a crowded system never stalls a frame.
// A key already waiting is replaced, not doubled, so a body that stays on screen is only queued once.
export function makeQueue() {
  var jobs = new Map();
  return {
    add: function (key, priority, run) {
      jobs.set(key, { priority: priority, run: run });
    },
    drop: function (key) {
      jobs.delete(key);
    },
    clear: function () {
      jobs.clear();
    },
    size: function () {
      return jobs.size;
    },
    runOne: function () {
      var bestKey = null;
      var best = null;
      jobs.forEach(function (job, key) {
        if (best === null || job.priority < best.priority) {
          best = job;
          bestKey = key;
        }
      });
      if (best === null) return false;
      jobs.delete(bestKey);
      best.run();
      return true;
    },
  };
}

// How much detail a body of this on-screen diameter gets. "point": a coloured dot, no texture. "disc": a flat tinted
// sphere. "low": the 128 x 64 texture. "high": the 256 x 128 texture. "shader": reserved for the close-up tier
// (R16, 9.1); nothing paints it yet, so it asks for the high texture.
export var TIER_DISC_PX = 3;
export var TIER_LOW_PX = 24;
export var TIER_HIGH_PX = 120;
export var TIER_SHADER_PX = 500;

export function tierForDiameter(px) {
  if (!(px >= TIER_DISC_PX)) return "point";
  if (px < TIER_LOW_PX) return "disc";
  if (px < TIER_HIGH_PX) return "low";
  return "high";
}

export var TIER_SIZES = { low: [128, 64], high: [256, 128] };

// FNV-1a over RGBA bytes: the picture's fingerprint for the determinism tests.
export function hashPixels(data) {
  var h = 0x811c9dc5;
  for (var i = 0; i < data.length; i++) {
    h ^= data[i];
    h = Math.imul(h, 0x01000193);
  }
  return h >>> 0;
}

export function parseHex(hex) {
  var m = /^#?([0-9a-f]{6})$/i.exec(hex || "");
  var n = m ? parseInt(m[1], 16) : 0x8a8f9c;
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}

// Reads the look record off a marker's `data-look` attribute value; null when absent or not valid JSON.
export function parseLook(text) {
  if (!text) return null;
  try {
    var look = JSON.parse(text);
    return look && look.v === 1 ? look : null;
  } catch (err) {
    return null;
  }
}
