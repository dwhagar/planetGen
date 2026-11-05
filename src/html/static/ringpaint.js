// html/static/ringpaint.js
//
// The painter for a planet's rings (MAP.275). A pure function from the ring entry of a look record (planetgen/web/maps/ringlook.py: the edges in
// planet radii, the optical depth, the colour, the peak opacity and the gaps a moon digs at its resonances) to a strip of RGBA pixels running
// from the ring's inner edge (pixel 0) to its outer edge: ringlets of a broad ring from seeded noise, a few thin ringlets for a narrow ring, a soft
// haze for a dusty one, and a gap where each resonance falls. The strip is mapped onto the ring mesh by `ringGeometry` (radial UVs).
//
// What is stored is the ring (its edges, depth, ice share and the moons); what is drawn here is a picture of it: the ringlet pattern is seeded noise
// with no measurement behind it (a stand-in, gap-filling-plan 6.170, rows RNG7 and RNG8).
//
// Same ring, same seed, same size: same pixels. Do not change the noise or the constants without a changelog note.

const VERSION_QUERY = new URL(import.meta.url).search;
const { fbm, hashPixels, makeRandom } = await import(`./bodytexture.js${VERSION_QUERY}`);

export var STRIP_PIXELS = 512;
// How much of the strip's width each kind fades over at its edges (a share of the width).
var EDGE_SHARE = { dusty: 0.35, narrow: 0.04, broad: 0.03 };
// Broad rings: how far the ringlets swing the opacity either way (a share of the mean).
export var RINGLET_SWING = 0.35;
// Narrow rings: the number of ringlets, their width as a share of the ring's, and how much more opaque their cores are than the ring's mean.
export var NARROW_RINGLETS = 4;
export var NARROW_WIDTH = 0.025;
export var NARROW_CONCENTRATION = 6;
var NARROW_PEAK = 0.85;

function smoothstep(lo, hi, x) {
  var t = Math.max(0, Math.min(1, (x - lo) / (hi - lo)));
  return t * t * (3 - 2 * t);
}

// The opacity profile of the ring at share `x` (0 inner edge, 1 outer edge) before the gaps: `{ level, tone }`, `level` 0 to 1 of the peak opacity,
// `tone` 0.85 to 1.15 brightening or dimming the colour with the ringlets.
function profile(ring, seed, ringlets, x) {
  var edge = EDGE_SHARE[ring.kind] === undefined ? EDGE_SHARE.broad : EDGE_SHARE[ring.kind];
  var fade = smoothstep(0, edge, x) * smoothstep(0, edge, 1 - x);
  if (ring.kind === "dusty") {
    // A haze that thickens toward the middle: a smooth bell, no ringlets.
    return { level: fade * (0.5 + 0.5 * Math.sin(Math.PI * x)), tone: 1 };
  }
  if (ring.kind === "narrow") {
    var core = 0;
    for (var i = 0; i < ringlets.length; i++) {
      var d = (x - ringlets[i]) / NARROW_WIDTH;
      core = Math.max(core, Math.exp(-d * d));
    }
    return { level: fade * (0.12 + 0.88 * core), tone: 0.9 + 0.2 * core };
  }
  var wobble = (fbm(x, 0.5, seed + 71, 24, 4) - 0.5) * 2;
  var fine = (fbm(x, 0.2, seed + 97, 90, 3) - 0.5) * 2;
  var level = Math.max(0, Math.min(1, 1 + RINGLET_SWING * (0.7 * wobble + 0.3 * fine)));
  return { level: fade * level, tone: 1 + 0.15 * wobble };
}

// The strip, `n` RGBA pixels wide (4 * n bytes), for the ring entry `ring` and a seed.
export function paintRing(ring, seed, n) {
  var width = n > 0 ? n : STRIP_PIXELS;
  var out = new Uint8ClampedArray(width * 4);
  var rand = makeRandom((seed >>> 0) ^ 0x3c6ef372);
  var ringlets = [];
  for (var r = 0; r < NARROW_RINGLETS; r++) ringlets.push(0.1 + 0.8 * rand());
  var gaps = ring.gaps || [];
  var peak = ring.kind === "narrow" ? Math.min(NARROW_PEAK, ring.alpha * NARROW_CONCENTRATION) : ring.alpha;
  for (var p = 0; p < width; p++) {
    var x = (p + 0.5) / width;
    var here = profile(ring, seed >>> 0, ringlets, x);
    var open = 1;
    for (var g = 0; g < gaps.length; g++) {
      var center = gaps[g][0];
      var half = gaps[g][1] / 2;
      var depth = gaps[g][2];
      // Empty across the gap's inner third, fading to the ring's full opacity at its edge.
      open *= 1 - depth * (1 - smoothstep(half * 0.35, half, Math.abs(x - center)));
    }
    var o = p * 4;
    out[o] = Math.round(Math.min(255, ring.rgb[0] * here.tone));
    out[o + 1] = Math.round(Math.min(255, ring.rgb[1] * here.tone));
    out[o + 2] = Math.round(Math.min(255, ring.rgb[2] * here.tone));
    out[o + 3] = Math.round(255 * Math.max(0, Math.min(1, peak * here.level * open)));
  }
  return out;
}

// The picture's fingerprint, for the determinism tests.
export function ringFingerprint(ring, seed, n) {
  return hashPixels(paintRing(ring, seed, n));
}

// The ring's mesh, a flat annulus from `innerR` to `outerR` (scene units), its texture coordinate `u` the
// share of the way from the inner edge to the outer one, so the strip paints along the radius. `THREE` is the caller's own.
export function ringGeometry(THREE, innerR, outerR) {
  var geometry = new THREE.RingGeometry(innerR, outerR, 96, 1);
  var position = geometry.attributes.position;
  var uv = geometry.attributes.uv;
  for (var i = 0; i < position.count; i++) {
    var radius = Math.hypot(position.getX(i), position.getY(i));
    uv.setXY(i, (radius - innerR) / (outerR - innerR), 0.5);
  }
  uv.needsUpdate = true;
  return geometry;
}
