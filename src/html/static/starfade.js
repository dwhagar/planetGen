// html/static/starfade.js
//
// MAP.153: stars fade in with the zoom instead of all at once at a tile
// level change. Boss (2026-10-09 22:30Z): "I want a very smooth transition
// where stars and objects are slowly added as one zooms in."
// docs/design/zoom-star-visibility.md (3.1 and 5, stage 1) has the rule and
// the evidence; this file is its arithmetic, free of the page, so the
// Galaxy Map's shader and its tests share one definition (the shader
// repeats `zoomShare` and `levelGlide` in GLSL; keep them in step).
//
// A star listed in a tile has a rank r in that tile's list (1 = most
// luminous, the order the server sends). Its BIRTH RADIUS is the camera
// radius at which it starts to appear,
//
//     R_b = R* 2^W (N0 / r)^(1/3),   R* = edge / (2 FETCH_RADIUS_FACTOR),
//
// R* being the closest camera radius that tile level serves, and it is
// fully shown W octaves of zoom later (R_b / 2^W). The exponent 1/3 keeps
// the count of drawn stars per view level (a tile one level finer holds an
// eighth of the volume). Over the octave in which a level is used its
// opacity cross-fades from the one its rank in the PARENT tile's list gives
// it (nothing when the parent did not list it) to the one its own rank
// gives it; the finest "detail" tiles, which list every star around the
// target, use their own rank only. The opacity depends on the camera
// radius alone, so a view always draws the same picture and zooming out
// removes stars as smoothly as zooming in adds them.

export const FADE_OCTAVES = 1;
export const REFERENCE_RANK_COUNT = 400;
export const FETCH_RADIUS_FACTOR = 1.6;
// A point object (black hole, neutron star, quasar) is not ranked in
// stage 1: a birth radius past any camera radius shows it fully.
export const ALWAYS_SHOWN = 1e12;

export function smoothstep(x) {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
}

// The camera radius (pc) at which a star ranked `rank` (1 = first) in the
// list of a tile `edgePc` across starts to appear.
export function birthRadius(edgePc, rank, factor = FETCH_RADIUS_FACTOR) {
  const closest = edgePc / (2 * factor);
  return closest * Math.pow(2, FADE_OCTAVES) * Math.pow(REFERENCE_RANK_COUNT / rank, 1 / 3);
}

// How much of a star with birth radius `birth` shows at camera radius
// `radius`: 0 at or beyond `birth`, 1 FADE_OCTAVES of zoom closer. A
// birth of 0 (not listed) shows nothing.
export function zoomShare(birth, radius) {
  return birth > 0 ? smoothstep(Math.log2(birth / radius) / FADE_OCTAVES) : 0;
}

// 0 at the zoomed-out end of the octave a tile level serves
// (radius = edge / factor) to 1 at the zoomed-in end (edge / 2 factor).
export function levelGlide(edgePc, radius, factor = FETCH_RADIUS_FACTOR) {
  const far = edgePc / factor;
  return smoothstep((far - radius) / (far - far / 2));
}

// A star's opacity from zoom alone. `fade` is [own birth radius, parent
// birth radius (0 = not listed there), tile edge, detail (0 or 1)].
export function starZoomOpacity(fade, radius, factor = FETCH_RADIUS_FACTOR) {
  const own = zoomShare(fade[0], radius);
  if (fade[3] > 0.5) {
    return own;
  }
  const glide = levelGlide(fade[2], radius, factor);
  return (1 - glide) * zoomShare(fade[1], radius) + glide * own;
}

// A point object (stage 1: never ranked).
export const POINT_FADE = [ALWAYS_SHOWN, ALWAYS_SHOWN, 0, 1];

// The 1-based ranks of a tile's stars, `{ b<id>: rank, s<id>: rank }`,
// from its bright list (`stars`) and its generated list (`generated`),
// each in the order the server sent them (most luminous first).
export function tileRanks(tile) {
  const ranks = new Map();
  (tile.stars || []).forEach((star, i) => ranks.set(`b${star.id}`, i + 1));
  (tile.generated || []).forEach((star, i) => ranks.set(`s${star.id}`, i + 1));
  return ranks;
}

// The sum of `starZoomOpacity` over `fades` at `radius`.
export function opacitySum(fades, radius, factor = FETCH_RADIUS_FACTOR) {
  let total = 0;
  for (const fade of fades) {
    total += starZoomOpacity(fade, radius, factor);
  }
  return total;
}
