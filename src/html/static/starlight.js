// html/static/starlight.js
//
// MAP.87: how much brighter a star's point of light is drawn than its
// luminosity alone would make it, on the Galaxy Map and the Sector Map
// alike. Boss: the dimmest red dwarfs "4 times as bright as they are now
// and when we approach the 1000+ sol lum mark it evens out to be not any
// brighter. That's just in how it's displayed." The boost is
// BOOST_AT_DIM at BOOST_LOG_LUMINOSITY[0] (1e-4 L_sun) and 1 from
// BOOST_LOG_LUMINOSITY[1] (1000 L_sun) up; in between its exponent falls
// as share^BOOST_TAPER_POWER of the log range, slowly enough that a
// brighter star is never drawn fainter than a dimmer one (starmap.py's
// `_BOOST_TAPER_POWER` says why). A Sun gets about 2.2.
//
// A boost of b means b times the halo's light: its strength grows by
// sqrt(b) and its area by sqrt(b) (its width by b^(1/4)); the core, sized
// and lit by the star itself, is left alone (boostLight).
//
// The browser twin of planetgen/web/maps/starmap.py's `star_light_boost` and
// `_boost_light`, which work out the Sector Map's points before the page
// is sent; tests/js/starlight.test.mjs checks the two agree.

export const BOOST_AT_DIM = 4;
export const BOOST_LOG_LUMINOSITY = [-4, 3];
export const BOOST_TAPER_POWER = 1.5;

// The boost for a star of `luminositySol` solar luminosities (1 when it
// has none on record).
export function starLightBoost(luminositySol) {
  if (!(luminositySol > 0)) {
    return 1;
  }
  const lo = BOOST_LOG_LUMINOSITY[0];
  const hi = BOOST_LOG_LUMINOSITY[1];
  const share = Math.min(1, Math.max(0, (Math.log10(luminositySol) - lo) / (hi - lo)));
  return Math.pow(BOOST_AT_DIM, 1 - Math.pow(share, BOOST_TAPER_POWER));
}

// A point of light's halo (`sizePx`, its full width, and `glow`, its
// strength) drawn `boost` times as bright.
export function boostLight(light, boost) {
  const root = Math.sqrt(boost);
  return { sizePx: light.sizePx * Math.sqrt(root), glow: light.glow * root };
}

// MAP.241: a star's own look at sector depth, one for the Galaxy Map and the Sector Map. Boss: "stars on the starmap should
// glow brighter but not be bigger", so a star's size stays in a narrow band and luminosity is carried by how strongly it
// glows (and by the core's opacity), its radius by a small core, its temperature by its blackbody colour. The sizes are
// screen pixels. planetgen/web/maps/starmap.py's `_star_light` has the same numbers (tests/test_starmap.py and
// tests/js/starlight.test.mjs hold both to one table).
export const LIGHT_LOG_LUMINOSITY = [-4, 6];
export const LIGHT_LOG_RADIUS = [-1, 3];
export const LIGHT_CORE_PX = [1.6, 2.4];
export const LIGHT_SIZE_PX = [8.5, 11];
export const LIGHT_GLOW = [0.1, 0.9];
// The glow grows as the luminosity share to this power: slowly for the faint stars (MAP.87's boost already draws them
// brighter) and steeply for the bright ones, which keeps the halo's light (strength times area) rising with luminosity
// everywhere, as the boost taper falls away between 1 and 1000 suns.
export const LIGHT_GLOW_POWER = 2.5;
export const LIGHT_BRIGHT = [0.75, 1];

function logShare(value, range) {
  const share = (Math.log10(Math.max(value, 1e-12)) - range[0]) / (range[1] - range[0]);
  return Math.min(1, Math.max(0, share));
}

function lerp(range, share) {
  return range[0] + (range[1] - range[0]) * share;
}

// The look of a star of `luminositySol` (L_sun) and `radiusSol` (R_sun, or none: guessed from the luminosity): `share`
// (where its luminosity falls on the ramp, 0..1), `corePx` (the core's diameter), `sizePx` (the halo's full diameter, never
// less than the core plus a pixel's rim each side), `glow` (the halo's strength) and `bright` (the core's opacity).
export function starLight(luminositySol, radiusSol) {
  const luminosity = luminositySol > 0 ? luminositySol : 1;
  const radius = radiusSol > 0 ? radiusSol : Math.pow(luminosity, 0.35);
  const share = logShare(luminosity, LIGHT_LOG_LUMINOSITY);
  const corePx = lerp(LIGHT_CORE_PX, logShare(radius, LIGHT_LOG_RADIUS));
  const halo = boostLight({ sizePx: lerp(LIGHT_SIZE_PX, share), glow: lerp(LIGHT_GLOW, Math.pow(share, LIGHT_GLOW_POWER)) }, starLightBoost(luminosity));
  return { share: share, corePx: corePx, sizePx: Math.max(halo.sizePx, 2 * corePx + 2), glow: halo.glow, bright: lerp(LIGHT_BRIGHT, share) };
}
