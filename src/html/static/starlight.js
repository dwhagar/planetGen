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
// The browser twin of lib/starmap.py's `star_light_boost` and
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
