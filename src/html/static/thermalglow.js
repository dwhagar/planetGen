// html/static/thermalglow.js
//
// MAP.275 part 3: a world hot enough to glow. The look record's `thermal` entry (planetgen/web/maps/bodylook.py: `thermal_glow`) names the blackbody
// colour of its temperature and a strength from 0 (at the Draper point, where a body starts to shine red) to 1 (at THERMAL_FULL_K). The System Map
// adds it as light the body gives itself (a material's emissive colour over its own picture) and the object page's picture adds it to the shade, so
// the night side of a hot giant or a molten world is not black.
//
// Display only. The constants are stand-ins (gap-filling-plan 6.170, rows RNG9 and RNG10).

const VERSION_QUERY = new URL(import.meta.url).search;
const { parseHex } = await import(`./bodytexture.js${VERSION_QUERY}`);

// How strongly a full-strength glow adds to a pixel of the body's own colour (1 = as bright as the colour itself).
export var EMISSION_GAIN = 1.1;

// `{ rgb: [r, g, b] (0 to 255), strength }` for a look, or null for a body that does not glow.
export function thermalOf(look) {
  var hot = look && look.thermal;
  if (!hot || !(hot.strength > 0) || !hot.hex) return null;
  return { rgb: parseHex(hot.hex), strength: hot.strength };
}

// The channel `value` (0 to 255) of a pixel after the glow adds `colour` (0 to 255) scaled by the pixel's own brightness (`own`, 0 to 255), so a
// dark band glows less than a bright one.
export function glowed(value, colour, own, strength) {
  return Math.min(255, value + colour * (own / 255) * strength * EMISSION_GAIN);
}
