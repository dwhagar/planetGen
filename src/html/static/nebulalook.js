// html/static/nebulalook.js
//
// MAP.175: how a nebula is painted on the Galaxy Map. A nebula is a translucent sprite on a transparent canvas, so what shows is its colour at its core
// opacity composited over the map's background, which is black (`--map-bg`, UX.255 removed the light theme's looks). Every family reaches CONTRAST_TARGET
// (3:1, WCAG's floor for graphics) against black; tests/js/nebulalook.test.mjs checks it against the real colour in style.css.
//
// A look is [colour, core opacity, edge opacity], as galaxymap3d.js reads them. The dark family is a dust colour, a pale umber; MAP.181 settles the final
// one. The hues are those the Sector Map uses (planetgen/web/maps/starmap.py's _NEBULA_TYPE_COLORS, MAP.175 leaves that map as it is).

export const CONTRAST_TARGET = 3;

export const NEBULA_LOOKS = {
  diffuse: ["#eabdd6", 0.47, 0.14],
  emission: ["#ff6f91", 0.69, 0.25],
  reflection: ["#6fa8ff", 0.63, 0.22],
  planetary: ["#5be8c9", 0.66, 0.24],
  dark: ["#b08a67", 0.7, 0.3],
  default: ["#c9a8e0", 0.56, 0.22],
};

// "#rrggbb" as [r, g, b], 0..255.
export function parseHex(hex) {
  return [1, 3, 5].map(function (i) { return parseInt(hex.slice(i, i + 2), 16); });
}

// `colour` at `opacity` over `background` ([r, g, b] each), as [r, g, b].
export function compositeOver(colour, opacity, background) {
  return colour.map(function (c, i) { return opacity * c + (1 - opacity) * background[i]; });
}

function channel(c) {
  const s = c / 255;
  return s <= 0.04045 ? s / 12.92 : Math.pow((s + 0.055) / 1.055, 2.4);
}

// WCAG relative luminance of [r, g, b].
export function luminance(rgb) {
  return 0.2126 * channel(rgb[0]) + 0.7152 * channel(rgb[1]) + 0.0722 * channel(rgb[2]);
}

// WCAG contrast ratio of two [r, g, b] colours (1 to 21).
export function contrastRatio(a, b) {
  const la = luminance(a);
  const lb = luminance(b);
  return (Math.max(la, lb) + 0.05) / (Math.min(la, lb) + 0.05);
}

// A look's core, as it shows over `background`, against `background`.
export function lookContrast(look, background) {
  return contrastRatio(compositeOver(parseHex(look[0]), look[1], background), background);
}
