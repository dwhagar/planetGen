// html/static/poiicons.js
//
// The point-of-interest symbols (MAP.320): seven vector icons, one SVG file each in static/poi/, drawn the same at map
// size in a page (`<svg class="poi-icon"><use href="poi/<id>.svg#poi"/></svg>`) and in the Galaxy Map's marker layer
// (galaxymarkerlayer.js draws the file into a texture). This file is the list of them and what each is called; it has
// no DOM at import and no three.js, so it runs under plain node for tests (tests/js/poiicons.test.mjs). The Python side
// is `planetgen/web/lib/poi_icons.py`, and a test holds the two lists and the files' own titles to each other.
//
// A symbol is a silhouette (colour is never the only sign) filled in the colour of its token:
//   id       the file's name without `.svg`, also the `data-poi` of an element drawing it;
//   label    its accessible name, also the file's `<title>`;
//   token    the CSS custom property (style.css) that colours it, `null` for the ones the habitability score colours
//            (`TIER_TOKENS`): the habitable world and the inhabited world.

export const POI_ICONS = [
  { id: "habitable-world", label: "Habitable world", token: null },
  { id: "uninhabitable-world", label: "Uninhabitable world", token: "--poi-brown" },
  { id: "inhabited-world", label: "Inhabited world", token: null },
  { id: "gas-giant", label: "Gas giant", token: "--poi-gas" },
  { id: "asteroid-field", label: "Asteroid field", token: "--poi-grey" },
  { id: "sector", label: "Sector", token: "--accent" },
  { id: "region", label: "Region", token: "--accent" },
];

// The colour of a habitability score (the PHI-6 tier, 0 Blue to 3 Red): the tokens the object pages already use.
export const TIER_TOKENS = ["--tier-blue", "--tier-green", "--tier-yellow", "--tier-red"];

export const TIER_NAMES = ["Blue", "Green", "Yellow", "Red"];

// The symbol with id `id`, or `null`.
export function poiIcon(id) {
  return POI_ICONS.find(function (icon) { return icon.id === id; }) || null;
}

// The CSS custom property that colours symbol `id` for a world of PHI-6 tier `tier` (0 to 3; ignored by a symbol with a
// token of its own, and `1`, Green, when a score-coloured symbol has no tier yet).
export function tokenOf(id, tier) {
  const icon = poiIcon(id);
  if (!icon) throw new Error("no point-of-interest icon named " + id);
  if (icon.token) return icon.token;
  const index = tier == null ? 1 : tier;
  if (!(index >= 0 && index < TIER_TOKENS.length) || Math.floor(index) !== index) throw new Error("no habitability tier " + tier);
  return TIER_TOKENS[index];
}

// The shapes of an icon file's text, for a canvas (which takes no colour from the page, and whose Content-Security-Policy
// takes no image from a data: or blob: URL): `{size, paths: [{d, evenodd}]}`, `size` the viewBox's side and each `d` the
// path data of a `<path>` for `new Path2D(d)`. The files are plain, one `<path>` per shape with the viewBox on the root.
export function shapesOf(text) {
  const box = /viewBox="0 0 (\d+(?:\.\d+)?) (\d+(?:\.\d+)?)"/.exec(String(text));
  if (!box || box[1] !== box[2]) throw new Error("a point-of-interest icon has a square viewBox from 0 0");
  const paths = [];
  const pattern = /<path\b([^>]*)>/g;
  let found = pattern.exec(text);
  while (found) {
    const d = /\bd="([^"]*)"/.exec(found[1]);
    if (d) paths.push({ d: d[1], evenodd: /fill-rule="evenodd"/.test(found[1]) });
    found = pattern.exec(text);
  }
  if (!paths.length) throw new Error("a point-of-interest icon has at least one path");
  return { size: parseFloat(box[1]), paths: paths };
}
