**Headline:** Every map is drawn on black in either theme: the Galaxy Map, the sector and System Maps (diagram and 3D), the Phenomena and navigation maps, the rendered picture of a world, and the nebula and phenomenon views.

### Changed
- **Black map backgrounds (Boss, 2026-10-31).** A new `--map-bg` token (black) is the background of every map viewport (`static/style.css`), where the maps used the page's themed panel colour, and two near-black navy fills (the 3D system view, the phenomenon picture) and two scene backgrounds (`nebulaview.js`, `phenomenonrender.js`) are black too. The rendered picture of a world on its page sits on black.
- **Readable on black.** The map viewports always take the dark theme's colours (text, controls, accent), so labels and controls are light on black in the light theme as well. Ring, axis, compass and orbit lines use a new `--map-line` token that is 3:1 or more against black. The scripts no longer treat a map as light (`mapcore.isLightBackground()` is always false) and read the map's colours from a viewport (`mapcore.mapVar`).

Display only: no generator output, schema or API change.

generation-output: unchanged

TODO items resolved: none. Follow-up for TODO to file: remove the light-theme branches in the map scripts that this leaves unused (`isLightBackground` and its callers, `nebulalook.js` light looks, the glow layer's `light` option, `sectorscene.js`'s `lightBackground` uniform).
