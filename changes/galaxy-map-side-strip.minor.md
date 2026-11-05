**Headline:** The Galaxy Map has an optional side strip that shows height above and below the galactic plane, linked to the map by selection.

### Added

- Menu, Side strip opens an edge-on strip under the map: the stars of a 60 degree wedge from the galactic centre (the wedge facing the map's focus) with the distance from the centre across and the height up and down, drawn 2.2 times taller than true so the thin disc can be read. Dashed lines mark one scale height of the thin and of the thick disc either side of the plane.
- Pick a star on the map and it is ringed in the strip, with its exact height and how many thin-disc and thick-disc scale heights that is. Pick a star in the strip and it is selected on the map, with its panel.
- On a phone or a touch screen the strip is a drawer fixed to the bottom of the window, drawn with fewer stars (2,000), and its Close button puts it away.
- Adds `galaxystrip.js`, `tests/js/galaxystrip.test.mjs` and `tests/test_web_browser_side_strip.py`.

generation-output: unchanged
