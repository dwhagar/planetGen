**Headline:** A new browser module works out how one star of the sky behind a system is drawn, ready for the system page's backdrop.

### Added

- `skylight.js` turns a star's apparent magnitude into the strength, width and colour of its point of light, with the same arithmetic written out in shader code for the sky layer (VIEW.33) to use. It builds on the Galaxy Map's star visibility law and smallest point size, so the sky and the map agree on what a faint star looks like.
- A star's colour is its blackbody colour (cool stars take the model's offset), read from a table of 128 temperatures, and fades toward white as the star dims. The limiting magnitude is a setting: 7.5 by default and 6.5 on a phone.
- Nothing on a page draws it yet. The constants (the peak knee, the point widths and the fade to white) are by-eye stand-ins, as the sky design says.
- Adds `tests/js/skylight.test.mjs`, which holds the colour to the Python colour code at fifteen temperatures and the shader text to the browser functions.

generation-output: unchanged
