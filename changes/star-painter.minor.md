**Headline:** Stars on the System Map are painted from their class: their true blackbody colour, a surface that matches the kind of star, and a glow that grows with luminosity.

### Added
- **The star painter (MAP.276).** `static/starpaint.js` paints a star's disc from its look record: O, B and A stars are nearly featureless, F and G stars finely granulated with a few dark spots, K and M stars coarser with more and larger spots (an M dwarf can show a bright flare patch), giants and supergiants show a few huge convection cells, and a white dwarf is one smooth disc. The colour is the blackbody colour of the star's effective temperature (`photometry.star_colour`, the one the sky-colour model assumes), so a Sun-like star now reads white with a faint warm tint rather than the spectral-letter yellow. The glow around a star is sized and strengthened by its luminosity (capped), where every star used to get the same one.

### Changed
- **The star look record carries the Yerkes class.** The system query (`get_system_detail`) now reads `stars.yerkes_class`, and `bodylook.star_look` puts it in `data-look` as `yerkes`.
