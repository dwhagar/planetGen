**Headline:** Whether a body formed beyond the snow line, and so is icy, now reads the same viscous-disk snow line as its water.

### Changed
- **One snow line for the water and the ices (GEN.1679).** `bulk_draw.origin_of` (icy or rocky make-up) and `planets.formed_beyond_snow_line` (the stored `born_beyond_snow_line` flag the rings and the make-up read) used the luminosity line while the water of a rocky body (GEN.1678) read the viscous-disk line, so a body could be icy by origin and dry by the water law. Both now read `formation_snow_line.water_snow_line_au`. The planet chain, the belts and the comets still read the luminosity line.

generation-output: changed
generation-affects: fill
