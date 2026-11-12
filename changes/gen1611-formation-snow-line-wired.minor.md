**Headline:** The line beyond which a planetary system's ices formed now follows the light the star's disk actually saw: white dwarfs and giants get the line of the star they were born as (tens of AU for a white dwarf, not 0.02 AU), red dwarfs get a line two to three times farther out than their present dim light gave (they were brighter while their disks were there), and the spread of the line is cut at the observed factor of three instead of four.

### Changed

- GEN.1611: a star's stored `snow_line_au` is drawn around the blackbody line of the luminosity its disk saw: its initial mass at 10 Myr (`physics/birth_luminosity.py`, the BHAC15 isochrones to 1.4 solar masses and the zero-age main sequence above), not its luminosity today. Everything that asks "formed beyond the snow line" (planet and moon ice shares, giant ices, belt type and origin, the disk edge at 18 snow lines) reads it. A white dwarf born at 4.36 solar masses stores about 40 AU where it stored 0.02 AU; a K giant born at 1.04 solar masses about 2 AU where it stored 22 AU. On 120 seeded systems the median line of a red dwarf (under 0.6 solar masses) moves from 0.16 to 0.42 AU, a Sun-like star's from 1.42 to 1.34 AU, and the share of an M dwarf's planets born beyond it falls from 0.36 to 0.29.
- `tuning.SNOW_LINE_SCATTER_MAX_SIGMA` is now the sigma of a factor of three (2.39) instead of three sigma (a factor of four, which put the line inside the habitable zone for 104 of 3,000 stars): the sourced spread is a factor of 2 to 3 ([R3193]). `tuning.SNOW_LINE_SPREAD_FACTOR` holds the factor (registry PA19).
- Compact remnants that are not planet hosts (black holes, neutron stars) keep a line from their own luminosity; no progenitor mass is stored for them.

generation-output: changed
generation-affects: fill
