**Headline:** A free-floating planet's shielding depth and radiation score now declare their true units, and a Wolf-Rayet star can carry bodies whose air is drawn.

### Fixed
- **Free planet field units (TEST.425).** The shielding column of a rogue planet was declared in a unit only the database column table knew, and its radiation score was declared in radians because its name ends in `rad`. The unit registry now has a column mass unit (kg/m² and g/cm², 10 kg/m² to the g/cm²), the column table reads it from there, and the radiation score is declared a plain ratio like the other scores. Values and stored columns are unchanged.
- **Wolf-Rayet hosts in the air recipes (TEST.425).** A Wolf-Rayet star (spectral type starting with W) had no host type in the mixture recipe tables, so drawing the air of a body around one raised an error. It takes the F rows, like the other stars hotter than the tables.
TODO items resolved: TEST.425 (TODO closes it after reading this). Added: none.
