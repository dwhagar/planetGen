**Headline:** A body's climate can now be worked out in seven bands of latitude (the equator, three bands on each side, the poles in the outermost): hours of day and night, the longest polar day and night, the day and night temperatures and the state of each; nothing stores or shows it yet.

### Added
- **Seven latitude bands (GEN.1108, part A).** A new `physics/latitude_bands` module is a pure function of the stored columns that gives, for each of seven equal-width bands, its share of the surface, mean temperature, day and night temperature, hours of day and night in a solar day, the longest polar day and night of the year, and a state word. Earth comes out at 12 hours of day at the equator, a polar day of 184 days at the pole (186 with refraction, which the model leaves out) and 179 days of night; Uranus has polar days and nights of years in every band; a world locked to its star gets seven rings from the day side to the night side. The registry has a new row (LB1). The storage part (columns, migration, the page) replaces the four bands `regional_temperature` stores and rides a later train.

generation-output: unchanged
