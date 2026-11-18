**Headline:** The first planet of an inner system keeps one orbital period for every star because the edge it sits outside is the co-rotation radius.

### Changed

- The notes on the first-planet period of an inner system (registry row PA19, `tuning.INNER_FIRST_PERIOD_MEDIAN_DAYS`, `planet_architecture.first_distance_au`) now say why it has no stellar-mass term: [R3196] puts the drop in planet occurrence at a distance that goes as the cube root of the stellar mass, the pre-main-sequence co-rotation radius, and that is a fixed period. A dust-sublimation edge would put F stars' first planets near 60 days and M2 dwarfs' near 2.5 days. A test pins the edge in AU to the cube root of the star's mass at one period. No generated value changes (GEN.1793).
