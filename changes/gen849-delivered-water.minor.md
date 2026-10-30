**Headline:** A rocky planet's own water now follows the same delivery law as its air: more water the farther from the star and the heavier it is.

### Changed
- **GEN.849: a rocky body's bulk water is drawn from the Bern delivery model (`physics/delivered_water`).** The body takes a water class (dry, low, high or rich) from the law the atmosphere recipes b1 to b4 already read, by its mass and the snow lines it formed at, then a log-uniform share inside the class (the rich class runs from 10 percent to a cap of 50 percent of its mass). It replaces the old 97 percent dry range and 3 percent of ocean-world draw, which did not look at the planet's mass or how near the snow line it formed. A rocky planet formed close to the star is now mostly dry and more of those formed near the snow line hold an ocean's worth or more. Bodies beyond the snow line keep their ice-and-rock range. `tuning.DRY_WATER_MASS_FRACTION` and `WATER_WORLD_*` are gone; the registry rows are CC3 (the icy range) and DW1 to DW3.

generation-output: changed
