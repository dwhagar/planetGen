**Headline:** The reductants that eat a young world's first oxygen are now counted at the heat the world gave off then, so the oxygen of a planet rises later than before, as the Great Oxidation did on Earth.

### Changed

- The surplus that fills the store before the oxygen rises (`biogenic_oxygen.rise_age_gyr`) took the reductant flux `R` at the heat flow the body has now, for every age. The heat of a rocky body was 2 to 5 times higher when its life began, and `R` goes with it, so the surplus `S - R` is summed with `R` at the heat of each age (`reductant_history_factor`, from `interior_heat.heat_at`, the body's own thermal history, in steps of 0.25 Gyr) and past the star's age `R` is what it is now. The steady oxygen of the air still follows the balance now.
- `EARTH_STORE_MOL` is 4e21 mol, not 5e21: it is set so that Earth's delay after the 2.7 Ga biomarkers is still 0.714 Gyr with `R` at 2.2 times today's value falling to 1.9 times (the model's heat history of an Earth-mass body), inside the 0.3 to 1 Gyr bracket. Registry rows BO3 (value) and BO10 (the history factor) say so.
- Effect on 100 seeded systems of a G2V star with a habitable-world request (34 worlds with oxygenic life candidates): the median wait from the first producers to the oxygen rises from 6.3 to 8.9 Gyr, 25 of the 34 rise ages move, the worlds that never rise stay 9, and the worlds whose oxygen rises inside the star's age fall from 13 to 10. This is the physical direction (more reductant when young, so a later rise) and not the lever the filing hoped for: the shortfall of oxygen worlds is made a little larger, not smaller.

generation-output: changed
generation-affects: fill
