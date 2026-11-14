**Headline:** Five main test reds fixed and three stale scripts brought up to date.

### Changed

- TEST.740 (Lane 21 part): an oxygenic world may end under an ice-covered ocean (its oxygen rise burns off the methane greenhouse and the sea can freeze over, as Earth's did in the Huronian glaciation), so `test_only_oxygenic_life_holds_oxygen_it_made_and_only_in_a_sunlit_sea` allows it; `activity.SUN_LOGG` carries its measured tag; `test_a_radius_gets_the_mass_its_composition_needs` passes the rock's FeO share, not the whole composition; the materials rows `tholin` and `insoluble_organic_matter` lose their vocabulary-only tags because the tuning table reads them; `build_fudge_table.py` reads the quantities the class conditions use (`class_fit.SHOULDERS` is gone), `measure_sky_sphere.py` draws a sphere with `luminosity_sphere_batches` and `backfill_cells` (`backfill_sphere` is gone), and `class_census.py` no longer patches the removed age raise.

generation-output: unchanged
