**Headline:** Eleven tests that were red on main are green again, and seven library modules no longer use maths whose last bit depends on the Python or BLAS version.

### Changed
- **Test expectations brought up to the generator (TEST.503).** `test_object_rows` (the star row holds `state` and `decay_age_yr`, which `Star.to_dict` does not carry; it reads seeds until every kind of body has been seen), `test_db_column_round_trip` (`moons.magnetopause_rp` is NULL where no moon has a field; `planets.porosity` is written now), `test_work_units` (the phenomena plan carries `outline_shells` and `min_mass_solar`; a luminosity fill ends with `refresh-block-lights`), `test_aquifers` (the life classes let an aquifer world hold animal grade A1), `test_climate_sweep` (the Earth twin is the first Earth-mass planet that stands at a temperate temperature, not a fixed seed), `test_gen_hosting_wiring` (a directed moon count keeps one moon where the orbits take no second).
- **Seven unread modules use `planetgen.util.detmath`** instead of `math.hypot`, `@` and `numpy.dot` (`rogue_density`, `band_winds`, `co2_far_wing`, `field_dispersal`, `mantle_stripping`, `planet_occurrence`, `quadruple_stability`); nothing reads them yet, so no stored number changes.
- A docstring in `scripts/bench_region_queries.py` names `src/planetgen/db/near.py` in full, which `test_deploy_scripts` needs.

generation-output: unchanged
