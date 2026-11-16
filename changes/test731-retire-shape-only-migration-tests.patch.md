**Headline:** 28 migration tests that only checked a column was added are retired, and the 14 stored old schemas they alone loaded go with them.

### Removed

- The `downgrade_to(N)` or `load_old_schema(N)` then migrate test of 28 `tests/test_migration_*.py` files (ozone_column, survival_equipment, regional_temperature, aerosols, radiogenic_content, bulk_composition, water_base_pressure, moon_tidal_flux, interior_heat, aquifers, rogue_habitability, surface_composition, atmosphere_colour, lid_regime, day_night_shares, rogue_rings, ice_species, phenomenon_sort_indexes, surface_liquid, ocean_tides, ring_age, dose_belt, surface_liquid_solutes, nuclear_arrival, surface_liquid_mix, temperature_swings, temperature_coverage, moon_origin). Each asserted only that a column exists, has a type or reads NULL after the migration. `tests/test_db_old_schemas.py` compares the shape (type, nullability, default, index, keys) of the baseline and the newest three old schemas migrated with a new database, and `tests/test_db_column_round_trip.py` (TEST.11) proves every column is stored and read.
- 12 of those files had no other test and are deleted. Each other file keeps its save-and-load round trip, and `test_a_stopped_run_can_be_started_again` stays in survival_equipment and regional_temperature (it proves a half-applied migration finishes, which is not a shape check).
- `test_migration_snow_line.py` is not in this set: its v131 test computes backfilled values, so it goes with TEST.733.
- 14 old-schema fixtures nothing loads by number any more (v137, v152, v168, v169, v177, v179, v185, v187, v189, v190, v194, v196, v204, v216).

generation-output: unchanged
generation-affects: none
