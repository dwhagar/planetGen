**Headline:** The tests of one-off data rewrites in migrations Boss's server has already run are retired, with the old schemas they loaded.

### Removed

- The old-schema run of every data-backfill migration test in 23 `tests/test_migration_*.py` files and the same kind of test in `test_system_interests.py` (v164 to v175, 7 tests) and `test_life_row_text.py` (v155) (snow_line, brine_seep, star_members, star_count, habitable_worlds, parent_link, albedo_greenhouse, block_ladder_root, block_lights, life_scores, giant_interiors, gas_tables, planet_rings, region_stats, stage_epochs, sector_mass, life_classes, american_spelling, charted_designation, dust_species, isotope_record, generator_epoch_record, sector_sort_keys): the row rewrite, the "a second run changes nothing" repeat and the CHECK read on a migrated database. All of those migrations are at or below 0235 and ran on Boss's server when it updated (it is a test server, wiped when stalled, and carries no backward compatibility); the end shape is proved by `tests/test_db_old_schemas.py`.
- 11 of those files had no other test and are deleted; the others keep their round trips and their checks on a new database (`test_a_planet_cannot_hold_two_layers_at_one_index_or_an_unknown_kind`, `test_a_ring_kind_the_check_does_not_know_is_refused`, `test_a_stale_stage_is_found_by_its_own_epochs` and the rest). `test_migration_helpers.py` stays. The helpers other tests import (`ADDRESS`, `CENTER_PC`, `_make_species`, `_system`) stay.
- 14 old-schema fixtures nothing loads by number any more (v132, v133, v134, v139, v146, v160, v167, v174, v178, v180, v210, v221, v231, v234), and the duration line of the deleted `test_migration_region_stats.py` in `scripts/data/gate_durations.json`.

generation-output: unchanged
generation-affects: none

A `downgrade_to(config, N)` call loads the same stored schema as `load_old_schema`; TEST.729's check had not counted those callers, so it removed fixtures these tests still loaded and they went red on main. This PR retires those tests (they are data-backfill tests of migrations 0165 to 0176 and 0156, all run on Boss's server) and every fixture number a test still loads is present.
