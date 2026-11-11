**Headline:** Tests that need a rare case search a bounded number of seeds for it instead of pinning the seed that once had it.

### Changed

- TEST.714: new `tests/fuzz_support.first_seed_with(build, wanted, seeds=range(64), what=...)`: builds each seed under `deterministic_entropy`, stops at the first that `wanted` accepts, and fails naming how many seeds it tried when none has the case (with its own tests in `tests/test_first_seed_with.py`).
- `test_describe_render.py::test_headings_follow_the_tree_in_each_format` searches the seeds for a described system with a moon (it took the first of three pinned seeds). `test_moon_history.py` finds the planet with its born moons (two tests) and the moon that moved by 5 percent (it skipped when none of 16 seeds had one) by the same search; the sample test runs its 16 seeds as before and goes on to 64 only while it has not seen a fallen moon and six moons.
- The Train 22 pin reds the Train Manager listed (`test_planet_zones`, `test_circumbinary_light`, `test_db_column_round_trip`, `test_species_limits`/`test_population`, `test_swept_hand_on`) are on `claude/train-22`, fixed there by their owners; they are converted to the helper once they are on main.

generation-affects: none
