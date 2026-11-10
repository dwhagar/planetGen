**Headline:** The belt-and-habitable-world system test builds two B2V systems, not fifteen, as it does for the other heavy stars.

### Changed

- `tests/test_systems.py::test_asteroid_belt_and_habitable_world_forced_together_both_succeed[B2V]` built fifteen all-flags systems at about 48 s each (716 s measured on a lane machine with four other jobs running) and ran into the 900 s limit on a slower one. B2V joins `HOT_STARS` (O5V, M3IA+, M50, PERF.241 and PERF.249 took the same step), so it builds `HOT_STAR_TRIALS` = 2 systems and keeps the belt, overlap and request checks for each. The generator is unchanged.

generation-affects: none
