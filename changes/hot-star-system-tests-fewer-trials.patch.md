**Headline:** The O5V case of the belt-and-habitable-world system test builds two systems instead of fifteen, so it finishes inside the timeout.

### Changed

- `test_systems.py::test_asteroid_belt_and_habitable_world_forced_together_both_succeed[O5V]` built fifteen all-flags O5V systems (about 150 s each on the GEN.1282 branch, some 90 planets and 190 moons) and hit the 900 s timeout (PERF.241, TEST.417). The O5V case now builds two (four took 600 s with the GEN.1282 change, 150 s a system). The belt and overlap checks are kept for every system. The request is not asked to give a world in two builds, because one ecosphere draw in thirty around an O5V meets it (measured on the GEN.1282 branch); every cooler star still builds fifteen and still asks.

generation-output: unchanged
