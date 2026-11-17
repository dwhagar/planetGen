**Headline:** `test_tidal_spin` makes its seeded systems four at a time, so the file runs in a third of the time.

### Fixed

- TEST.788: `tests/test_tidal_spin.py` made the whole systems of three seeded sweeps one after another (the giant-tilt search up to 300 seeds, 12 seeds of the tide-state check, 24 default systems) and took 1,583 s alone, 1,662 s in a gate. The systems are now made by forked workers, each on its own seed, and the giant-tilt search reads them in seed order and stops at the same seed, so every sample and every bound is the one the one-by-one run used.
