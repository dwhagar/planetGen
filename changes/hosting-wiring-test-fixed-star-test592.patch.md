**Headline:** The directed-moon-count test builds its star from a fixed seed, so it no longer fails now and then.

### Changed

- `test_gen_hosting_wiring.py` built its `sun_like` star outside `draw.bound`, a different star every run, and `test_a_directed_moon_count_still_goes_through_the_gate` demanded that every planet asked for two moons keeps at least one. About 1 planet in 24 draws only moons its parent cannot hold (`MOON_REJECTION_LIMIT`) and keeps none, which the docstring of `generate_moons` allows ("fewer than requested"), so the test failed on some stars. The star is now drawn inside `draw.bound(1)`, and the test allows a planet with none while still demanding that no planet gets more than asked and that at least 80 percent get both. Twenty runs in a row passed. Tests only.

generation-output: unchanged
