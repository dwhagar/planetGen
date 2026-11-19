**Headline:** The life-air tests no longer fail on a bigger fixture world: the fixed methane fluxes are scaled to the world's air, and the oxygen check is the ledger's own share.

### Changed

- `test_life_air.py`: the fixture world now holds an air ten times larger than the 3.1e20 mol the fixed methane fluxes were chosen for, so the methane step fell under 1 Pa and did nothing (7 tests). The fluxes are scaled to the fixture world's air (`_scaled`), as the history test already did.
- `test_the_generated_worlds_keep_one_air_and_methane_beside_life_leaves_no_free_oxygen`: a 608 kPa carbon-dioxide world holds 0.6 kPa of oxygen beside methane, 0.099 percent of its air and inside `oxygen_max_share` (0.1 percent). The check is now that share instead of 0.5 kPa.

generation-output: unchanged
generation-affects: none
