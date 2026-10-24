**Headline:** The gas giant density test no longer fails on about one run in six, and the slow scatter-pass tests get a timeout that fits a loaded machine.

### Fixed
- **`test_gas_giant_density_is_its_mass_over_its_volume` failed in 7 of 40 runs, with or without load (the load flakes noted on TEST.344).** Since GEN.460 a giant's radius comes from the interior tables, and a hot young giant swells past 1.3 times the median mass-radius relation (one draw gave 46,695 km against a 32,741 km ceiling), so the test's 3-sigma band around the median was wrong. The test now checks the radius stays in the class's radius range and that the density is the mass over the volume; the 60 reruns after the change all pass.
- **`test_star_scatter_passes.py` runs the 900 s default timeout close on a loaded machine.** Its three scatter runs build the bright-star tables in pure Python (60 to 130 s alone), so the module now carries the same 3600 s timeout `test_budget_queries.py` uses.
TODO items resolved: none (these are the two load flakes noted on TEST.344, whose own item is the serial browser phase). Added: none.
