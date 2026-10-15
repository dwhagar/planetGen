**Headline:** The tests that are not about the Galaxy Map now check the spherical sector grid

### Changed
- Nine tests that only needed an address, a position or a neighbor now use `galaxy/spherical.py` instead of the old cylinder `galaxy/geometry.py`.
- The sector-exhaustion, galaxy-boundary and line-sector tests are rewritten for the shell, band and wedge cells, and the grid fuzz tests (`test_fuzz_galaxy_grid.py`) check layout, positions, cells, neighbors (against a brute-force scan), designations, enumeration and `describe_sector_cell` for the new grid.
- The fuzz tests found that a sector designation cannot be made for a cap or high band past about shell 1,300 (the band field is 12 bits); the tests pin that as a clean `ValueError`.
- The Galaxy Map, drill, viewport and prism tests stay on the old ladder until the Maps rebuild replaces it.
