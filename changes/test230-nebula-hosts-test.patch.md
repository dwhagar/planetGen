**Headline:** The O-star nebula test no longer fails on the hotter cluster stars that now live as long as real ones do.

### Changed
- **`test_o_stars_sit_in_h_ii_regions_and_cool_stars_light_nothing` (TEST.230).** It assumed the stars of a class D or E region begin with the nebula's own O5V. Since the stellar bundle (#2112, GEN.467) a star's life follows the sourced stellar tables, so at 4 million years the cluster behind the host still holds O3V and O5V stars (an O3V or O5V lives about 4 million years; the old law gave them under 1) and the host list may begin with a hotter type or count several O5V. The test now asks that the host is among the stars, and that a class C nebula (the host alone) lists exactly `O5V x1`. The generator is unchanged.

TODO items resolved: TEST.230. Added: none.
