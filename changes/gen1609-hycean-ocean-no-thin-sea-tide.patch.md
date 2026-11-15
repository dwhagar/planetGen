**Headline:** The reason a hycean ocean has no thin-sea tide is written down and a test pins it.

### Changed

- GEN.1609: `physics/ocean_tides.compute` gives a hycean ocean no tide because the formula and its loading factor are for a thin layer on the surface, and a hycean ocean lies hundreds of kilometres deep under a hydrogen envelope (the fluid Love numbers of GEN.1035 are what it needs). The module text says so and `test_a_hycean_ocean_has_no_thin_sea_tide_by_design` pins it. No stored value changes.

generation-output: unchanged
generation-affects: none
