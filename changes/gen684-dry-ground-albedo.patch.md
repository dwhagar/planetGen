**Headline:** A dry rocky world's Bond albedo now comes from its ground, so airless rock is as dark as the Moon and Mercury instead of twice as bright.

### Changed
- **Dry bodies take their Bond albedo from their ground (GEN.684 part A).** A rocky body with less than `tuning.ALBEDO_WET_MIN_WATER_FRACTION` of water (or no starlight) used to draw its Bond albedo uniformly from 0.12 to 0.35 before its air or ground existed. It now takes `ground_albedo.fresh_bond_albedo`: the mean visible reflectance of the dominant dry ground in the materials table (basalt powder 0.215; carbon rock 0.018; a stripped metal body the 0.15 stand-in), and a body that ends with no air is darkened by space weathering (`ground_albedo.weathered_albedo`: tau_w 1 to 2 Gyr, d 0.3 to 0.7 of the body's age). On 40 census systems the drawn airless dry rocks now sit at 0.075 to 0.143 (the Moon is 0.11, Mercury 0.088) where they were drawn at 0.12 to 0.35, so their equilibrium temperature rises by about 8 percent (25 K at 280 K), and classes move with it. Wet and frozen bodies, giants and bodies with an air keep their albedo laws. Registry rows GA1 to GA3 (section 6.287).
- **`ground_albedo.band_albedo`, `row_albedo` and `BandAlbedo` are retired**: `surface_composition.band_albedo` is the same area-weighted mean and is what the colour reads; `ground_albedo` keeps the weathering law and gains the dry-ground albedo.

generation-output: changed
generation-affects: fill
