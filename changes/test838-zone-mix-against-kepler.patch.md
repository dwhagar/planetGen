**Headline:** The zone mix test now holds the generated planets against Kepler's counts instead of the model's own old output

### Changed

- TEST.838 part 2 / GEN.1858: `test_zone_mix_is_spread` held the share of hot, ecosphere and cold planets of all stars to bands (hot 0.20 to 0.48) that recorded the model's own output of 11-05; GEN.1709 (Train 40, 11-16, sourced to [R2020]) moved the first inner-system planet from 5 to 25 days on purpose and the share fell to 0.18 to 0.22. The share cannot be observed (82 percent of a sample's stars are M dwarfs, whose hot zone lasts a few days, and the cold zone has no small-planet census), so the test now holds what Kepler measures: the planets of 1.25 to 4 Earth radii per FGK dwarf from 8 to 256 days, and the part of them inside the habitable zone's inner edge (the hot zone), against Hsu et al. 2019 Table 2 ([R2020], read from the paper's PDF) with the cells' errors and Poisson at the sample size, 3 standard errors of both. The zones still all occur. 0.79 +- 0.14 generated against 0.92 +- 0.09 observed (8 to 256 d), 0.69 hot against about 1.0. The 256 to 500 day bin runs low (0.14 against 0.52, 2 sigma) and is GEN.1871's. No product change.

generation-output: unchanged
