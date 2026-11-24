**Headline:** The biggest spiral galaxies now sit in halos like the biggest spirals really do, ten times lighter than before.

### Changed
- GEN.618: a spiral of 1.5 to 2e11 solar masses of stars drew a halo of 1.4e13 to 2.8e13 (and up to 6e13 at the high end of the scatter), because the stellar-to-halo relation of Moster and others is fitted to all galaxies and at that mass its halos are a group's. The rotation curves of 110 spirals put the largest ones, 1 to 3e11 solar masses of stars, in halos of 2 to 5e12 ([R5352], read), and "virtually none" of the late types sit above 5e12 ([R5353], read); Andromeda's is 0.6 to 2.8e12. `mass_chain` now bends the spiral's median halo onto 2.5e12 (`SPIRAL_HALO_CEILING_SOLAR`, row HF5 of 6.131, tier T3, tweak L1) with a sharpness of 4 (HF6, tier T9, tweak L3) and squeezes the scatter above the median so no halo passes 5e12. Measured on this change, halo mass at the median (16th to 84th percentile): 5e10 solar masses of stars 1.57e12 to 1.51e12 (0.9 to 2.75e12), 1e11 5.7e12 to 2.5e12 (1.2 to 4.4e12), 2e11 2.8e13 to 2.5e12 (1.1 to 4.5e12). Lighter galaxies change by under 4 percent. Spiral is the only family built; a later family reads its own relation.
- The references list gains [R5352, R5353] (the two papers the ceiling rests on).

generation-output: changed
generation-affects: plan
