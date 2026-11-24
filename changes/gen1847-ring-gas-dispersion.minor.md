**Headline:** A barred galaxy's central gas ring now has its own gas velocity dispersion, drawn from the central gas of barred galaxies, not the Milky Way's 5 km/s for every ring.

### Changed
- GEN.1847: the ring of `galaxy/core_gas.py` took the Central Molecular Zone's 5 km/s for every galaxy's ring, and that dispersion sets both the marginally stable gas of the ring (GEN.1846) and the feedback time of the most massive cluster it can form (GEN.1195). It is now drawn with the galaxy's seed: log-normal with a median of 11.5 km/s and 0.30 dex ([R5360: Table 2], PHANGS-ALMA CO(2-1) of 30 barred galaxies without an AGN), held to two standard deviations (2.9 to 46 km/s). The sample is of the central 1.5 kpc at 150 pc resolution and finds no change with ring size, so every ring size takes the same law, and its bar's non-circular motions and beam smearing make it an upper bound on a ring's own turbulence (row RG1 of 6.334, tier T2, tweak L0). The Milky Way's 5 km/s stays the reference structure's test shape. The ring of a typical galaxy holds about 2.3 times the gas it did and forms heavier clusters at the top; the old constant `CORE_RING_GAS_DISPERSION_KMS` is gone.

generation-output: changed
generation-affects: plan
