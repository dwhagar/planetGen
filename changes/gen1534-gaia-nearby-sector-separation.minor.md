**Headline:** Systems in the galaxy may now lie as close together as real stars do (30 percent of real systems have a neighbour within a parsec), because the placement rule that forbade it is fitted to the Gaia stars within 25 parsecs.

### Changed

- GEN.1534: new setup dataset `gaia_nearby` (the Gaia DR3 sources within 25 pc whose parallax is measured to 10 sigma, 5,323 stars, 767 kB, fetched with the others into the cache folder by `update.sh`; `required` False, only the fit and a test read it). `datasets/gaia_nearby.py` merges the stars of a system (within 0.1 pc) and gives the nearest-neighbour distribution and density of the 958 systems within 15 pc. Licence line added to DATA-LICENSE.md.
- `tuning.SECTOR_HILL_SEPARATION_SHARE` (stand-in SEP1) is now 0.05, was 0.5. `scripts/fit_sector_separation.py` places systems at the sample's density by the galaxy's own rejection and compares nearest-neighbour distributions: the old half of the sum of two Hill radii (1.3 pc between two median systems) allowed no neighbour under 1 pc and moved the median nearest neighbour from the real 1.33 pc to 1.65 pc; shares up to 0.075 are inside the sample's 95 percent band, 0.05 is the middle. Registry SEP1 goes from tier T9, tweak L3 to T1 for the statistic, L2; burn-down line set to replaced.
- Systems the old share rejected are now accepted, so the random draws of placement shift and the galaxy's layout changes with the epoch (no system is dropped: a smaller share has more room).

generation-output: changed
generation-affects: scatter
