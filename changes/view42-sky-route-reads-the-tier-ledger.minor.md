**Headline:** The system sky page now knows when charting has built its stars: "built to" follows the spheres charting drew instead of staying at zero.

### Changed
- The sky route (VIEW.31) took "built to X pc" from the cells' luminosity levels alone, which stop at the 50 Lsun of the luminosity bands, so after a full chart it said 0.0 and "filling" stayed on. `sky_charting.built_to_pc` (VIEW.42) now reads the VIEW.25 `sky_tier_draws` ledger: when the spheres that hold the observer reach from the scatter's luminosity level down to 50 Lsun and from the mass limit down to the lightest ring, the sky is built out to the largest of those spheres' margins; otherwise it is 0.0. The route reports the larger of that and the cell-level figure.

generation-output: unchanged
