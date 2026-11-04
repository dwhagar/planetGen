**Headline:** The resonant-chain tide test calls the tide functions that exist now.

### Changed
- **TEST.559.** `test_the_forced_eccentricity_feeds_the_tide_as_it_does_for_a_moon` called `tidal_heating.molten_power_w` and `solid_power_w`, which the tides chain (GEN.955 and GEN.778) replaced. Io's forced eccentricity now goes through the box's melt-limited tide (`interior_heat.heat_flow`, 1e13 to 5e13 W, the band `test_tidal_heating` holds) and TRAPPIST-1 d through `peale_power_w` with the rigid-sphere `k2` over a Q of 100, with its zero for a circular orbit. Test only.

generation-output: unchanged

TODO items advanced: TEST.559.
