**Headline:** The most water an icy body can hold now follows the disc's budget, where carbon monoxide and dioxide hold part of the oxygen, so it falls from 0.54 to 0.34 of the solids for the Sun (GEN.1363 part B).

### Changed

- `bulk_draw._icy_water_range` takes the top of an icy body's water from the new `disk_condensation.water_share_of_solids`, at the body's formation temperature and with its star's carbon and oxygen over its magnesium. For the Sun that is 0.34 beyond the water line and 0.23 beyond the CO line (Lodders's table gave 0.54 and 0.38); a star with more oxygen over its magnesium than the Sun's holds more water. The water of an icy body is a draw between a floor and this top, so every icy body's water ceiling falls by about a third for a Sun-like disc.
- `condensate_shares.water_share_of_solids` is deleted (nothing else read it); the jump of the solid surface density at the snow line stays on Lodders's table (`WATER_ICE_SOLID_BOOST`, 2.17) until the disc's solids are weighed against the minimum-mass nebula.
- With part C (the C/O of a giant) this is the third and last part of GEN.1363; the three share one epoch.

generation-output: changed
