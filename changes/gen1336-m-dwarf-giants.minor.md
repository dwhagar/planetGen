**Headline:** Giant planets are now rarer round lighter stars, so M-dwarf hot Jupiters come at 0.3 percent of stars and not 0.5 percent (GEN.1336).

### Changed

- `planet_mass.gas_outcome`: the share of runaway-capable cores that find gas is also multiplied by (star mass / the Sun's mass) for a star lighter than the Sun (`tuning.GIANT_GAS_SURVIVAL_STAR_MASS_EXPONENT` = 1.0, the survey's linear law, [R2172]; a heavier star is not raised). The disk's gas goes with the star's mass, and an M dwarf's snow line lies so close in that its cores reached the critical mass as easily as a Sun-like star's.
- Measured: stars of 0.55 solar masses with a giant within 2.5 AU 7.2 to 4.0 percent (survey 3.5 percent at 0.5 solar masses); hot Jupiters 0.5 to 0.32 percent over 12 seeds of 2000 stars (survey 0.27 +- 0.09 percent, [R1504]). Sun-like stars are unchanged.
- The test of the hot Jupiter rates now pools three seed blocks for the red dwarf (it passed on seed 405 only). Stand-in row PA2 of the registry gains the new constant.
- Every giant of a star lighter than the Sun can change, so the planets of those systems change (epoch).
