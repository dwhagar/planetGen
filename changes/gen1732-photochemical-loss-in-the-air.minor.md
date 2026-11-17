**Headline:** The air of a warm world now loses its ammonia and methane to starlight over the star's age, so the thick ammonia and methane airs of old worlds thin to what light leaves standing.

### Changed

- GEN.1732: `surface_air.solve` now applies the photochemical loss of NH3 (to N2) and CH4 (to CO2) after the hydrate hold, over the age of the star with the star's far ultraviolet, for every non-giant body, by a share that rises from 0 to 1 across 215 to 245 K of the surface temperature (a ramp, so the solved temperature has no step in it; stand-in row 6.301 PC5, tier T7, tweak L1). Reuse: `photochemical_loss.lose` (PR 6702), `surface_air` damped fixed point unchanged.
- The conversion between pressure and column in `photochemical_loss` uses the mean molar mass of the air (partial pressures are Dalton, proportional to moles), where the merged module used each gas's own molar mass. Moles are now conserved to the product.

generation-output: changed
