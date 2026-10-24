**Headline:** The hot-Jupiter chemistry regimes follow computed equilibrium lines at the pressure of the layer, not fixed temperatures.

### Changed
- A gas-phase equilibrium module (Gibbs minimisation with the NASA TM-4513 polynomials) computes the methane to carbon monoxide line and the water and hydrogen dissociation temperatures at a pressure and a metallicity; it agrees with an independent computation within 0.6 percent. Atmosphere recipe rows h1, h2 and h3 take their edges from it as functions of the photospheric pressure, in place of 1500 and 2200 K. The recipes are not yet read by the generator, so no generated value changes.
