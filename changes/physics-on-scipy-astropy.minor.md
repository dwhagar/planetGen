### Changed
- Kepler's equation is solved with scipy (Brent's method, and a vectorised solver the comet orbit update uses for all comets at once); the hand-written Newton and bisection solvers are gone.
- Physical constants (G, c, the Boltzmann and Stefan-Boltzmann constants, solar, Earth and Jupiter figures, the AU, parsec and light-year) come from astropy. Solar mass is now 1.98841e30 kg (was 1.989e30) and solar luminosity 3.828e26 W (was 3.82e26), so derived star figures shift by up to 0.2%.
