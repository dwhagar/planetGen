**Headline:** A rocky rogue planet now turns at the rate its building impacts left it, and the mixture draw's default light is the orbit's mean light.

### Changed
- **Rogue rocky planets use the spin law (GEN.842, GEN.731).** A rocky free planet drew a day log-uniform over 8 to 48 hours; it now has the birth spin the giant impacts that built a planet leave (`spin.accretion_spin`: a rate scaled to its breakup rate and held under the breakup fraction, and an angle from the galactic pole that is isotropic), with no tide since there is no star. A rogue giant keeps the gas-giant day.
- **The mixture draw's default distance is the circle that gives the orbit's mean light (GEN.842, GEN.727).** `mixture_subject` took a planet's semi-major axis when no distance was handed in; the live generation path already hands in the mean-light circle, so only a direct call changes.

generation-output: changed
