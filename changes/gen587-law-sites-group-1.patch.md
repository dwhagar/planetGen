**Headline:** Ten more places that wrote a physical law by hand now call the one declared law, with the same numbers.

### Changed
- The Schwarzschild radius (black holes, quasars, the phenomena plausibility check), surface gravity (planets, rogue planets, the plausibility range, the planet query), escape speed and mean density (collision outcomes) and the volume of a sphere (comet nuclei, gas clouds, pulsar planets) come from `physics/relations.py` instead of a typed formula (GEN.587, first file group). Every replaced expression gives the same float bit for bit; the sites whose operation order would change the last digit (`multistar._schwarzschild_au`, the `mathcheck` Schwarzschild helper, `lensing`) are left for the epoch-changing pass.
