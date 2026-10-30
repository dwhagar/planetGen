**Headline:** A star type you ask for (G2V, K2III, B0V) is now a star the model makes, with the mass, brightness and age such a star really has.

### Changed
- **GEN.941: forced star types come from the star model.** A forced type (`+star_type=K2III`) was built from hand-set tables: a luminosity drawn where two range tables overlapped, a mass drawn evenly across a range per class and an age from a lifespan fraction. It is now drawn from the star model that makes random stars: the mass function and the time a star spends in each state decide which mass and age a type has, so a K2 giant is a 1 to 2 solar mass star on its giant branch, 8 to 11 billion years old, not a star of any mass in 0.8 to 8. A type the model never makes (a B giant, an M hypergiant) is the model's nearest star of that class by temperature, at the temperature you asked for. `IA+` and `D` are the hypergiant and white dwarf classes under their own names. A subdwarf (VI) is a metal-poor halo star at least 8 billion years old or, for the hot letters, a hot horizontal-branch star. Registry rows PS10 and PS11.
- The class reference pages for the spectral letters and luminosity classes read their mass and luminosity ranges, and a share of all stars, from the same model.

### Removed
- The luminosity, mass and spectral-share range tables the old forced-type and random-star branches read, the lifespan-fraction age draw, the white dwarf and subdwarf age windows, the evolved-mass resampler, the dead random branch of the star generator and the `mass_override` option nothing passed.
