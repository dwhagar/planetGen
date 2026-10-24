**Headline:** The brown dwarf type windows take their edges and softness from the measured temperature of each spectral type.

### Changed
- Atmosphere recipe rows g1 to g5 use the type boundaries of the Kirkpatrick et al. 2021 temperature map (2238, 1273, 686 and 460 K with a scatter of 134 or 79 K) in place of the design edges. A new module `brown_dwarf_types` gives the temperature of a spectral type and the chance of each class at a temperature. The recipes are not yet read by the generator, so no generated value changes.
