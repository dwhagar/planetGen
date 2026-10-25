**Headline:** A brown dwarf's L, T or Y type is now drawn with the scatter of the 20 pc census instead of read from one mean temperature.

### Changed

- **Brown dwarf spectral types from the census map (GEN.760).** The generator draws every brown dwarf's type through the Kirkpatrick et al. 2021 temperature map (`physics/brown_dwarf_types.py`, GEN.517): the type is the map read at the dwarf's temperature plus the census scatter (134 K among L dwarfs, 79 K among T and Y dwarfs), so dwarfs of one temperature no longer all share one type, and the flat L/T transition spreads over several subtypes. The same type shows on a free brown dwarf's page and in a system that holds one. The older inversion of the mean relation (`tuning.BROWN_DWARF_SPECTRAL_TYPE_FITS`) is removed.

generation-output: changed
