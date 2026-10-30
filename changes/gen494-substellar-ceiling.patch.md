**Headline:** The hottest a brown dwarf can be at each age is now read from the Baraffe 2015 table (not yet used by the generator).

### Changed

- New library `physics/substellar_ceiling.py`: the temperature of the hydrogen-burning limit mass by age from the Baraffe et al. 2015 isochrones (about 2950 K at 10 Myr, 2050 K at 1 Gyr, 1630 K at 10 Gyr), with the age at which a temperature stops being possible for a brown dwarf. Nothing in generation reads it yet; the wiring step (the hot edge of the young late-M-like row and the 2800 K clip) follows with the Baraffe brown dwarf tracks.
