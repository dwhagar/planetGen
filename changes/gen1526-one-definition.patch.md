**Headline:** Physical constants that several modules each defined for themselves (Avogadro, the gas constant, the densities of ice and water, the melting point of water, the latitude band count) are now defined once in the constants module.

### Changed

- GEN.1526 (from the GEN.945 audit): `AVOGADRO`, `GAS_CONSTANT`, `ICE_DENSITY_KG_M3`, `LIQUID_WATER_DENSITY_KG_M3`, `WATER_MELTING_K` and `BAND_COUNT` are defined in `physics/constants.py` and the modules that each held a copy (isotope record, air column, elements, interior heat, species thermo, core heat, tidal rheology, rings, clathrate brine, aquifers, hydrosphere, regional temperature, planet classes, latitude bands, class domains) read them from there under the same names. The values are the ones the modules already held, to the last digit (`GAS_CONSTANT` keeps its ten-digit CODATA value; astropy's full-precision `R` stays a separate name), so nothing generated changes. `test_astropy_constants` now holds `AVOGADRO`, `GAS_CONSTANT` and `WATER_MELTING_K` against astropy and fails when a module gives one of the six names a number of its own again.

generation-output: unchanged
generation-affects: none
