**Headline:** A planet's or moon's air, pressure and surface temperature now come from its interior, its sinks and its greenhouse together, not from a range set by its class.

### Changed
- **The air is solved, not drawn from a class range (GEN.411, GEN.422).** `physics/surface_air.py` repeats, until the pressure and the temperature stop moving: the mixture draw at the body's own temperature and pressure over the gas its interior has degassed (GEN.726), less what the carbonate cycle and an ocean take; the freeze-out of every gas above its saturation pressure, with water standing at a relative humidity of 0.77 over an ocean or ice; the hydrate cages of an icy cold body (GEN.762); the pressure as the weight of what is left; the density and scale height; and the greenhouse and temperature of that air. A giant, a brown dwarf or an envelope is its drawn mixture at the 1 bar level.
- **What this changes on a map.** Every rocky air, surface temperature, pressure, scale height and density is different, and the class no longer sets any of them: a world in the habitable zone can land on a thin cold air, a temperate one or a hot thick steam-and-carbon-dioxide one, as its interior and its sinks give.
- **The sulfur an ocean takes** is the solubility of GEN.757 and a drawn share for the sulfate and sulfide the ocean and crust lock up (a stand-in until a sulfur balance replaces it). The water-loss history (GEN.639) now reads the star distance before the air is made, so the surface water a young M dwarf stripped counts.

### Removed
- **The class density, molar density and greenhouse ranges (GEN.683).** `atm_density_range`, `atm_molar_density_range` and `greenhouse_multiplier_range` of every class in `tuning.py`, `CLASS_MIXES` and the mixing-ratio draw in `physics/atmosphere.py`, `ATMOSPHERE_DENSITY`, `ATMOSPHERIC_MOLAR_DENSITY` and `CO2_BASE_MOLAR_DENSITY` in `physics/constants.py`, and the tests and tuning options that only measured them (`test_climate_tuning.py`, the retention-factor and greenhouse-formula tests, the molar-density and greenhouse options of `climate_tuning_cli.py`). The class albedo ranges stay.

generation-output: changed
