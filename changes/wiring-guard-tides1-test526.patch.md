**Headline:** The eight tide and surface-liquid functions that only tests called are wired, registered with the item that will call them, or deleted.

### Changed
- **TEST.526.** `orbit_elements` now takes the frequency of the eccentricity tide from `tidal_rheology.tidal_frequency_rad_s` (the orbital rate, the same number as before). `moon_orbits.tidal_power_w` and `surface_liquids.liquid_ids` had no caller and are deleted. `moon_orbits.drift_rate_m_s` is registered `by design` (the tests hold the Earth-Moon, Phobos and Triton recession to it), `surface_liquids.stable_liquids`, `ideal_sea` and `sea_density_kg_m3` are registered to GEN.1218 (the second liquid and mixed seas) and `habitability_world.harshness` to POP.38. None of the eight names TEST.520 in the wiring register now.

generation-output: unchanged

TODO items advanced: TEST.526.
