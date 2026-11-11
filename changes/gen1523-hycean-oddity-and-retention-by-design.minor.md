### Changed

**Headline:** Planets inside the Hycean box (1 to 10 Earth masses, 1.4 to 2.6 Earth radii, 150 to 500 K, 1.5 to 3 g/cc) are now named as Hycean candidates (oddity W26) and the page says why they may carry an ocean under a hydrogen envelope.

- GEN.1523: `physics/hycean_candidate.py` is now read by the generator. Oddity W26 is a naming row on a new lever L8 (a name on the finished body; no input moves), applied at the end of `calculate_atmospheric_conditions`, after the envelope share. The gate quantities `density_g_cm3` and `equilibrium_temperature_k` are now computed in `oddities.quantities_of` (the temperature from the light of every star at the body's distance, at the Earth's Bond albedo of 0.3). It draws nothing outside its own stream, so no other draw moves; the oddity ids stored on planets change, so the golden and the epoch are Lane 02's to regenerate.
- GEN.1523: `atmosphere_retention` and its rows SH1 to SH4 are marked by design in the wiring register: the pressure stage reads the Berta-Thompson logistic of the atmosphere recipes (the fit to observed planets), not the module's K-window fit.
