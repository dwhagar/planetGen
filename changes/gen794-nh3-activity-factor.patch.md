**Headline:** The ammonia active-Sun factor is computed from measured solar spectra instead of taken from an unchecked value.

### Changed
- The NH3 active-Sun factor of the photolysis table (GEN.794) is 1.07, computed by `scripts/nh3_activity_factor.py` from the Leiden NH3 photodissociation cross section and the daily solar irradiance of SORCE SOLSTICE, UARS SUSIM and UARS SOLSTICE (strongest year over weakest year of each record: 1.066, 1.083 and 1.060). It replaces the unchecked 1.078; the registry rows NA1 and NA2 give its tier, tweak level and bracket (1.05 to 1.15). No generated world changes: no generation code reads the active-Sun factor.
- The wiring register names GEN.850 for the acid cloud's column modules (GEN.1020 was retired), so the register check passes again.
