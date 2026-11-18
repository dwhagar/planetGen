**Headline:** hydrogen chloride now condenses on a liquid curve above its triple point, so it boils at one atmosphere where the gas table says it does.

### Fixed
- GEN.1630: the HCl vapour pressure above its triple point (159.07 K) was the solid's curve carried on, 36 percent high at the normal boiling point (137.6 kPa at 188.17 K) and 6.9 MPa at 273.15 K, and `condensation_temperature_k("HCl", 1 atm)` returned nothing. HCl now has a Wagner liquid curve in `saturation_pressure.json` whose four coefficients are fixed by values the gas table holds (the triple point, the normal boiling point, the acentric factor and the Edmister slope at the critical point; registry row AM104, tier 6, tweak L0). It gives 101.3 kPa at 188.17 K and 2.61 MPa (25.8 atm) at 0 C against the handbook's 25.6 atm, which the fit did not use.
- Stored values move only where HCl is in an air above 159 K: 12 of 1,226 bodies in 80 sample systems changed, by under 0.02 K and 0.1 percent of the pressure, no class or state changed.

generation-output: changed within the solvers' tolerance (no golden, no epoch, lock item)
generation-affects: fill
