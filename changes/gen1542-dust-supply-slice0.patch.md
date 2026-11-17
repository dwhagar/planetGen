**Headline:** A planet's cosmic dust supply is worked out as a pure function of its system's belts and comets, ready for the V8 sulfur-cloud species to read.

### Added
- **The carbon a planet's dust supplies the V8 layer (GEN.1542 slice 0).** `physics/dust_supply.py` turns the dust a planet takes in each year (`influx.planet_influx_kg_per_yr`, GEN.1543) into the carbon it lays on each square meter (`carbon_flux_kg_m2_yr`) and says whether it reaches the floor the V8 layer of the sulfur-cloud set needs (`v8_supplied`). It is a pure function with tests and is not wired in: `sulfur_cloud_set` still takes `carbon_dust` from its caller, so no stored value moves. New registered stand-ins 6.306 DUS1 (carbon share of the dust, 0.1) and DUS2 (the floor, 1e-9 kg per square meter a year), both T6 and L4 with burn-down lines, and DUS3 (the dust a planet lifts itself is left out and why). Slice 1 will read each planet's influx at its fill and pass `carbon_dust`.

generation-affects: none

generation-output: unchanged
