**Headline:** The V8 layer of a planet's sulfur cloud now reads the cosmic dust its own system sends it, where before nothing told it whether any arrived.

### Changed
- **A rocky planet's cloud chemistry reads its dust supply (GEN.1542 slice 1).** `generation/cosmic_dust.apply` runs once the system's belts are final, right after the late accretion: for each rocky planet that holds a sulfur cloud it works out the dust the system's belts and comets send it (`influx.planet_influx_kg_per_yr`, GEN.1543, less the giants in the way), asks `dust_supply.v8_supplied` whether that is enough carbon (6.306 DUS1 and DUS2), and makes the cloud chemistry again with `carbon_dust` set to the answer (`sulfur_cloud_set.fill(..., carbon_dust=)`). V8 and V9 are single-source hypotheses and off by default (`sulfur_cloud_set.INCLUDE_HYPOTHESES`), so the pass does nothing while they are off and no stored value moves. With them on, a census of 150 default systems found 147 rocky bodies with some of the set and one with the acid cloud (a body at 0.2 AU, which also takes in 5e-8 kg of dust carbon per square meter a year, above the floor), and none that gains V8: the other gates of V8 (iron, chlorine) were not met there, so the dust is not what holds it closed.

generation-affects: none

generation-output: unchanged
