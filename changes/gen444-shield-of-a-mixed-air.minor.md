**Headline:** a world's cosmic-ray and stellar-particle dose now follows the gases of its own air, not an air of Earth's make-up.

### Changed
- The dose a body takes under its air (`radiation.generate_dose`, GEN.444) counted every atmosphere as Earth air of the same mass. It now counts each gas by its own nuclear interaction length, worked from its elements (`1 / lambda = sum w / lambda_i`, the Particle Data Group's element lengths, which the rule reproduces for CO2, water, methane, ammonia and dry air to 0.2 percent), so a CO2 air shields 1.3 percent more per gram, a hydrogen air counts 1.7 times thicker per gram of column and a sulfur dioxide air 0.9 times; and a stellar proton is stopped by Bragg's rule over the same elements (`surface_dose.stopping_mixture`, `proton_range.energy_mev`). The reference is the rule's own dry air, so an N2 and O2 sky changes by under 0.5 percent. The mixed column that was written and read by nothing (`surface_dose.effective_air_column_g_cm2`) is now what the dose reads. Stand-ins SHL1 and SHL2 in section 6.291.

generation-output: changed
generation-affects: fill
