**Headline:** A planet's or moon's ocean and ice are now as deep as compressed water really is: a world with a tenth of its mass in water has a 746 km layer, not the 1,047 km of a fixed density.

### Changed
- **The compressed water layer (GEN.461).** `planetgen.physics.water_layer` integrates the water over the covered area from the surface down: the density of water rises with the pressure through the liquid, ice VI, ice VII and ice X (the materials table's 300 K points), and the gravity at depth follows the mass below it. It matches the integrated reference of the design (part C 2.4) within 1 percent from a ten-thousandth to half of the mass. The liquid ends at the liquidus of the bottom temperature (the melting curves of ice III, V and VI, recomputed from the Journaux splines); under it the water is high-pressure ice and the ocean is sealed from the rock.
- **Real kilometres.** `ocean_depth_km`, `ice_shell_km` and `hp_ice_km` are the compressed thicknesses now, not water-equivalent ones. At Earth gravity the deepest liquid is 65 km at 280 K and 88 km at 300 K (it was 73 and 111), and a column is sealed from about 5e-3 of the mass. Rogue planets read the same module.
- **Stored (schema 151).** `planets` and `moons` gain `water_base_pressure_pa`, the pressure under the water (NULL for a gas giant, a dry body and a row made before this release).
- The old `max_liquid_depth_km`, `split_column` and their fixed densities (1050 and 1300 kg/m3) are gone.
- Generator epoch bump (numbered at merge): the ocean depth, the lid, the high-pressure ice and the ocean class of every watery planet and moon change, and with the ocean class the pH and water-activity draws that follow.

generation-output: changed
