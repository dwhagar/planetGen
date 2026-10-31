**Headline:** Every planet, moon and comet now keeps how round it is, how strong its material is and how tightly it is bound, and every belt and asteroid field keeps the mean density and ice share of its members.

### Added
- `physics/body_strength.py` (GEN.1159, GEN.1178): from a body's mass, radius, density and composition it works out `roundness_ratio` (the mean radius over the round radius, 1 or more is round), `yield_strength_pa` (the whole-body strength of its material mix) and `binding_energy_j` (two uniform shells for a rocky body with a core, the n = 1 polytrope for a gas giant), and sets them when a planet, a moon or a comet is made. A gas giant keeps only the binding energy. The strength is drawn from the range each material has (ice 0.4 to 3 MPa, rock 4 to 30, metal 50 to 250) by a quantile from a stream of the body's own (its mass, radius, distance and whether it is a moon), so the run's stream does not move and the line between round and not round is fuzzy as it is in the Solar System. The strength is the cold strength: the engine's temperature is the one a body took its shape at, and no body stores that yet.
- Schema v190 (migration 0190): `roundness_ratio`, `yield_strength_pa` and `binding_energy_j` on `planets`, `moons` and `comets`; `density_g_cm3` and `ice_mass_fraction` on `comets`, `asteroid_belts` and `asteroid_fields` (the belt and field columns are the mean of the members, from the drawn makeup). The property sheet shows them; the API gains optional fields and its version stays. Rows made before v190 get them when the galaxy is regenerated.

### Fixed
- The olivine and ice creep laws returned an overflow error at a few kelvin; the stress is now infinite there, so cold rock keeps its whole-body strength.

generation-output: changed
