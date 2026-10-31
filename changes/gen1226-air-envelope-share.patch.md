**Headline:** A rocky planet with a hydrogen-rich air now stores the true mass share of that air as its gas envelope, where it stored zero.

### Changed
- **The envelope share of a hydrogen-air body (GEN.1226).** `planetgen.physics.air_envelope` gives the air's mass over the body's mass for a rocky body whose air is mostly hydrogen and helium; `planets.py` stores it in `envelope_mass_fraction` (any other rocky body keeps 0). The class T "gas dwarfs" hold 1e-7 to 1e-3 of their mass in air, so the body type of GEN.1219 still types them terrestrial; the number is true and the type reads it.
- generation-output: changed
