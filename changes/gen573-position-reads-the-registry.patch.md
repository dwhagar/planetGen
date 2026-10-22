**Headline:** A point in space now takes its length unit by registry key and states no unit factor of its own.

### Changed
- **`SpatialPosition3D` takes `length_unit` ("m", "au", "ly", "pc"), a length key of the unit registry, instead of
  `length_unit_m`** (GEN.573). The factor comes from `unit_registry`; a key that is no length unit is refused. The point
  class loses its own `SPEED_OF_LIGHT_MS` (it reads `constants.SPEED_OF_LIGHT_M_S`), the bare 1000 between m/s and km/s
  (the registry's km/s factor, here and in `MassObject`), and `mu_of` and `mass_of` (a mass and its `mu` are the
  registry's declared "mass-mu" relation). Every stored and drawn value is bit-identical.
