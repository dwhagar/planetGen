**Headline:** One function now moves every body with a mass, whether it is a planet, a moon, a comet, a star or a wanderer.

### Changed
- `MassObject` (`generation/mass_object.py`) holds the state a body has relative to its parent (position and velocity, the pull of what it orbits) and `advance(seconds)`: a planet, moon or comet is stepped on its state vector with the Kepler step, and a star, remnant or wanderer on a circular galactic orbit has its phase advanced by its period. Boss's rule of 5:20 PM EDT: all motion is built into the one set of mass object functions.
- The `HoldsOrbitPosition` mixin in `physics/position.py` is gone; `Planet` and `Comet` take its methods (`set_position_au`, `set_velocity_kms`, `orbit_from_vector`, `projected_orbit_au`) from the core.
- No stored value or field changes. The per-class fields that become derived once the stored record is the core's (speed, axial tilt, rotation period, galactic period and phase) are removed in the DB.73 revision (Core 1), not here.
- `docs/design/mass-object-class.md` section 11 records it.

generation-output: unchanged

TODO items resolved: none (GEN.565 stays open) / Added: none
