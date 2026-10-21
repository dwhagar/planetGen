**Headline:** Every body with a mass can now report its whole kinematic record: where it is, how it moves, how it spins, its mass and when it is next due.

### Changed
- `MassObject.kinematic_record()` returns a `KinematicRecord` (position and velocity relative to the primary, spin axis and rate in rad/s, mass and mu, epoch and next-due), the record Boss's S16 card gives every mass object. `speed_kms()`, `axial_tilt_from_state()` and `spin_rate_rad_s` read what is stored today as speed, axial tilt and rotation period from the state and the spin, and equal the stored values; a non-rotating black hole's rate is zero.
- `spin.orbit_normal_from_state` gives the orbit normal `r x v`.
- No stored value or field changes; the per-class columns that become derived are dropped in the DB.73 revision (Core 1). `docs/design/mass-object-class.md` section 12 records the fields.

generation-output: unchanged

TODO items resolved: none (GEN.565 stays open for the family classes) / Added: none
