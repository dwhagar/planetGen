**Headline:** A planet or moon is now its state vector, moved by Kepler's universal step, with no stored orbital phase.

### Changed
- **Planets and moons are stored as a point with a parent (GEN.289).** The nightly update steps each due body's position and velocity with `kepler.universal_step` (a compiled batch step, `physics/state_kernel.py`) from its own epoch, taking `mu` from the stored distance and period, and sets the next due time from the speed the new state has. `orbital_phase_deg` is gone from `planets` and `moons` (migration 0109, schema v109); the phase a page needs is read from the position. `insert_planet` and `insert_moon` now fill `parent_type` and `parent_id`, and the migration points a circumbinary planet at its system. The system scene API carries each body's `position_km`, `velocity_kms` and `epoch_unix` in place of `phase_deg`, and `orbitpositions.js` steps that state with a port of the same algorithm, so the browser and Python positions still agree. Tests hold the compiled step to the Python one, a quiet system to the cost of the old update, and the migration to a stopped run being started again.

TODO items resolved: GEN.289. Added: none.
