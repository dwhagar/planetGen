### Added
- `SpatialPosition3D` keeps a velocity in the galactic frame and relative to the nearest star (the system frame), the star's own velocity, and an `epoch_unix`: when the position and velocity hold. A body's speed relative to its star is read separately from the star's speed round the galaxy.
- `planetgen.physics.state_vectors`: `state_from_elements` and `elements_from_state` convert between a position and velocity and orbital elements for any ellipse, parabola or hyperbola, with `mean_anomaly_from_true` and `closed_orbit_points` for drawing the ellipse; `orbits.circular_orbital_velocity_au_per_year` gives the velocity on the circular orbits planets and moons follow.

### Changed
- `SpatialPosition3D.get_time_to_observable_movement` measures a system-scale or planetary move by the speed relative to the star, and a galactic one by the galactic speed.
