**Headline:** A giant given a polar orbit at the final roll of the orbit rows is turned onto that plane with its position, spin and moons, not left lying off it.

### Fixed
- `planets.finish_orbit_oddities` rolls the orbit rows again once the orbit is final, and the polar orbit (W4, `orbit_oddities.polar_orbit`) writes a fresh inclination, after the secular exchange or the Kozai pump has moved the plane. With the eccentricity unmoved nothing recomputed the position, so a giant (a class I at 92.4 AU round an M50 hypergiant, where 92 AU is under 10 stellar radii) lay 40 AU off the plane its inclination and node name; `test_systems[M50-True-COMETS]` failed on that about 2 times in 270 on the Train 62 head.
- New `planets.turn_orbit_plane` takes the planet to the new plane the way the exchange and the pump do (same argument of latitude and periapsis, `release_moons`, `carry_spin_with_orbit`, `update_orbital_position`); `finish_orbit_oddities` calls it when the roll changed the inclination or the node. Reuse: the steps of `kozai_pump_planets` and `secular_exchange`. Regression test in `test_orbit_settle_vector.py` (fails without the fix: 0.07 AU off the plane).

generation-affects: fill
