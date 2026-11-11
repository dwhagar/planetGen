**Headline:** the sphere a planet's moons must stay inside now counts the tide of the companion stars and the other planets of the system, not only the star's, so a planet with a heavy neighbour or a close companion holds its moons a little closer.

### Changed

- `physics/tide_sources.py` (GEN.1449): `perturber_ratio` is the sum of `m g3` of every source besides the host (the companion stars, the sibling planets, the belts of the list) over the host's `M / a^3`, the same means the ocean tide and the ledger use. `planets.update_hill_sphere` shrinks `hill_radius` by `(1 + ratio)^(-1/3)`; `min_orbit_distance` (the clearance the next planet keeps) stays the star's Hill radius alone. `StarSystem.validate_system` sets the ratio on each planet before its moons are made, so the moon limit (`moon_bounds_km`, the escape in `moon_orbits.evolve`) uses it.
- Registry TS8 (the sum of tides in the Hill balance, T7 L1) and TS9 (moon-moon resonance not built) in 6.260.
- Census, 719 planets of 80 systems (253 in binaries), with and without the ratio: the sphere shrinks by a median of 5e-5 and at most 2.1e-3 (7 planets over 1e-3, none over 1 percent). Moons come out different in 228 of 330 planets that have moons by a median of 1e-5 of their distance (4 over 1e-3, none over 10 percent), because the limit scales the draws of their orbits; 2 planets end with a different number of moons.

generation-affects: fill
