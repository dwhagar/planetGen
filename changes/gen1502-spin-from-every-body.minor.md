**Headline:** the tidal locking and spin-down of a planet or moon now count every body that pulls on it, not only its host and its own moons.

### Changed

- `tide_sources.register` tells each planet and moon who else pulls on it (the other planets of its host, and a moon's planet) when the final-age spin history is run; `body_spin.torques_of` adds the tide of every companion star and sibling planet (and, for a moon, the host star at its planet's place and the other moons of its planet) as extra `spin_torques.Raiser`s (GEN.1502, GEN.1446 part 2). The first raiser stays the host, so spin-orbit resonances are measured against it as before. Over 40 systems (148 bodies) one rotation period moves, by 5e-9 of itself: the extra pulls are far below the host's.
- Registry row TS6: the torque of such a source uses the square of its mean `1 / d^3`.

generation-affects: fill
