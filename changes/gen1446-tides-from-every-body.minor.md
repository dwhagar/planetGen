**Headline:** the ocean tide of a planet or moon now counts every body that pulls on it: the companion stars, the other planets of its star and the other moons of its planet, and each pull is noted in the effect ledger.

### Changed

- `physics/tide_sources.py` (GEN.1446, first part): lists every source of tide on a planet or moon with the time-averaged `1 / d^3` it gives (exact two-body mean for the host and a moon's planet, an exact elliptic-integral mean for two circles, an orbit average for an eccentric source) and `ocean_tides` takes the sources it did not have as extra raisers. Spring and neap ranges move by the sibling's share, about 1e-5 to 1e-4 of the host's for a planet and more for a moon with a large sibling moon. Each pull is added to the effect ledger as a tide-raising acceleration. A lone star with no other body gives the numbers it gave before.
- Registry section 6.260 (TS1 to TS5): the coplanar upper bound, the planet with its moons as one point, the receiver on its circle, the synodic rate, and belts not yet sources (GEN.1464).

generation-output: changed
generation-affects: fill
