**Headline:** a planet or a moon can now be passed to a new parent without moving in space: it keeps its absolute position and velocity, its orbit about the new parent is worked out from them, and its climate, class and spin are made again for the new place. Nothing in a generated galaxy changes yet.

### Changed

- `physics/reparent.py` (GEN.583, reused: `state_vectors.elements_from_state`, `kepler.gravitational_parameter_au3_yr2`): the new relative state of a body about a new parent, the new pull by the generator's convention, whether the orbit is still an ellipse, and the centre of mass, linear momentum and angular momentum that the change must not alter.
- `generation/reparent.py`: `reparent_planet` and `reparent_moon` apply it and make again only what the new orbit changes, through `validation.reconcile_moved_planet` (the code GEN.1597 fixed to keep the spin a planet was born with) and the new `validation.reconcile_moved_moon` (the moon loop of that function, unchanged in behaviour); a body that would not be bound is raised as `Freed` before anything is changed. No generation path calls them yet: the unbind, capture and merger events (GEN.291) and the orbital update of 8.3 will; the wiring register says so.
- Tests (TEST.244 folded in): absolute position and velocity kept, the period is the generator's Kepler law, there and back returns the state, a moon's velocity is the sum of its hops, a body that would leave is reported, barycentre, momentum and angular momentum are kept, a generated planet and moon keep what they are, and a locked moon's day follows its new planet.

generation-affects: none
