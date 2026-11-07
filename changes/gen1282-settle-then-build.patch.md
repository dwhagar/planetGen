**Headline:** A system settles the spacing of its planets' orbits before any air is worked out, so a planet that is pushed outward is built once, and a crowded system of hot stars is several times faster to make.

### Changed

- **Orbits first, climate after (GEN.1282).** `StarSystem._generate_planets` builds each plain planet with `defer_climate=True`: the body and the shape of its orbit are drawn, but not its air, class, motion or moons. `validation.space_orbits` pushes planets outward by mutual Hill exclusion as before, and a deferred planet only has its zone, Hill sphere and period brought up to date (`reconcile_moved_planet`). `validate_system` then calls `Planet.finish_climate` on each planet once, at the distance it has. Before this a planet pushed by the spacing had its air, class, spin and moons worked out again for each push, and the all-flags O0V profile (PERF.241) spent 7 ms in the gas equilibrium solver for each of 29,800 solves. A system that bars habitable worlds from its ecosphere (`HABITABLE_WORLD` false) keeps the old order, since its redraw loop needs the finished body. Measured on O0V seed 6: 7 times faster.

generation-output: changed
