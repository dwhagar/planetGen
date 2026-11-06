**Headline:** The class-by-class sweep that the consistency tests read holds moons again, so the rules that read a moon have moons to read.

### Changed

- `tests/audit/sweep.py` (TEST.641): after GEN.962 a rocky planet draws a massive moon from a giant impact in 1 planet in 12 and its small moons only by capture, so the sweep's 126 planets (14 in each of three zones round each of three hosts) held 6 moons where they held 91 before it, and `test_the_sweep_holds_planets_of_many_classes_in_every_zone_and_some_moons`, which asks for 20 moons, found 7. The floor is right (the tidal, parent-shadow and moon-air rules need moons to read) and the draw is the one GEN.962 meant, so the sweep now also draws three planets in each zone round each host with four moons asked for (27 planets, 33 moons, each from a stream of its own, so the earlier 126 planets are the planets they were). The sweep holds 153 planets and 39 moons. No rule found a violation among the added bodies; the sweep tests that were red on main (`test_sweep_bulk` envelope helium, and the stand-in row tests) are red on it for the same reasons as before.

generation-output: unchanged
