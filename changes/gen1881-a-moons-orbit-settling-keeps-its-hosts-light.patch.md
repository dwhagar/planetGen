**Headline:** A moon's orbit settling no longer sets its star exposure from its distance to its planet.

### Fixed

- GEN.1881: `planets.settle_orbit`, run for every planet and every moon at the star's final age, ended by setting the exposure (the star's XUV, flares and the stored star distance) of the body and its moons from `flux_distance_au(body.distance, ...)`. For a planet that is its distance from the star; for a moon, `distance` is the distance from its planet (0.000148 AU for Ekausor Ia, the sweep's seed 15), so a moon whose exposure existed already was given its own orbit as its star distance. Moons are normally built bare then (no star distance yet) and the line did nothing, but PERF.290 finishes a promoted belt planet's moons whole, and a system built with `HABITABLE_WORLD` off makes whole moons too: their stored star distance, exposure and the radiation read from it were wrong, and `test_radiation_sweep` found a moon at 17.63 K whose albedo and greenhouse at the stored distance gave 3,133 K. A moon now returns once its own eccentricity and place are settled: its light is its host's, set when the host's orbit moves (the loop that already covers the host's moons). `settle_warm_icy_orbit` reaches the same function and is covered. New test in `test_moon_frame`.

generation-output: changed
generation-affects: fill
