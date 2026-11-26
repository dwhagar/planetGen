**Headline:** The planet physics tests' host star no longer changes with the test order

### Fixed

- TEST.839: `test_surface_temperature_never_drops_below_the_cosmic_background[20.0]` failed under `pytest-randomly` for some orders (seeds 17 and others on main). The cause was the module's `host_star` fixture: it drew the star outside any stream, so it took whatever the process run stream held when the first test of the file asked for it. A different age and luminosity changed the heat a 20 Earth-mass planet radiates (326 K for a young star against 43 to 85 K for an older one) and whether the draw is a giant ("g") or a smaller body with an envelope ("t"). The star is now drawn on its own named stream at the Sun's age and mass. No product change, no bound loosened, no retry.

generation-output: unchanged
