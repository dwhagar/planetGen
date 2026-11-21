**Headline:** A moon history that no integrator finishes is no longer left as born: a moon the tide brings in within the age is lost to its planet, and the planet is marked.

### Fixed

- Over 233 generated systems (3,174 moon histories) two ended the `moon_orbits.STRATEGIES` ladder (0.06 percent), one with a `ZeroDivisionError` and one with an `OverflowError` inside the right-hand side. A trial step carried past the Roche limit or the Hill radius evaluated a distance of zero or of `exp(800)`; the right-hand side now holds that distance at the edge of what the moon can reach, so the step is rejected and the event found rather than the rung ended.
- What was left (a history that is still too stiff) used to keep every moon at its birth orbit with no account of the tide. `moon_orbits.inspiral_times_s` now gives each moon's `a / |da/dt|` at birth and `moon_history.apply` loses the moons that fall within the age of the system (recorded as `FallenMoon`s, the ring material), keeps the others as born, and sets `build.moon_history_unfinished` on the planet so that a census can count them (GEN.1816, TM8).

generation-output: changed
generation-affects: fill
