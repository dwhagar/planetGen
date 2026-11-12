**Headline:** A star's light while its disk was there can now be worked out for any mass, with the contracting stage the main-sequence value leaves out; nothing reads it yet.

### Changed

- GEN.553 (birth step, for GEN.1611): new `physics/birth_luminosity.py` (`birth_luminosity_sol(initial_mass_sol, fe_h, age_gy)`, a pure function). From 0.01 to 1.4 solar masses it reads the isochrones of [R0254] at the age, so a Sun-like star at 10 Myr has 0.46 of a solar luminosity and a 0.3 solar-mass star 0.059 against 0.02 on the main sequence; above 1.4 solar masses it reads the zero-age main-sequence luminosity from `stellar_evolution.evolve_star`, so a star that is dead at 10 Myr (30 solar masses) is read as it began. The model's own contraction track is brighter than [R0254] at 10 Myr from 0.5 solar masses up (0.75 against 0.46 at 1 solar mass), so the two meet in a step at 1.4 solar masses (1.95 against 3.5 solar luminosities); the function does not hide it and the step is reported to the owner of the contraction track.
- Nothing in the generator reads it yet (wiring register): Lane 21 wires it into the formation snow line in GEN.1611 slice 1.

generation-output: unchanged
