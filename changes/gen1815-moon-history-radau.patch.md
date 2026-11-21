**Headline:** A planet with a big moon close in no longer gets its tide heat from a moon history that failed to integrate: the moon-orbit ladder has a fourth integrator that finishes the stiff cases.

### Fixed

- `moon_orbits.evolve` tried LSODA and then three BDF runs and raised `IntegrationFailed` when none finished; `moon_history.apply` then kept the moons as born and spun the planet without them. A rocky planet of 8.9e24 kg with a 5.0e23 kg moon born at 3.8 radii over a star-locked spin (seed 9101 of the aquifer tests) failed all four, and the moon, still inside the corotation radius of a spin that had never felt it, put 2.5e7 W/m2 of tide heat into the planet at every age (Lane 02's 7306 then asked for 6.5e10 flood basalts and ran Train 51 out of memory). `STRATEGIES` now ends with Radau (implicit, order 5, RTOL, up to 60,000 evaluations), which finishes that history in 1.4 s: the moon falls to its planet in 1.4 Gyr, and the planet's tide flux is 8e-08 W/m2. Systems whose history the first four integrators finished are untouched (GEN.1815).

generation-output: changed
generation-affects: fill
