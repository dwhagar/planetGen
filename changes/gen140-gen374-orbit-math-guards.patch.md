**Headline:** The orbit stepper and the browser's orbit positions keep their digits over long times, and Python and JavaScript agree.

### Changed
- **Stumpff functions by series for `|z| < 1`** (GEN.140): `kepler.stumpff_c2_c3`, the compiled `state_kernel._stumpff` and `stumpff` in `orbitpositions.js` sum ten terms of the series below `|z| = 1` (it was `1e-3`), where the closed forms cancel and a hyperbolic `c3` was off by 4e-13. Body positions the universal step gives on the maps change in the thirteenth digit and below; nothing stored by a generate changes (the step runs on read, `physics/body_positions.py`).
- **Orbit phases reduce first** (GEN.374): `kepler.wrapped_turn_deg(elapsed, period)` is `360 fmod(elapsed, period) / period`, used for circular bodies and elliptical comets (`physics/body_positions.py`), the nodes of a multi-star system (`physics/system_tree.py`) and `wrappedTurnDeg` in `orbitpositions.js`; the plain ratio is off by about 1e-4 degrees for a 0.3 day period at 1e9 days. The JavaScript Kepler solver's tolerance is 1e-12 (it was 1e-10), as the Python one.

### Added
- Tests (GEN.140, GEN.374): the Stumpff functions against exact rational arithmetic, the series-limit join, DOUBLE period and phase columns, radial (h = 0) motion, the wrapped turn against exact arithmetic, and Python and JavaScript agreeing on the Stumpff functions, the wrapped turn and the positions of every body at 3,000, 27,000 and 100,000 years. `docs/design/orbital-solvers-and-integrators.md` section 7 lists which guards are built and which wait on code not yet written.
