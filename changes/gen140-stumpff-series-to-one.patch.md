**Headline:** The orbit stepper keeps all its digits for near-parabolic and slow-bending orbits.

### Changed
- **Stumpff functions by series for `|z| < 1`** (GEN.140): `kepler.stumpff_c2_c3` and the compiled `state_kernel._stumpff` sum ten terms of the series below `|z| = 1` (it was `1e-3`), where the closed forms cancel and a hyperbolic `c3` was off by 4e-13. Body positions the universal step gives on the maps change in the thirteenth digit and below; nothing stored by a generate changes (the step runs on read, `physics/body_positions.py`).

### Added
- Tests of the other guards in `docs/design/orbital-solvers-and-integrators.md` section 7: the Stumpff functions against exact arithmetic, DOUBLE period and phase columns, and radial (h = 0) motion. The section lists which guards are built and which wait on code not yet written.
