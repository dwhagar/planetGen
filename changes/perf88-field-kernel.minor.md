**Headline:** Saving the sector paths is about twenty times faster: the integrator now runs as one compiled numba call.

### Changed
- The sector-path integrator (mass selection, integration, exit and spline simplification) is a numba kernel in `planetgen/physics/field_kernel.py`, 40 microseconds a path instead of 840. The saved paths are bit-identical to the old integrator's. The old Python loop is gone, with no fallback.
- `field_kernel.accelerations` and `field_kernel.potentials` give the Plummer-softened pull and potential of a `MassTable` at any points, the shared force every later gravity feature uses.
- `sector_path.integrate_path` takes a `Cell` (`sector_inside`, `box_cell`, `ball_cell`) for the region instead of a Python callable, because the compiled loop cannot call back into Python.
