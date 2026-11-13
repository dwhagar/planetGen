**Headline:** The air solver no longer keeps an earlier step's carbon sink, and an air that settles on its last step still gets the escape pass.

### Fixed
- **`surface_air.solve` (GEN.1594).** The carbon sink (`Air.co2_sink`) was written only on steps that weather, so an air whose last step did not weather kept an earlier step's sink for the life stage to balance against; the record is now cleared at the start of each step. And when the no-escape fixed point first settled on the last allowed step (`MAX_STEPS`), the escape pass never ran and the air was called converged; the loop now allows `ESCAPE_MIN_STEPS` more steps from the step that sets the escape.

generation-output: changed
generation-affects: fill
