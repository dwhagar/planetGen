**Headline:** A trojan region or planet is no longer taken for a belt or planet out of order or too close to its giant, which failed the wide binary check.

### Fixed
- GEN.1150, GEN.1151: a trojan region is put at the end of its star's list and stands at its giant's distance, and a trojan planet stands there too, so `check_orbit_spacing` called a region "listed after the belt at 24 AU but orbits inside it" and "too close to the belt" (`test_wide_binaries_never_overlap`), and the final respacing pass could push them off the giant's orbit. `validation.slot_planets` now leaves out any body with a `lagrange_point`, as it does a promoted belt member: they share the giant's orbit by design, in step with it, so they are no slot of the spacing rules. The check itself is not loosened for any other body (a test pins that a belt listed out of order is still caught). Where the final respacing used to move a trojan, it no longer does; main has no trojan, so only the Train 52 head moves.

generation-output: changed
generation-affects: fill
