**Headline:** a test that expected a belt dose of 0.0 on every planet now accepts the bodies that store no dose at all.

### Changed

- `test_a_generated_system_stores_the_belt_dose_in_its_surface_dose` expected `dose_belt_msv_yr == 0.0` on every body that is not a giant or a belt, but a class T body under its envelope stores no dose at all (like a giant: `surface_dose_msv_yr` and `dose_belt_msv_yr` are both `None`, `physics/radiation.update_ground_dose`). The test now expects `None` for a body with no stored dose and 0.0 for one with a dose. Nothing stored changes.
