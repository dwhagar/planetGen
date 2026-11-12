**Headline:** the generated-system belt dose test follows today's generated systems: a body with no stored dose has no belt dose, and the giant it needs is found by search.

### Changed

- `test_a_generated_system_stores_the_belt_dose_in_its_surface_dose` failed on main for two reasons, both in the test: (1) it expected a belt dose of 0.0 on every planet that is not a giant or a belt, but a class T body can store no dose at all (`surface_dose_msv_yr` and `dose_belt_msv_yr` both `None`, like a giant); it now expects `None` for a body with no stored dose and 0.0 for one with a dose. (2) Seeds 5 to 59 now hold one giant in all and it has no moon, so its fallback ("put a moon of the first giant in the belts") found none; it now searches systems forced to hold a class I planet for a giant with a moon (`tests/fuzz_support.first_seed_with`), and gives that giant a field, since a generated giant often has none. Nothing stored changes.
