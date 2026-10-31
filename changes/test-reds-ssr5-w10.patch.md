**Headline:** Two tests that were red on main are green again, with no change to what the generator makes.

### Fixed
- `test_solar_system_standin_rows` (SSR5): a registered row may show what replaces it in plain words (TEST.487 took the item code off the page), so the test now checks the code only where a row still shows one.
- `test_oddity_plausibility` (W10): the harness has an effect for the mantle-stripping row, so every shipped oddity row can be forced and held to the invariants.

generation-output: unchanged
