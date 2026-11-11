**Headline:** The squeezed-heliopause system-sheet test names its nebula, so the text it checks matches what the sheet says for a named cloud.

### Fixed
- **`test_a_system_inside_a_dense_nebula_reports_its_squeezed_heliopause` named an unnamed nebula.** `_placed_nebula` inserts a nebula with an empty name, and the sheet then says "the gas of the cloud around the system presses it in" while the test looked for "the gas of  around ..." with the empty name. The test now gives the nebula the name "Test Cloud" with the same update that sets its density and temperature, so the sentence it checks carries a name. A test-side fix only; no product code changed. Reused: the existing `_placed_nebula` helper and the sheet's own `or 'the cloud'` fallback.

generation-output: unchanged
generation-affects: none
