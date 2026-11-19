**Headline:** The life-constant registry test knows `habitability.GAS_ORDER` is the names of four gases, not a value.

### Fixed
- `test_every_life_constant_has_a_row_or_a_reason_to_have_none` failed on `habitability.GAS_ORDER` (CO2, CO, H2S, SO2; GEN.1772): it is the order the four gas bands are read in, a structure like `DOMAINS` and `TIERS`, so it gets a reason to have no stand-in row in `NEEDS_NO_ROW`. The test is not loosened and the generator is not touched.

generation-output: unchanged
