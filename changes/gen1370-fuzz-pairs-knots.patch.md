**Headline:** The tuning fuzz check no longer mistakes the white dwarf core table for low-high ranges.

### Fixed
- **`test_every_numeric_constant_is_finite_and_every_pair_ordered` failed on `WD_CENTRAL_OXYGEN_BY_MASS_SOL` (GEN.1370).** The check reads every two-number tuple as (low, high); the table's knots are (white dwarf mass, central oxygen fraction), so the check now lists it with the other knot tables. Test only. generation-output: unchanged
