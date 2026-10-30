**Headline:** The inner-system share anchors of the planet architecture are two ordered tuples, so the constants test passes again.

### Fixed
- **`test_fuzz_utils_and_config` failed on `INNER_SYSTEM_FRACTION_ANCHORS`.** The `(mass, share)` pairs were not in ascending order; the masses and the shares are now separate tuples with the same values. No generated output changes.
