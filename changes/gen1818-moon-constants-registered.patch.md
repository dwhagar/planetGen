**Headline:** The constants `moon_orbits` gained in GEN.1816 and GEN.1818 are registered stand-ins (TM10, TM11), so `test_constants_completeness` stays green on Train 52.

### Fixed

- `moon_orbits.MAX_LOG_DISTANCE`, `FAR_DISTANCE_M` (numerical guards, TM10, T9, L0) and `HEAT_ATOL_SHARE` (TM11, T6, L1) have rows in section 6.25 and a tag in their docstrings; the baseline is untouched.

generation-output: unchanged
