**Headline:** Four tests that check the Solar System stand-in rows pass again.

### Fixed
- TEST.648: the stand-in registry row SSR9 still named the ground type `oxidised` after DOC.103 respelled it `oxidized`, so `test_solar_system_standin_rows` could not match it to the registered row. The row's key is respelled, and the generated stand-in page is rebuilt.
- TEST.648: `tests/data/solar_system_rows.json` still listed the ozone rows `m1` and `m2` of table m, which GEN.524 retired, so `test_calibration_fence` failed on rows that no longer exist. Both entries are removed.
