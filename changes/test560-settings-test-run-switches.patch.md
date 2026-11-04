**Headline:** The settings test knows the two test-run cache switches, so it passes on main again.

### Changed
- **Test fix, no code change (TEST.560).** `test_the_code_reads_no_variable_the_model_does_not_name` fails on any `PLANETGEN_*` variable the code reads and the settings model does not name. `PLANETGEN_BAND_TABLE_CACHE` (`generation/band_table_cache.py`, PERF.218) and `PLANETGEN_BATTERY_DIGEST_CACHE` (`generation/battery_digest_cache.py`, PERF.209) are read in code and set only by the test suite, so they are not `config.json` options; they are now listed with the other non-option variables, with their meaning (test-run switches).
generation-output: unchanged
