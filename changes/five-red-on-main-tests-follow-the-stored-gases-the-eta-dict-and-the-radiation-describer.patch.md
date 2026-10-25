**Headline:** Five tests that went red on main when the air moved into the gas tables, the queue ETA became a dict and the Radiation describer came back now match the code.

### Fixed
- **`test_age_needs.py` (2 tests): the test bodies kept their gases as `p_<gas>_kpa` keys.** Since DB.71 and GEN.424 (#3110) the rules read the air from `gas_partials_pa`, so a carbon-dioxide or reducing atmosphere counted as no air and came out `watered`. The test's `body()` now folds the keys into that field with `tests/air_support.fold_air`.
- **`test_object_rows.py` read-back test: its fake connection had no `planet_gases` or `moon_gases` rows,** so `body_gases` failed with `KeyError: 'partial_pa'`. The fake now answers those tables with each body's own gases, and the read-back compares them too.
- **`test_eta_estimator.py` queue page test passed a float where `_eta_text` takes the ETA dict** (`{"seconds", "stalled"}` since PERF.99). It now passes the dict, and a stalled and a seconds-less ETA are covered.
- **`test_protection.py::test_there_is_one_describer_of_the_six_categories` asserted `physics/habitability_explain.py` did not exist.** GEN.445 (#3052) brought that module back to name the composition behind the Radiation term only. The test now holds the rule it stood for: the module's public names are just the radiation ones, so a second describer of the six categories cannot appear in it.
- Not changed: `test_web_pages.py::test_in_process_transport_skips_http_inside_a_request` fails here with `DatabaseLoginError` for the `planetgen` account, the OPS.62 gap on lane machines, not a code bug.
TODO items resolved: none (these six tests were red on main with no item of their own). Added: none.
