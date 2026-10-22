**Headline:** The generate tests no longer depend on the free space of the machine that runs them, and dead copies in two modules are gone.

### Fixed
- **Five tests failed on a machine whose drive was nearly full (TEST.268).** A generate run refuses to start when it would leave under 5 GB free (PERF.3) and `--strict` turns that into exit status 1, so `test_an_address_named_outside_the_outline_is_generated_from_its_profile` (3 cases), `test_neighborhood_around_a_sector_outside_the_outline_takes_the_reach_and_says_how_big_it_is` and `test_estimate_only_writes_nothing` failed wherever the free space was under 5 GB (a test MySQL server's binary logs had filled one machine to 97%). The test suite now shows the code a drive with at least 200 GB free, unless a test sets its own.
- **Dead code left by earlier moves is gone.** `web/lib/xfetch.py` still held copies of `is_cacheable` and `shared_key` that named a constant and a module it never imported (the live ones are in `pagecache.py`), and `cli/orbits.py` defined `_stage_label` twice.

TODO items resolved: TEST.268. Added: none.
