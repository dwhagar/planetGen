**Headline:** A work-unit test that waited for ever when Redis or MySQL was missing now skips like the tests beside it.

### Changed

- **`test_a_unit_runs_with_the_switch_its_run_was_started_with` asks for the Redis fixture (TEST.682).** It started a run through the queue wrapper but, unlike every test beside it, did not take the `conn` fixture, so with no test Redis configured the inline queue was not in place, no worker took the units and the test waited in `follow_stage` for ever (both parameters; the hang guard fails it after the time limit). A run of `test_cache_warmup.py`, `test_warm_phases.py`, `test_warm_progress.py` and `test_work_units.py` with Redis only reproduced it. It takes `conn` now: it skips without Redis and passes with it. The test named in the report, `test_a_tiered_pass_travels_as_json_and_comes_back_equal`, passes in every combination. Test only.

generation-output: unchanged
