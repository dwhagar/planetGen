**Headline:** The waiting-in-line tests hold the machine with a run another process follows, so the guard against a run waiting for itself stays as it is.

### Changed

- `test_work_units.py`: the three tests that make a second run wait for the first (`..._waits_in_line_...`, `..._do_not_depend_on_the_place_it_started_from`, `..._gives_up_after_the_wait_bound_...`) started the first run in the test's own process, which the ADM.136 guard (`nested_guard` in `runline.start_when_free`) refuses at once with `Busy`, as it should: a run started inside a run of the same process would wait for itself for ever. The first run is now held by a run that a live second process follows (`other_process`, `_held_by`), which is what a run from another terminal or job looks like. The guard is unchanged.
- New test: a run started inside a run of the same process is refused at once, naming that it would wait for itself, even when the wait bound is thirty minutes.

generation-output: unchanged
