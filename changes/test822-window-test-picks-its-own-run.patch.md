**Headline:** The window test that carries a run on picks the run it started, and the fold of a finished stage keeps the same slowest units each time.

### Fixed
- TEST.822 (a): `test_window_spawn.py::test_a_run_carried_on_from_one_that_stopped_in_a_later_window...` took the newest run id of the list as "its" run. A run of an earlier test that the run store restores from the control database (OPS.284) after the test's snapshot joins the list in the same second, and its id sorts by its random end, so in about half the runs the test finished the wrong run ("did not record the galaxy and the build"). It takes the new run that recorded its galaxy seed now.
- `planrows.fold_finished` kept "the `FOLD_KEEP` slowest" units in an order the database chose among equal times; folding a stage twice could keep different rows. It breaks ties by the lower job id, so a second fold finds nothing left to fold (the `folds_the_history` case failed on main in about half the runs).

generation-output: unchanged
