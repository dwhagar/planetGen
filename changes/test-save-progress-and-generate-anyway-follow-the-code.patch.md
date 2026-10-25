**Headline:** Two tests that were red on main (the sector save report that took over 15 minutes, and the Generate anyway log) now test what they meant to.

### Fixed
- **`test_save_progress.py::test_a_sector_save_reports_each_system_it_inserts` timed out at 900 s.** It stored a galaxy shape first, so the sector was filled from the galaxy's density and the `num_systems = 3` it asked for no longer limited it (817 to 1,532 systems and 7,000 to 13,000 phenomena at three addresses, over 15 minutes at the address it used). The test is about the save's progress report, which is the same for any sector, so it no longer stores a shape: the sector has 3 systems and the test takes 2 seconds.
- **`test_web_generate.py::test_generate_anyway_starts_the_job_and_logs_it` and `test_a_plain_confirm_is_not_logged_as_an_override` counted every activity-log event.** The page's in-process API calls are logged as `call` events by the same function, so the list was `["call", "call", "job.start", ...]`. The two tests now keep only the `job.` events they assert on.
- Not changed: `test_db_column_round_trip::test_every_column_is_written` passes alone (a load failure in a long run), and the `test_api_auth_sweep` failure was fixed by the run log download change (SEC.36's PR).
TODO items resolved: none (TODO closes TEST.375 after reading this). Added: none.
