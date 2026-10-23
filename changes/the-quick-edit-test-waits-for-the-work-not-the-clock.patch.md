**Headline:** The quick-edit test waits for its job instead of the clock, and two flaky browser and run tests now say what they saw.

### Fixed
- **`test_a_quick_edit_answers_in_the_same_response` failed once in a gate (TEST.289).** A route that queues quick work waits `SHORT_WAIT_SECONDS` (8 s) for the worker and answers `202` after that, so on a machine busy with a whole gate the test saw `202` and not `201`. The test now gives the route a long wait, and its message carries the status and body. `test_a_neighborhood_job_runs_its_sectors_through_the_wrapper`, which failed once with no text, now puts the answer, the sector count and the runs with their states in each assertion message.
- **`test_galaxy_breadcrumb_always_matches_the_view_its_url_names` timed out waiting for the steps of a fresh page (TEST.287).** The helper that opens the page again now reports its address, title, text and logged errors when no steps appear, so the next sighting names the cause.

TODO items resolved: none. Added: none.
