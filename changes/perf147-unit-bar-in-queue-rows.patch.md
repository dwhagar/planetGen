**Headline:** The Running units table of the Queue run page shows each unit's own progress bar and time left under its counter (PERF.147 in the OPS.139 rows).

### Changed
- **The unit rows draw the bar macro.** A unit that is expected to take over 15 seconds, or has run that long, shows its bar and a time left ("30%, about 120 s left", with "(from earlier runs)" while the time left is still the recorded rate's) under its counter; a unit that ends sooner shows none. Reused: `partials/task_bar.html` and the `progress` of each unit in `GET /api/admin/work/<id>/units` (`taskprogress.as_dict`), the unit rows of OPS.139 (`queue_page._unit_rows`).

generation-output: unchanged
generation-affects: none
