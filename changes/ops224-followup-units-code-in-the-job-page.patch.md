**Headline:** Three places that wrote a unit by hand (the job page's memory claim, a retry hold's seconds, the unit beat's debug line) now use the shared units formatters.

### Fixed

- OPS.224 follow-up: `queue/jobview.py` (the retry hold's time left), `web/templates/admin_queue_job.html` (a job's memory claim) and `queue/unitbeat.py` (a restarted unit's wait and its debug beat line, OPS.223) wrote a number followed by its unit; they use `format_duration_seconds`, `format_megabytes` and the page's `bytes_text` (UX.318), which `test_hand_formatted_units` guards.

generation-output: unchanged
generation-affects: none
