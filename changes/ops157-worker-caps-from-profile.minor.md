**Headline:** The number of workers a stage may run follows what that kind of unit has really used on this machine, and a stage starts with the worker limit where the last run of it settled.

### Changed
- **OPS.157.** `workers.stage_worker_cap` and `stage_worker_mb` take the connection and the run, and size a stage's workers by the claim the gate makes from the job profile (OPS.156; the typed table with no history), so the cap and the claim cannot disagree and a run holds one figure. The adaptive worker control remembers, per kind of stage, the limit it settled at after a memory wait (`adaptive._remember`, forgotten when memory has been free up to the cap), and a later stage of that kind starts there (`adaptive.start_limit`) instead of at the configured cap; a limit from another machine shape counts half. Both are logged with the kind, old and new value and the reason, and kept for the Queue page. `PLANETGEN_SELF_TUNING=0` restores the typed caps and the configured starting limit.
No drawn value changes (the work seed does not depend on the worker count), so there is no epoch or golden.
generation-output: unchanged
