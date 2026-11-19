**Headline:** A generation run now shows as running with its start time in the job ledger, and the bright-star backfill stage writes a status line with real counts and a time left every five minutes.

### Changed

- OPS.281: a run's row in `engine_runs` stays `waiting` with no start time while its units run; it is now `running` with `started_at` set when its first job starts (`ledger.begin_attempt`), a held or ended run keeps its state, and a run that finishes done is `done` with `finished_at` in the ledger too (it followed a cancel and a hold before, not a finish).
- OPS.282: the backfill stage (the mass scatter from the neighborhood plus the sky tiers) took hours with one log line. `BackfillFollower` now writes a line every five minutes: sectors visited of those to draw, sky parts done of all, stars placed by the finished parts, what the working parts count inside themselves, and the time left as hours, minutes and seconds.
- A test used the old name `rq_units` after the module became `unit_helpers`; it uses the new one.

generation-output: unchanged
