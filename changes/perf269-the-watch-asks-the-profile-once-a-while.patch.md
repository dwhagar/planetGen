**Headline:** The supervisor's watch no longer opens a database connection on every beat to ask the job profile for a kind it has too few samples of.

### Changed
- PERF.269: `profiletime.unit_seconds` kept a figure for the run but not the answer 'too few samples', so the heartbeat's judge read the profile table through a new control connection on each beat (457 connections in a 130 s run; 10 of the 14 s they cost was building the SSL context). The answer 'nothing yet' is now kept for `NO_FIGURE_RECHECK_SECONDS` (30 s) and then asked again, so a kind still picks its figure up once the run's own units add the samples. Nothing a run computes changes.
