**Headline:** The launcher that hands out a forked run's units now claims for each unit what the memory gate has learned of its command, and holds each stage to its share of the host's memory, instead of claiming the 150 MB minimum.

### Fixed
- A unit row had no memory figure, so the launcher claimed `MIN_CLAIM_MB` for every unit and the gate would have let as many 3 GB units start as there were cores. A unit now claims the typed figure raised or lowered to the job profile's (as the RQ path's `claims.claim_for` did), its cores and named resources (`unit_tokens`), and an opening claims the cut's whole lane (OPS.228 slice D2b-2a).
- A stage runs at most `workers.stage_worker_cap` units at once (PERF.220's share of the host's memory in units of the stage's figure), besides its run's `--workers` count; the supervisor gives its launcher its Redis to read the learned figures from.
