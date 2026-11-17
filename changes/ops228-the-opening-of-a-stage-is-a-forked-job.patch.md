**Headline:** With forked jobs on, the opening of each later stage of a run (its cut into units) is handed out by the supervisor's launcher too, as a job of its own, instead of an RQ job.

### Changed
- The opening of a stage is a ledger row (`opening`) that claims the cut's memory; the launcher forks it like a unit, writes its attempt (with the process's id) and marks it done, canceled when the run was no longer live, or retried once when its process died (OPS.228 slice D2a).
