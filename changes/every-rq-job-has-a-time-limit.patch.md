**Headline:** Every job queued to RQ now has a time limit taken from how long the system expects it to take, so a job whose worker was killed expires instead of staying on the queue for ever.

### Fixed
- **The Generate page's jobs, the API's jobs and the work queue's tasks were queued with no time limit (PERF.144).** Each one had `job_timeout=-1`, so RQ never expired the entry of a job whose worker died. Now every site gives a finite limit, ten times the best estimate, never under 15 minutes and never over its ceiling: a Generate page job from the recorded time of its steps (a floor of one hour, a ceiling of 30 days, which is also the limit of a job with a step that was never recorded), an API job from its route's estimate (the neighbourhood's `estimate.seconds`, a regeneration's recorded sector time, a slow page's expected seconds; seven days at most, and one hour for an edit nothing can estimate), a work queue task from its submitter's estimate or the recorded seconds of its kind (6 hours at most and when nothing is recorded), and the tile prefetch and rebuild from the size of the chunk. A units' limit is as before.
- **A Generate page job that RQ ends for its limit now stops its step.** The runner stops the step it was waiting on and reports the job failed with the reason, where the step used to run on.

### Changed
- **`redisqueue.queue` returns a `TimedQueue`.** It refuses a job queued with no limit, a limit of -1, zero or `None`, so a new kind of job cannot skip `redisqueue.job_timeout`; a test also fails on any `enqueue` call that does not name its limit and on any module that makes an RQ queue or job of its own. `units.UNIT_TIMEOUT_SECONDS`, `units.TIMEOUT_MULTIPLE` and `units.MIN_UNIT_TIMEOUT_SECONDS` are now `redisqueue.MAX_TIMEOUT_SECONDS`, `TIMEOUT_MULTIPLE` and `MIN_TIMEOUT_SECONDS`.

TODO items resolved / Added: PERF.144 / Added: none
