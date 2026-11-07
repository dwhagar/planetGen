**Headline:** The two queue time-limit tests that were red on main pass again: `requires.enqueue` names the job limit as a keyword, and the over-limit test builds its worker the way the real one is built.

### Changed

- **`requires.enqueue` takes `job_timeout` as a required keyword (TEST.667).** The guard `test_every_call_that_queues_a_job_names_its_limit` reads the source for a `job_timeout` at each queuing call and could not see one behind `**options` in `queue/requires.py`. The limit is now a named keyword and is passed on by name; a call without one fails with a plain `TypeError` before it reaches the queue (which refuses a job without a limit already). No caller in the program changes: they all name it.
- **The over-limit test reads the queue's JSON jobs (TEST.667).** `test_a_job_over_its_limit_is_ended_and_reported_failed` ran a worker with RQ's default pickle serializer and plain job class against a queue that holds JSON jobs, so the worker could not read the job (`DeserializationError`) and never reached the limit. The worker now has the JSON serializer and the guarded job class, as `cli/worker.py` builds one, and `time.sleep` is added to its list of functions the way the suite does.

generation-output: unchanged
