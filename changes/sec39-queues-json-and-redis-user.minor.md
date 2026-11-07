**Headline:** No queue carries pickle any more, and the example Redis configuration ships a password-protected user that update.sh asks for.

### Changed

- `queue/redisqueue.py` (SEC.39, audit finding S-M3): `serializer_for` returns the JSON serializer for every queue name, `JSON_QUEUE_PREFIXES` is gone, `fetch_job` and `cancel_now` default to JSON, and a worker of any queues uses it (the check that one worker serves only one kind of queue goes with it). A run's units already went through the shared JSON queues since PERF.95 to PERF.98; the pickle path that remained served a queue nothing writes to.
- `examples/redis/planetgen-acl.conf` ships the Redis user `planetgen` (default user off, no dangerous commands, INFO and CONFIG GET allowed) with the steps to install it and put its password in `redis.url`; `examples/redis/planetgen.conf` points to it. The queue, worker, update-guard and stop-all tests pass as that user.
- `scripts/check_redis_settings.py` (run by `update.sh`) warns, with the lines to add, while the connection is the passwordless default user. A job still names the function a worker calls, so the password, not the serializer, keeps another local user from queuing one.
- Tests: every queue name serializes as JSON, nothing that builds an RQ queue, worker or job fetch leaves the serializer to RQ, the queue package never imports pickle, and the warning and the example file agree.

generation-output: unchanged
