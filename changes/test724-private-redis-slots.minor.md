**Headline:** Two test runs can share one machine: `--slot N` gives a run its own three Redis databases, so they no longer empty each other's queues.

### Added
- **`--slot N` on `scripts/test.py` and `scripts/train_gate.py run`** (TEST.724): five slots of three Redis databases each (`tests/private_redis.py`, `PLANETGEN_TEST_REDIS_SLOT`). A slot run takes only its slot's databases, has at most three workers and ends only its own stray queue workers (`stray_workers.reap_all(databases)`), so two gates on one machine never meet in Redis (the `RunStopped: the run None` failures). A run with no slot is unchanged. MySQL needed nothing (every test database is named at random), and the gate already reads `scripts/data/known_reds.json`. Documented in `docs/reference/testing.md`.

generation-output: unchanged
generation-affects: none
