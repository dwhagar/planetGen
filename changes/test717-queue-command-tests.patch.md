**Headline:** tests only: the tests that run a queue command in a process of their own no longer read the machine's data folder through a leftover worker.

### Changed

- TEST.717: `datasets status`, `check-math` and `plan` run through the queue, which hands the work to any worker already serving the Redis database, in that worker's environment. A worker left on the test database by an earlier test answered with the machine's `/var/cache/planetgen/datasets` instead of the test's folder, so `test_datasets_gate` and `test_math_gate` failed together on main. The new `queue_command_env` fixture ends such workers and empties the database before the command runs (and after), and four tests take it. A guard test (`tests/test_queue_command_env.py`) fails for any test that runs one of those commands in a process without it, and checks that a decoy worker is ended. No program code changes.
