**Headline:** The tests of the worker count no longer disturb each other through the one Redis key that holds the admin's wanted number of workers.

### Fixed
- Tests: `planetgen:workers:wanted` is one key on the Redis server, read by every `ensure_workers` and set by the admin's page and `set_wanted`. Tests that set it (`test_workers.py`, `test_work_queue_admin.py`) ran in other test processes beside tests that count the workers a queue should have, so a count of 2 came out as 3 (or the key was deleted mid-test) and `test_a_batch_of_jobs_starts_no_more_workers_than_wanted_and_not_twice` and `test_a_worker_that_is_killed_is_replaced_and_the_queued_job_runs` failed in loaded gates (TEST.255). Each test now has a key of its own (`conftest.py`), emptied after; `test_wanted_workers_key.py` pins it (it fails with the shared key).

TODO items resolved: TEST.255 / Added: none
