**Headline:** The test that interrupts a parallel galaxy run no longer fails as "Busy" when another test worker has a run going.

### Fixed
- **`test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector` failed under load with "run ... is running" (TEST.266).** The test starts `planetgen generate` as a process, which takes the machine's one live-run lock (`planetgen:run:active`) on the shared Redis server; a process cannot be given the key prefix the in-process tests get, so a second test worker running at the same time was refused. Reproduced by starting the test twice at once (before: one of the two failed Busy; after: both pass). Each worker now takes the Redis database numbered after it for the processes it starts (`tests/private_redis.py`).

TODO items resolved: TEST.266. Added: none.
