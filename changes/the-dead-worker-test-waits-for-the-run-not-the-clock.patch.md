**Headline:** The test that a dead worker never hangs a run no longer fails when the machine is merely slow.

### Fixed
- **`test_a_dead_worker_never_hangs_the_run` timed out once for two of its three attempts under heavy load (TEST.277).** It gave the whole run 60 s, and on a machine at load 14 to 20 starting the two worker processes, losing one and starting a fresh one took longer. It now waits for the run and gives up only when the run's record (the number of tasks in each state) stands still for 60 s, which is what a hung run does; a slow one keeps changing it.

TODO items resolved: TEST.277. Added: none.
