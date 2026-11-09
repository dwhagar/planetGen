**Headline:** A queue test can no longer replace its own test process with a worker, which took a parallel test node down.

### Changed

- `tests/test_work_units.py`: `test_a_unit_survives_a_warm_that_outlasts_an_execution_record` ran the real queue worker in the test process. The worker starts itself over (`os.execv`, PERF.226) for a unit that needs no star tables when its process holds them, and a pytest-xdist node that an earlier test had left holding the tables was replaced by a worker ("node down" on every parallel run, passing alone). The test now says no tables are held.
- `tests/conftest.py`: an autouse fixture makes `os.execv` raise in a test process unless the test patches it (the tests of the worker's start-over already do), so the next such case fails with a message and not with a lost node.

generation-affects: none
