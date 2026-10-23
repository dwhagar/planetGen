**Headline:** Three four-worker draw tests start their workers fresh instead of forking a threaded test process.

### Fixed
- **`test_a_draw_is_deterministic_from_its_seed_at_one_and_four_workers` of the comet population hung until the 900 s timeout once under load (TEST.293).** The test, and the same test in `test_atmosphere_recipes.py` and `test_snow_line.py`, ran its four draws in a `ProcessPoolExecutor` that forks. A pytest-xdist worker has threads of its own, and Python warns that forking such a process "may lead to deadlocks in the child" (a lock another thread held at the fork is never released in the child). The pools now use the `spawn` start method, as the other multi-process tests and `sector_paths` already do. The stellar-composition test of the same kind still forks, because its children need the fixture datasets that test installs in the parent.

TODO items resolved: TEST.293. Added: none.
