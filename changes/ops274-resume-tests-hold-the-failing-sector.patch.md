**Headline:** the resume tests now say what a failing sector does (the run holds it, finishes the rest and ends held), instead of expecting the run to stop at it.

### Changed
- OPS.274: `tests/test_gen_resume.py` (six cases, three modes each of a sector that fails before its save and a save that fails in its transaction) asserts what OPS.221 made true: the run raises `RunStopped` with "The run is held", the interrupted run holds every sector of the uninterrupted run but the one that failed, that sector has no row and no system, and a new run ends with the same sectors as the uninterrupted one. The uninterrupted run is made first, so the interrupted run's sectors are compared with it.

generation-output: unchanged
generation-affects: none
