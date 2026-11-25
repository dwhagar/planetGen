**Headline:** The interrupted galaxy run test builds a sector of 4 systems instead of 16, so the first sector is saved well inside the test's 480 s wait.

### Fixed
- TEST.822 (b): tests/test_interrupted_galaxy_run.py failed on main and on a lane box with "no sector was saved" because a sector of 16 systems took longer than START_WAIT_SECONDS (480 s) on a box running the whole gate. The test is about interrupting a run, not about generation size; a sector of 4 systems is saved in about a quarter of the time and the run still has hours of sectors left when the signal lands (shell 5 plans 296 sectors). The 480 s limit is unchanged.

generation-output: unchanged
