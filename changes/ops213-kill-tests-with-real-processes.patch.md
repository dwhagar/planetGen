**Headline:** The supervisor is tested with real processes and real kills.
 
### Changed
- **OPS.213, output unchanged.** `tests/test_supervisor_kill.py` starts a real supervisor process on the suite's Redis and kills things with SIGKILL: the supervisor while a run waits in the line (the next ask starts a new supervisor, which starts the run), a second supervisor (it refuses the host's lock and leaves the first), and the process that started a run (the supervisor says once that it carries the run on, and the run stays live and is not taken over). `tests/test_supervisor_galaxy_kills.py` kills the unit that is ending a stage (the holder of the stage's end lease) in a real galaxy run and compares the galaxy with an undisturbed run's; it is skipped unless `PLANETGEN_ENGINE_KILL_TESTS=1`, as a small box cannot finish a galaxy.
generation-output: unchanged
generation-affects: none
