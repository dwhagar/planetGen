**Headline:** The supervisor is tested with real processes and real kills.
 
### Changed
- **OPS.213, output unchanged.** `tests/test_supervisor_kill.py` starts a real supervisor process on the suite's Redis and kills things with SIGKILL: the supervisor while a run waits in the line (the next ask starts a new supervisor, which starts the run), a second supervisor (it refuses the host's lock and leaves the first), and the process that started a run (the supervisor says once that it carries the run on, and the run stays live and is not taken over). The kill of a unit's process in the conclusion and the comparison with an undisturbed galaxy need a galaxy run and are left for the gate machine.
generation-output: unchanged
generation-affects: none
