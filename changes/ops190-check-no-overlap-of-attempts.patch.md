**Headline:** Test that a retry never overlaps a live attempt, closing the old OPS.188 and OPS.190 audit bugs against the new engine.

### Added

- A real-process test: a job whose process is alive but silent (stopped) is not handed out again by `recover` or `turn`, keeps its first attempt, and finishes once with every piece done once when it goes on. With the fence tests of `test_attempts.py` (a replaced attempt's late commit rolls back) this shows OPS.190 fixed; OPS.188 (the unit watch, the top-up and the run's alive key outside the starter) is shown fixed by the supervisor carrying a run whose starter was killed (`test_supervisor_kill.py`).
