**Headline:** `scripts/test.py` ends the queue workers an earlier run left behind before it starts, not only after it.

### Fixed
- **A queue worker left by a stopped test run took the units of the next run.** The runner ended a run's leftover workers when the run finished, so a run that was killed (a spend-limit stop, a lane restart, a phase stopped from outside) left them running, some on old code, and the next run's units went to them. The runner now ends every test worker (those that name a test Redis database) before its first phase and says how many it ended (TEST.398).

TODO items resolved: TEST.398.
