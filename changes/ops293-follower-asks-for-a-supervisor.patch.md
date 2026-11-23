**Headline:** A command that follows a run (check-db and every other command on the work queue) now asks for a supervisor while it waits and stops with the reason after a minute when none comes up, instead of waiting for ever with no word.

### Fixed
- OPS.293: the supervisor exits when idle and the follower never started one, so `check-db` against a finished run (and a run after its supervisor was killed, drill D1 of the OPS.247 trial) waited with no message, 15 minutes in the trial, until someone started a supervisor by hand. Every follower (`wrapper.run`, `wrapper.follow`) now calls `wrapper.supervisor_pump` every `units.ENSURE_SECONDS`: it asks for a supervisor (one Redis read while one is alive), and when none has been alive for `wrapper.SUPERVISOR_WAIT_SECONDS` (60 s) it stops with `units.RunStopped` saying so, why the start failed (`supervisor.last_failure`, the error `ensure_started` swallowed before) and how to start one by hand. Tests: `test_supervisor_pump.py`.

generation-output: unchanged
