**Headline:** the supervisor kill test that has been red since OPS.251 passes again: its stand-in for starting a supervisor now ends the claim to be starting one when the process ends, as the real start does.

### Fixed
- OPS.251 (bug, test only): `tests/test_supervisor_kill.py` started its supervisors with its own `start` function, which never ended the claim to be starting one (`supervisor._release_claim`) when the supervisor was killed; since OPS.251 slice 2 that claim has no time limit, so the next caller found it held by a live process (the test's) and started nothing. The product is unchanged: its own start (`supervisor._spawn`) ends the claim through the exit notice.

generation-output: unchanged
generation-affects: none
