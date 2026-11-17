**Headline:** The test suite now runs a generation run's units and openings through the launcher, the way the product does by default.

### Changed
- Tests set `queue.forked_jobs` on by default (`PLANETGEN_FORKED_JOBS=1` in `tests/conftest.py`); the inline test queue hands a ledger run's ready rows to the launcher with an inline forker whenever the follower asks for workers, and leaves a test's own time limit running for the command, not the unit (OPS.228 slice D2b-2b). Tests that cover the RQ path set it off for themselves.
