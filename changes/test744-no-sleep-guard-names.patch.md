**Headline:** The guard that keeps sleeps out of command code names the functions that sleep today, so it is green again.

### Fixed
- `test_no_sleep_in_commands.py` listed `queue/wrapper.py` `_await_plan` and `_await_conclusion` and `generation/run_galaxy.py` `_take_to_the_deepest_floor`. The layered-fixes cleanup (#6454) moved the wait of the first two into the one `_await`, and moved the retry of the third into `contention.run_with_retries`, which the guard already declares as a blocking primitive. The list now holds `_await` (the follower of a run, waiting for a stage's plan or conclusion) and no longer holds the two names that stopped sleeping. The guard itself is unchanged: every sleep still has to be a listed wait. Part of TEST.744.

generation-output: unchanged
generation-affects: none
