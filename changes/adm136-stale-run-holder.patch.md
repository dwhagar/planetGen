**Headline:** A New galaxy job no longer waits for ever on a generation run that nothing follows any more: a run left behind by a process that is still running is freed, and a run whose set-up fails ends instead of holding the machine.

### Changed

- **A run nobody follows is freed, even when its process is still running (ADM.136).** The run line freed a holder only when its process id was gone, or, for a run on another host, when its `alive` key had run out. A run left behind by a process that lives on (a worker that stays up, or a process id the system reused) counted as live for ever, so a New galaxy job that came after it waited behind it until its time ran out. A run now counts as dead when its record is over, when its process id is gone, or when its `alive` key has run out, which a followed run never lets happen.
- **A followed run keeps its `alive` key while its follower is busy (ADM.136).** The follower said so only at each poll, so a callback that took longer than two minutes (saving what a stage found) would have let the key run out and the run look dead. The process that follows a run now says so from a thread of its own every 30 seconds until the run ends.
- **A run whose set-up fails ends as failed (ADM.136).** The set-up after the machine is taken (making sure workers run, the debug relay, the run's log) was outside the block that ends the run on an error, so a failure there left the run live for a process that went on; the next run waited for it.
- `src/tests/test_run_control.py` (found on the way): `test_a_started_unit_is_sent_the_stop_command` passed the fetched job object to `queue.remove`, which does not take it for a job when the queue's job class differs from the fetched one; it passes the id.

generation-output: unchanged
