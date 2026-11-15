**Headline:** Every attempt at a queue unit now leaves a full row in the job ledger (what it cost, how it ended, its stack), and a unit that runs again writes a debug log of its own.

### Added

- OPS.223 (job engine stage 3, L9; absorbs OPS.146): an attempt's row in `engine_attempts` keeps its processor seconds and peak memory, and for a failure its error kind, its traceback and a stack file in the log folder's `unit-debug` folder. An attempt nobody could close (its horse killed for memory or `kill -9`, a unit stopped as stuck, a process gone) is closed by whoever saw it go as `killed_memory`, `lost` or `stuck`, with the stack the follower had the horse write; the next attempt closes any left open as `lost`; the time limit is `timed_out`.
- A unit that runs again (the second and later attempt: a retry after a kill, a stuck stop or a failure, a resume) runs in debug mode: every debug area, plus the beat's counter, memory and processor lines and each wait, goes to `job<job>-unit<unit>-attempt<n>.log` in the `unit-debug` folder, which opens with how the attempt before it ended and its stack; the run's log names the file. The new setting `queue.debug_restarts` (on by default) turns it off; each file is capped at 20 MB and the folder keeps the newest 400 files. `log.captured` is the one place a block of code is written to a file of its own at debug level.

generation-output: unchanged
generation-affects: none
