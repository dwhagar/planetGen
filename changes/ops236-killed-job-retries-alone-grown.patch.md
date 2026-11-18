**Headline:** a job the kernel kills runs again alone on a claim one and a half times as large, as the design says (a bug since RQ was removed).

### Fixed
- OPS.236 (bug): when OPS.235 deleted `units.horse_killed`, nothing called `memorygate.grow` any more, so a job process killed by a signal (the kernel's out-of-memory kill) was handed out again with the same claim and beside other jobs, and could be killed again. `Launcher._grow` now records the grown claim of a retry after a kill, and `Launcher.claim` gives that retry the grown figure and the `alone` token (read from the earlier attempt's error kind, no new column). A unit that failed on its own, or was lost with its supervisor, is not grown.

generation-output: unchanged
generation-affects: none
