**Headline:** A job process that ends on a plain SIGKILL is no longer recorded as a memory kill; only a rise of the kernel's memory-kill count while the job ran makes it one.

### Fixed
- OPS.294: the launcher recorded every SIGKILL of a job process as `killed_memory`, so a person's `kill -9` of a worker (or a test's) counted in the pressure report as a memory kill and fed the memory gate's self-tuning (the gate's edges move down after 3 kills) and a bigger retry claim. A SIGKILL is now `killed_memory` only when the kernel's own count of memory kills (`oom_kill` in the cgroup's `memory.events` and in `/proc/vmstat`, `memorygate.oom_kill_count`) rose between the job's start and its end. Any other SIGKILL, or one with no readable count, is recorded as `lost` with the text "ended on SIGKILL; the kernel's memory kill count did not rise while it ran, so it was not killed for memory". The job still goes again; only the memory record stops learning from it. Tests: `test_launcher.py` (a kill the count did not meet is `lost`, with and without a readable count).

generation-output: unchanged
