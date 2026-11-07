**Headline:** A warm-up that starts while another one runs now waits for it instead of failing, and the phase lines give the time the phase took.

### Fixed
- **update.sh's warm-up part 2 failed with `run ... is running` when the epoch's warm-up had started the moment the update ended.** A warm-up that finds another maintenance run holding the machine now says which run it waits for and waits for it to end (up to four hours), then goes on and finds what the first stored; the waiting writes nothing to the phase record, so the running warm-up's phases on the admin page stay as they were.
- **The `Phase N of 5` lines added up the units' own seconds** (5338 s for a phase that took 9:39 with nine chunks at once). They now give the time the phase took; the admin page keeps the worker seconds as `worker_seconds`.
- **The phenomena line said `chunks of up to 8,543,000 phenomena`** because the first slice's share of the objects the cap leaves out was counted into the chunk size. A slice now reports its own size.
- **A work unit whose process was ended by SIGTERM was reported as a probable out-of-memory kill.** The failure now names the signal: SIGKILL is the kernel's out-of-memory killer, SIGTERM is something asking the process to end.
