**Headline:** Two warm-ups never run at the same time: the update's warm-up has priority and cancels any other, a change's warm-up waits, one started by hand refuses; the phase lines give the time the phase took.

### Fixed
- **update.sh's warm-up part 2 failed with `run ... is running` when the epoch's warm-up had started the moment the update ended (OPS.117).** The update's own warm-up now has priority: it cancels any other cache warm-up that holds the machine (the epoch's, the one after an orbital update, one started by hand) through the normal cancel path, waits for the canceled run's units to stop, and goes on. The canceled warm-up ends with one plain line, not as a failure.
- **A warm-up that is not the update's never starts while another runs.** The one a change starts waits for the running one to end (nothing written to the phase record meanwhile) and then warms what is left; one started by hand says which run holds the machine and exits without starting.
- **A canceled warm-up could write into the phase record of the run after it.** The record belongs to the run that began it; late writes from any other process are dropped, so the admin page never shows two runs at once.
- **The `Phase N of 5` lines added up the units' own seconds** (5338 s for a phase that took 9:39 with nine chunks at once). They now give the time the phase took; the admin page keeps the worker seconds as `worker_seconds`.
- **The phenomena line said `chunks of up to 8,543,000 phenomena`** because the first slice's share of the objects the cap leaves out was counted into the chunk size. A slice now reports its own size.
- **A work unit whose process was ended by SIGTERM was reported as a probable out-of-memory kill.** The failure now names the signal: SIGKILL is the kernel's out-of-memory killer, SIGTERM is something asking the process to end.

generation-output: unchanged
