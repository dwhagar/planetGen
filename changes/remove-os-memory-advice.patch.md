**Headline:** The installer and updater no longer print advice about systemd memory settings or the MySQL redo log; the memory problem is solved in the program, not in the operating system.

### Changed

- `install.sh` and `update.sh` no longer run the redo log check (`scripts/redo_log_advice.py`, PERF.222) or the out-of-memory protection check (`scripts/oom_advice.py`, PERF.223), and both scripts, their tests and the INSTALL.md text that told you to add the `OOMScoreAdjust` drop-ins, the `MemoryHigh`/`MemoryMax` limit on `apache2.service` and the redo log capacity are removed (Boss, 2026-10-31, 6:07 PM and 6:09 PM EDT: no systemd drop-ins and no Apache memory limit, fix the memory use in code). Memory is held down by the worker caps (PERF.220, and PERF.226 for every stage). generation-output: unchanged
