**Headline:** The installer and updater now check that MySQL and Redis are protected from the kernel's out-of-memory killer, and offer a memory limit for the generation workers.

### Changed

- `scripts/oom_advice.py` (run by `install.sh` and `update.sh` after the redo log check, advice only) prints the `OOMScoreAdjust=-900` drop-in for the MySQL and Redis units that are loaded and not at -900 or lower (INSTALL.md said -800), and an optional `MemoryHigh`/`MemoryMax` drop-in for `apache2.service`, where the workers run, sized at 55 and 65 percent of the host's memory. INSTALL.md says so. Nothing on the server is changed. generation-output: unchanged
