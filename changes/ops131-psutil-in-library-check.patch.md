**Headline:** The update's library check now installs `psutil`, so the work queue's workers start again after the memory gate (PERF.254).

### Changed
- **`psutil` is on update.sh's list of required libraries (`scripts/install-python-deps.sh`).** PERF.254 added `psutil` to `setup.py` and the lock files but not to the script's list, so an existing install never got it and every worker died at its start on `import psutil` (exit status 1), which stopped the math check at stage 4 of the update. `--check` now sees it missing and installs it from the lock.
