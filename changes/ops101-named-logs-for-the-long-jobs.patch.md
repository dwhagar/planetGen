**Headline:** The nightly whole-map build, the nightly orbit update and update.sh each leave a log of their own in the log folder.

### Changed

- `planetgen tiles` (the nightly whole-map build) appends each run to `tiles.log` and `planetgen orbits` (the nightly orbit update) to `orbits.log`, in the log folder next to `warmup.log` (OPS.101). A run is a start line that says what started it (`by hand` or `systemd timer`), the lines it printed with their time, zone, severity and component, and an end line (`finished after 17:44`, `ended with exit status 75`, `ended by Ctrl+C`, `failed: ...`). The shared writer is `util/joblog.py` (`JobLog`); `warmup.log`'s `RunLog` is built on it; `logpaths.job_log_path` names the four logs.
- `update.sh` appends every line it prints to `update.log`: the shell keeps the records from the start in a temporary file and writes them to the log folder once the Python environment says where it is (`log_file_begin`, `log_file_to`; `python -m planetgen.util.logpaths update.log` prints the path).
- All four are `*.log` files of the folder, so the one logrotate rule covers them; the macOS newsyslog file lists the three new ones. A log that cannot be written is said once and never stops the job. `smoke` keeps its JSON output file. The maintenance units needed no change (journald still holds their output).
- `test_settings.py` names `PLANETGEN_LOG_TIMES` (OPS.98's switch), which the check of variables the code reads had not been told about, and the generated settings files carry the new `log_dir` text.

generation-output: unchanged
