**Headline:** Every cache warm-up now writes its own record to `warmup.log` in the log folder, in local time with its zone.

### Changed

- `planetgen warm` appends each run to `warmup.log` itself (next to `planetgen.log`, in `log_dir`): a line that says what started it (update stage 1 or 2, by hand, or a change such as an orbit update) and the time limit, the lines the run printed, each with the local time and its zone, and a last line that says it finished, was ended (by a signal or Ctrl+C) or failed, and after how long. Runs are added one after another, and the background warm-up that starts after a change, which wrote nowhere before, is in it too. The file is rotated with the other logs of the folder.
- `update.sh` no longer pipes the warm-up through `tee`, so the progress bars draw on its terminal as its comments say, and the old overwritten `/tmp/planetgen-warm-map.log` is gone. A warm-up whose log file cannot be written says so once and goes on.
