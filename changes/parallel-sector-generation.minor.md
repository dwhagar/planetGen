### Added
- **Sectors are generated several at a time (PERF.8).** `generate.py
  sector` and every `generate.py galaxy` mode, from the command line or
  the Generate page, fill sectors in parallel worker processes: by
  default 80% of the machine's cores (one fewer when MySQL runs on the
  same machine), at a lower priority than everything else on it.
  `--workers N` (or `PLANETGEN_WORKERS`) sets the count, and
  `--workers 1` works one sector at a time as before. On a 4-core
  machine with MySQL local, 60 sectors took 9.8 s instead of 16.5 s.
  Only one run's workers use the machine at a time: a second run waits
  for the first, using a lease in the control database (schema v5,
  `work_jobs`/`work_tasks`/`work_lease`).

### Changed
- **Linking a sector to its neighbors is about four times faster.**
  The nearest-systems search no longer walks a mostly empty grid for
  every object next to a newly saved sector; five neighboring sectors
  went from 13.9 s to 3.4 s.

Run `update.sh` (or `update.ps1`) after updating: it adds the work
queue's tables to the control database.
