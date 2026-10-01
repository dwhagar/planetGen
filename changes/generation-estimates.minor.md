### Added

- Every bulk generation now shows its size and time before it writes
  anything, and refuses one the database disk can't hold (PERF.3).
  `generate.py galaxy` (every mode) and `generate.py sector
  --num-sectors` print the expected sectors, star systems, database
  growth (+10%) and time, and on a terminal ask `Generate these N
  sectors? [y/N]` (`--yes` skips it). A run that would take more than a
  quarter of the database disk, or leave under 5 GB free, is refused
  with what it needs. `--estimate-only` prints the estimate and stops.
  The Generate page, the Galaxy Map's and Sector Map's Generate buttons,
  and a sector page's "Generate more sectors around this one" show the
  estimate and ask before starting.
- The server records how fast it generates, per star density (PERF.10):
  every filled sector adds its time to a log-scale density bucket (two
  per decade from 0.01) as a decaying average, and each galaxy's bytes
  per star system are measured after every run. The estimates use them;
  the admin Stats page and `GET /api/admin/generation-stats` show them.
  Control schema v6 (`generation_stats`, `generation_size`): run
  `update.sh`.
