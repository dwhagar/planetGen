**Headline:** The admin queue page shows the cache warm-up phase by phase: which are waiting, running, built, unchanged or stopped, with each one's result.

### Changed
- **A "Cache warm-up" panel on `/admin/queue` (PERF.141).** It shows the last warm-up run (its state, when it started and ended, what asked for it), whether a start after a change is waiting or running and for which changes, and one row for each of the five phases with its state and the line the run printed (what it built, what was already stored, what the budget left out). A run that died reads as interrupted, and a failed run marks the phase that failed. The same state is at `GET /api/admin/work/warm?database=NAME` for admins.

TODO items resolved / Added: PERF.141 (slice 3: the admin page shows each phase) / Added: none
