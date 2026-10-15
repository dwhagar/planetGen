### Added
- A run now has one record in Redis: its stages, how many units each has finished, an event stream that serves as the live log, and an estimate of the time left. When a phase reports no updates, the estimate comes from earlier runs' stage timings, or from the time per task so far in this run, so the whole bar always shows a time remaining.
- A public view of the active run (stage, percent, time left, whether the time comes from earlier runs) for pages to show.
