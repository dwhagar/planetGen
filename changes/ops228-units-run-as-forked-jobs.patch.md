**Headline:** The units of a generation run are now handed out by the supervisor's launcher and run as jobs forked from a warm template, instead of going to RQ workers.

### Changed
- `queue.forked_jobs` (environment `PLANETGEN_FORKED_JOBS`) is on by default. A run whose stages are cut into ledger rows opens each stage by adding the arguments to its rows and releasing them to the launcher, and the launcher forks each unit, begins and ends its attempt, and settles the row. Off, every unit goes to a worker as before (OPS.228 slice D1).
- A run carried on after an update hands its remaining units to the launcher the same way.
