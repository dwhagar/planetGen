**Headline:** The job engine has its ledger tables: a run, its jobs and every attempt are rows in the control database, and each galaxy holds the table a job's done mark commits into.

### Changed
- **Ledger tables and the module that reads and writes them (OPS.215).** The control database gains `engine_runs`, `engine_jobs` and `engine_attempts` (control schema v23); each galaxy gains `engine_items` (schema v260, migration 0260), the row a job's own commits update `WHERE attempt = mine` so its data and its done mark are one transaction. `planetgen/queue/ledger.py` adds a run and its jobs (repeating a plan adds nothing), raises an attempt at dispatch and writes it to the galaxy before the job starts, closes an attempt with how it ended and what it cost, and holds the fence (`commit_item` raises `Replaced` for a job another attempt took). A `Writer` runs ledger calls on a thread of its own so the supervisor's loop never waits on MySQL. Nothing reads or writes the rows yet: the supervisor (OPS.207) and the plan-as-rows item (OPS.221) build on it, so no run changes. A galaxy made before v260 gets the empty table at its next migrate.

generation-output: unchanged
