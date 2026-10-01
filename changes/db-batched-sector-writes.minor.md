### Changed
- **Sectors save about 2.5 times faster, with 30 times fewer database
  statements (PERF.12, PERF.13).** A filled sector's rows are now written
  as one multi-row INSERT per table instead of one INSERT per star,
  planet, moon, belt, comet and child row, and generation checks the
  schema once per process instead of on every connection. Measured on a
  40-system sector (median of five saves): 1.35 s and about 5,650
  statements before, 0.54 s and about 186 after. New ids come from a
  small `id_blocks` table (schema v45), so each row knows its parent's id
  before anything is written.
- **Several sector writers can run at once (PERF.14).** A sector reserves
  all its system and phenomenon names in one locked statement instead of
  a `SELECT ... FOR UPDATE` per name, which deadlocked when four
  `generate.py sector` runs started together. A sector save now runs at
  READ COMMITTED and is retried from the start on a deadlock or lock
  wait timeout.

### Fixed
- **A deadlock no longer leaves a half-saved transaction carrying on.**
  The connection pool used to re-run a statement that hit a deadlock on
  a fresh session, so the rest of the save continued in a transaction
  MySQL had already rolled back and failed later with a foreign key
  error (1452). The pool now raises the error instead, and the save
  retries cleanly.

Run `update.sh` (or `update.ps1`) after updating: it migrates the
database to schema v45.
