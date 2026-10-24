**Headline:** A work unit no longer runs the whole database schema again when it connects.

### Changed
- **A connection skips the schema statements when the database is already at the current version (PERF.174).** Each queued work unit is a new process, and the first connection of every process ran `schema.sql` (ten-odd `CREATE TABLE IF NOT EXISTS` statements and the views, under the schema lock) to find nothing to do: 861 times, about 130 s of server time, in the 617 s band step of the benchmark suite, and about 32,000 such statements in one four-minute scatter. A connection now reads the stored schema version first and applies the schema only when it is not the code's own. A new database, and one behind the code, are handled as before.

TODO items resolved: PERF.174. Added: none.
