**Headline:** The code, the tests and the docs now agree that planetGen runs on MySQL 8.0.19 or later and not on MariaDB.

### Removed
- **The MariaDB branches.** A statement time limit is set with MySQL's `MAX_EXECUTION_TIME` only (the MariaDB
  `max_statement_time` and its error 1969 are gone), the Generate estimate measures the database disk from the
  server's data directory only (the `information_schema.DISKS` fallback was MariaDB's), and the table-size reads no
  longer guard against a server without `information_schema_stats_expiry`.
- **MariaDB from the tests, the CI legs and the docs.** The reserved-word list, the CHECK constraint helpers, the
  disk tests and the pool cache no longer name MariaDB; the two MariaDB legs of the CI matrix run on MySQL 8.0; the
  `mysql.host` help text, `docs/config.md`, `docs/cli.md`, `docs/testing.md` and `docs/ci-runners.md` say MySQL.
