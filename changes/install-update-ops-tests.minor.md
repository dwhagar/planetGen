### Added

- **Install and update set up and check both logs on every OS (OPS.5).**
  On Linux and macOS, `install.sh` and `update.sh` (through
  `examples/apache/setup-debug-log.sh` and the new
  `examples/apache/log-locations.py`) now prepare the debug log whether or
  not `debug` is on, and the activity log's folder, wherever
  `PLANETGEN_LOG_FILE`/`log_file` and `PLANETGEN_LOG_DIR`/`log_dir` put
  them, and check that the web server's user and its group can really
  write each one. On Windows, `install.ps1` and `update.ps1` check that the
  app's account can write both logs' folders. Anything they can't set up
  (no rights, a read-only or missing drive, a folder in the way that is a
  file, a user that doesn't exist yet) only warns, with the exact
  `mkdir`/`chown`/`chmod` or `New-Item`/`icacls` commands that fix it, or
  how to point the setting somewhere writable; a log never stops an
  install or update.
- **CI runs `install.sh` and `update.sh` against a live database (TEST.62).**
  A new `linux-update` job runs both for real on Linux: a fresh install,
  an update with nothing new, a database that needs migrating, one newer
  than the code, an unreachable server and a failed migration. The shared
  database step is also tested on every database engine in pytest.
- **Command-line tests for the admin scripts (TEST.60):** `queryDb`,
  `adminStats`, `checkRenderParity`, `dedupeNames`, `resetDb`,
  `updateOrbits` and `loginLockouts`, including a bad port, an unknown
  database and an empty password against the environment for each.

### Changed

- **A database newer than the code is refused (TEST.10).** `migrateDb.py`
  (and so install and update), `migrateDb.py --status` and every
  read-write connection now stop with a clear message when the galaxy
  database or the control schema is at a higher version than this code
  knows, instead of carrying on (and before this code's older
  `schema.sql` touches it). `/api/health` says the database is newer than
  the code rather than telling you to run `migrateDb.py`.
- **`migrateDb.py --status` no longer changes anything:** it reads the
  version without laying down the schema, so it works with a read-only
  account, and reports a new empty database as current.

### Fixed

- `set-permissions.sh`, `create-cache-dir.sh` and `setup-debug-log.sh`
  stopped silently (exit 1, no message) on a Debian or Ubuntu server with
  Apache installed: reading `/etc/apache2/envvars` under `set -u` hit its
  unset `$APACHE_CONFDIR` and killed the user/group lookup. It is read
  safely now, and a lookup that still fails says so.

- `checkRenderParity.py`, `dedupeNames.py`, `resetDb.py` and
  `loginLockouts.py` print `error: ...` and exit 1 on a database they
  can't reach instead of a traceback.
- `updateOrbits.py` no longer fails every run with "elapsed_years must be
  >= 0" after the database server's clock went back; it moves nothing and
  restarts the clock from now.
- `loginLockouts.py --ip ""` is a usage error instead of quietly listing
  the lockouts.

### Removed

- `src/migrateSqliteToMysql.py` (TEST.61). It only accepted a SQLite file
  already at the current schema version, which no SQLite database ever
  reached (SQLite stopped at v5 and the MySQL migrations start at v8), so
  it could never import anything.
