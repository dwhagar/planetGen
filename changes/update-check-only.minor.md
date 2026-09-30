### Changed
- **`update.sh` no longer reinstalls anything.** It used to re-run
  `install.sh` after every pull that brought new commits, which
  force-reinstalled the Python package with pip, rebuilt the fallback
  venv and re-fetched the NLTK corpus. Now it checks instead: each Python
  library is imported with the system Python and compared with
  `setup.py`'s floor (`scripts/install-python-deps.sh --check`), and only
  one that is missing, too old or broken is installed, the same way
  `install.sh` would on that host. The NLTK corpus and mod_wsgi are
  installed, and Apache's modules enabled, only when missing. It prints
  one line per library (present, installed, upgraded, repaired or failed)
  with where it came from, and finishes by importing the web app as
  Apache's user, so a library www-data can't use fails the update instead
  of the site.
- **No more venv: libraries go into the system Python.** On an
  externally managed Python (Ubuntu 24.04+), everything apt packages at
  or above `setup.py`'s floor comes from apt. Only what apt lacks or ships
  too old is pip-installed system-wide into `/usr/local`, alongside apt's
  copy and never over it. The report says which ones and why. An existing
  `/opt/planetgen/venv` and its `planetgen-venv.pth` are removed on the
  next run, and their libraries are installed system-wide.
- The NLTK and Apache module steps now live in `scripts/deploy-common.sh`,
  shared by `install.sh` and `update.sh`, and `install.sh` installs
  mod_wsgi when apt can instead of only warning about it.
- `/usr/local/bin/planetgen` is now the checkout wrapper on every host
  (pip's console script ran the pip-installed copy, which would go stale
  once updates stopped reinstalling it).
- **Checking which Python Apache really uses.** `install.sh`/`update.sh`
  now warn when mod_wsgi is built for a different Python version than the
  one they set the libraries up for, and both take `PYTHON=` to pick
  another interpreter. The admin Stats page shows the web app's Python
  prefix and the directory it imports its libraries from. See
  `docs/apache-deployment.md`.

### Added
- **Migrate or delete the database on update.** When the database is
  behind the current schema, `update.sh` (and `install.sh`) asks whether
  to delete the galaxy data instead of migrating it: y/N, with a
  30-second timeout that defaults to N (keep the data and migrate it).
  The same happens with no terminal to ask on. Deleting wipes every
  generated sector and system (`resetDb.py`) and keeps admin logins.
- **Progress bar and ETA for database migrations.** `migrateDb.py` shows
  each migration step as it runs, with the elapsed time and an estimate
  of the time left. `migrateDb.py --status` reports the current and
  target schema versions without changing anything.
