### Added
- **Windows installers.** `install.ps1` and `update.ps1` do what
  `install.sh` and `update.sh` do, with the Windows guide's layout: a venv
  with the libraries and waitress from `requirements-server.lock`
  (checked by hash), the NLTK corpus with `NLTK_DATA` set machine-wide,
  `config.json` from the Windows example when there is none, the
  migrate-or-delete question (y/N, 30 seconds), the runtime folders, and
  `icacls` permissions for the app's account. `update.ps1` pulls and
  installs only what's missing. `examples/maintenance/install-maintenance-task.ps1`
  schedules the monthly orbit update and `update.ps1` with Task
  Scheduler.
- **The bash installers run on macOS.** `install.sh`, `update.sh`,
  `install-maintenance-timer.sh` and their helpers run under macOS's
  bash 3.2: a venv from Homebrew's `python3` with gunicorn, `_www`
  ownership, `newsyslog` for the debug log, and gunicorn, the orbit update
  and `update.sh` as launchd daemons (`examples/macos/org.planetgen.update.plist`
  is new). `install.sh --skip-database` leaves out the database step.
- **`requirements-server.lock`** and a `server` extra in `setup.py`
  (gunicorn, or waitress on Windows) for those venvs.
- CI runs `install.sh` on macOS and `install.ps1` on Windows.
