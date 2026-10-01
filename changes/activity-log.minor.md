### Security
- **An always-on activity log (SEC.28).** planetGen now always writes a
  short record of who did what to `planetgen.log` in the platform's
  standard log folder (`/var/log/planetgen/` on Linux,
  `/Library/Logs/planetgen/` on macOS, `logs\` under the checkout on
  Windows; `"log_dir"` in `config.json` moves it): sign-ins, logouts,
  credential and API key changes, refused requests (missing, expired or
  forged sessions, bad API keys, admin-only actions, failed form tokens),
  every database write made through the web interface or API, each
  schema migration step, and each `generate.py` run or Generate page job
  as one start and one finish line with its counts. One fixed line format
  (`2026-10-01T08:00:00Z planetgen[1234]: AUTH login.failed ip=203.0.113.5 user="admin"`),
  documented in `docs/config.md`, with the address always before any
  text a visitor typed, which is quoted and escaped. It is separate from
  the debug log, which keeps working as before and gets a copy of each
  line while it is on. `install.sh`/`update.sh` create the folder and
  install logrotate (`/etc/logrotate.d/planetgen-log`) or newsyslog
  (`/etc/newsyslog.d/planetgen-log.conf`) rotation, daily or past 100 MB,
  30 copies; on Windows, or without those files, the program rotates the
  file itself.
- **Failed and locked logins are recorded with their address (SEC.20).**
  Each one is an `AUTH login.failed`/`login.locked` line in the activity
  log and a row in the control database's audit log (kept 90 days), as is
  a wrong current password on `/account`. The admin Stats page lists the
  newest ones under "Failed sign-ins" (new `GET /api/admin/login-failures`).
  A wrong password on the `/login` page now answers 401 instead of 200,
  so the web server's access log shows it too.
