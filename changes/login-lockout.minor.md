### Security
- **Per-address login lockout (SEC.1).** Three failed logins from one
  address lock it for 5 minutes; each further lockout doubles that, up
  to a day. Checked before the password, so a guess during a lock learns
  nothing. A successful login clears the address's count, but its
  doubling level only halves for each day without a lockout. IPv6
  addresses are counted by their /64. Loopback and the new
  `login_allowlist` in `config.json` are never locked. If a private
  address gets locked while `proxy_fix` is off, the app warns in the
  error log and on the admin Stats page, because behind a reverse proxy
  every visitor shares that address.
- **Lockouts are shared and survive restarts (SEC.21).** The
  per-username backoff (10 free failures, then 1 s doubling to 15
  minutes) and the new per-address lockout live in the control
  database's new `login_throttle` table (control schema v2) instead of
  each worker's memory. Until `update.sh` creates the table, the counts
  are kept in memory as before.
- **Lifting a lockout.** The admin Stats page lists current lockouts
  with Lift buttons (`GET /api/admin/lockouts`, `POST
  /api/admin/lockouts/lift`), and `python3 src/loginLockouts.py` lists or
  lifts them from a shell (`--ip`, `--user`, `--all`) for an admin locked
  out of the site itself.
