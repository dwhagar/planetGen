#!/usr/bin/env bash
#
# examples/maintenance/install-maintenance-timer.sh
#
# Sets up the periodic `src/updateOrbits.py` maintenance job as a systemd
# timer -- the "once a month or so" cadence updateOrbits.py's own header
# comment and docs/database-schema.md's cron example both call for,
# packaged the Ubuntu/Debian-native way (systemd units + journald logging)
# instead of a hand-added crontab line.
#
# Usage:
#   sudo examples/maintenance/install-maintenance-timer.sh [database ...]
#
# Each argument is a database name to enable a monthly update timer for
# (a deployment with more than one game database -- see
# `PLANETGEN_MYSQL_DATABASE_PREFIX` in docs/database-schema.md -- gets one
# timer instance per database). Defaults to $PLANETGEN_MYSQL_DATABASE, or
# "planetgen" if that's unset too, matching updateOrbits.py's own default.
#
# What it does:
#   1. Installs planetgen-orbits@.service/.timer into /etc/systemd/system/.
#   2. Installs /etc/planetgen/maintenance.env (mode 600, since it holds a
#      DB password) from maintenance.env.example, but only if it doesn't
#      already exist -- never overwrites credentials a previous run or a
#      human already set up here.
#   3. `systemctl daemon-reload`, then `enable --now` each requested
#      database's timer instance.
#
# Linux only (systemd) -- same scope as install.sh/update.sh, which don't
# call this themselves since the orbit update runs on its own schedule,
# not on every deploy (see docs/database-schema.md).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $EUID -ne 0 ]]; then
    echo "error: must be run as root, e.g.:" >&2
    echo "  sudo $0 $*" >&2
    exit 1
fi

if ! command -v systemctl >/dev/null 2>&1; then
    echo "error: systemctl not found -- this script is for systemd-based" >&2
    echo "  systems (Ubuntu/Debian's default). Use docs/database-schema.md's" >&2
    echo "  plain crontab example instead if this host has no systemd." >&2
    exit 1
fi

databases=("$@")
if [[ ${#databases[@]} -eq 0 ]]; then
    databases=("${PLANETGEN_MYSQL_DATABASE:-planetgen}")
fi

echo "== 1/3: Installing systemd unit files =="
install -m 644 "$SCRIPT_DIR/planetgen-orbits@.service" /etc/systemd/system/planetgen-orbits@.service
install -m 644 "$SCRIPT_DIR/planetgen-orbits@.timer" /etc/systemd/system/planetgen-orbits@.timer
echo "  installed planetgen-orbits@.service and planetgen-orbits@.timer"

echo
echo "== 2/3: Setting up credentials file =="
ENV_FILE=/etc/planetgen/maintenance.env
if [[ -f "$ENV_FILE" ]]; then
    echo "  $ENV_FILE already exists -- leaving it alone."
else
    mkdir -p /etc/planetgen
    install -m 600 "$SCRIPT_DIR/maintenance.env.example" "$ENV_FILE"
    echo "  wrote $ENV_FILE (mode 600) from maintenance.env.example -- edit it"
    echo "  now with this database's read-write MySQL credentials before the"
    echo "  timer's first scheduled run (defaults match a freshly installed"
    echo "  MySQL server's typical planetgen account, not a production one)."
fi

echo
echo "== 3/3: Enabling the timer for: ${databases[*]} =="
systemctl daemon-reload
for db in "${databases[@]}"; do
    systemctl enable --now "planetgen-orbits@${db}.timer"
    echo "  enabled planetgen-orbits@${db}.timer"
done

echo
echo "Done. Check scheduled runs with:"
echo "  systemctl list-timers 'planetgen-orbits@*'"
echo "Check a specific database's last run with, e.g.:"
echo "  journalctl -u planetgen-orbits@${databases[0]}.service"
echo "Run one immediately (without waiting for the schedule) with, e.g.:"
echo "  sudo systemctl start planetgen-orbits@${databases[0]}.service"
