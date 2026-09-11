#!/usr/bin/env bash
#
# examples/maintenance/install-maintenance-timer.sh
#
# Sets up planetGen's periodic maintenance as systemd timers -- the
# "once a month or so" orbit update updateOrbits.py's own header comment
# and docs/database-schema.md's cron example call for, plus (by default)
# the same `sudo ./update.sh` a human would otherwise have to remember to
# run -- packaged the Ubuntu/Debian-native way (systemd units + journald
# logging) instead of hand-added crontab lines and logfiles.
#
# Usage:
#   sudo examples/maintenance/install-maintenance-timer.sh [--skip-update-timer] [database ...]
#
# Each remaining argument is a database name to enable a monthly orbit
# update timer for (a deployment with more than one game database -- see
# `PLANETGEN_MYSQL_DATABASE_PREFIX` in docs/database-schema.md -- gets one
# timer instance per database). Defaults to $PLANETGEN_MYSQL_DATABASE, or
# "planetgen" if that's unset too, matching updateOrbits.py's own default.
#
# --skip-update-timer: only install the orbit-update timer(s), not
# planetgen-update.timer. Use this if this deployment's branch should
# only ever be updated by a human deliberately running `sudo ./update.sh`
# -- see planetgen-update.service's own header comment for what
# unattended updates mean (whatever's on the tracked branch gets pulled
# and deployed automatically, with no human review in between).
#
# What it does:
#   1. Installs planetgen-orbits@.service/.timer (always) and
#      planetgen-update.service/.timer (unless --skip-update-timer) into
#      /etc/systemd/system/.
#   2. Installs /etc/planetgen/maintenance.env (mode 600, since it holds a
#      DB password) from maintenance.env.example, but only if it doesn't
#      already exist -- never overwrites credentials a previous run or a
#      human already set up here.
#   3. `systemctl daemon-reload`, then `enable --now` planetgen-update.timer
#      (unless skipped) and each requested database's orbit timer instance.
#
# Linux only (systemd) -- same scope as install.sh/update.sh, which don't
# call this themselves since maintenance runs on its own schedule, not on
# every deploy.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

install_update_timer=1
args=()
for arg in "$@"; do
    if [[ "$arg" == "--skip-update-timer" ]]; then
        install_update_timer=0
    else
        args+=("$arg")
    fi
done

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

databases=("${args[@]}")
if [[ ${#databases[@]} -eq 0 ]]; then
    databases=("${PLANETGEN_MYSQL_DATABASE:-planetgen}")
fi

echo "== 1/3: Installing systemd unit files =="
install -m 644 "$SCRIPT_DIR/planetgen-orbits@.service" /etc/systemd/system/planetgen-orbits@.service
install -m 644 "$SCRIPT_DIR/planetgen-orbits@.timer" /etc/systemd/system/planetgen-orbits@.timer
echo "  installed planetgen-orbits@.service and planetgen-orbits@.timer"
if [[ "$install_update_timer" -eq 1 ]]; then
    install -m 644 "$SCRIPT_DIR/planetgen-update.service" /etc/systemd/system/planetgen-update.service
    install -m 644 "$SCRIPT_DIR/planetgen-update.timer" /etc/systemd/system/planetgen-update.timer
    echo "  installed planetgen-update.service and planetgen-update.timer"
else
    echo "  --skip-update-timer given -- not installing planetgen-update.service/.timer"
fi

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
echo "== 3/3: Enabling timers =="
systemctl daemon-reload
if [[ "$install_update_timer" -eq 1 ]]; then
    systemctl enable --now planetgen-update.timer
    echo "  enabled planetgen-update.timer (fires 03:00 on the 1st of each month)"
fi
for db in "${databases[@]}"; do
    systemctl enable --now "planetgen-orbits@${db}.timer"
    echo "  enabled planetgen-orbits@${db}.timer (fires 03:30 on the 1st of each month)"
done

echo
echo "Done. Check scheduled runs with:"
echo "  systemctl list-timers 'planetgen-*'"
if [[ "$install_update_timer" -eq 1 ]]; then
    echo "Check the last deployment update with:"
    echo "  journalctl -u planetgen-update.service"
fi
echo "Check a specific database's last orbit update with, e.g.:"
echo "  journalctl -u planetgen-orbits@${databases[0]}.service"
echo "Run one immediately (without waiting for the schedule) with, e.g.:"
echo "  sudo systemctl start planetgen-orbits@${databases[0]}.service"
