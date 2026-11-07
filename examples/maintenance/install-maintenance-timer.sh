#!/usr/bin/env bash
#
# examples/maintenance/install-maintenance-timer.sh
#
# Sets up planetGen's periodic maintenance as systemd timers -- the
# "once a month or so" orbit update planetgen.cli.orbits's own header comment
# and docs/reference/database-schema.md's cron example call for, plus (by default)
# the same `sudo ./update.sh` a human would otherwise have to remember to
# run -- packaged the Ubuntu/Debian-native way (systemd units + journald
# logging) instead of hand-added crontab lines and logfiles.
#
# Usage:
#   sudo examples/maintenance/install-maintenance-timer.sh [--skip-update-timer] [--skip-update-button] [--skip-tiles-timer] [--with-smoke-timer] [database ...]
#
# Each remaining argument is a database name to enable a monthly orbit
# update timer for (a deployment with more than one game database -- see
# `PLANETGEN_MYSQL_DATABASE_PREFIX` in docs/reference/database-schema.md -- gets one
# timer instance per database). Defaults to $PLANETGEN_MYSQL_DATABASE, or
# "planetgen" if that's unset too, matching planetgen.cli.orbits's own default.
#
# --skip-update-timer: only install the orbit-update timer(s), not
# planetgen-update.timer. Use this if this deployment's branch should
# only ever be updated by a human deliberately running `sudo ./update.sh`
# -- see planetgen-update.service's own header comment for what
# unattended updates mean (whatever's on the tracked branch gets pulled
# and deployed automatically, with no human review in between).
#
# --skip-update-button: do not install planetgen-update.path (SEC.43), the root-owned unit that runs update.sh when the admin
# pages' Update button drops <state folder>/update.request. The button never gets a root shell: the web user only creates that
# empty file, and the unit starts planetgen-update.service, whose command takes nothing from the file. install.sh and update.sh
# keep the unit's paths right; this installs it for a host that only ran this script (Linux (systemd) only).
#
# --skip-tiles-timer: do not install planetgen-tiles@.service/.timer (PERF.233), the nightly rebuild of the whole Galaxy Map
# tile cache at 05:00, after the orbit update (`planetgen.cli.tiles`; Linux only). It needs room on the cache's disk for the new
# cache beside the live one.
#
# --with-smoke-timer: also install planetgen-smoke.service/.timer (PERF.183), the nightly New galaxy smoke run at 04:30
# that fails when a step of a small New galaxy peaks above its memory budget. Linux (systemd) only; the account in
# /etc/planetgen/maintenance.env must be allowed to create and drop databases.
#
# This script and install.sh/update.sh share one function, `install_maintenance_units` (scripts/deploy-common.sh): every unit is
# rendered for the checkout this script lives in (OPS.110). The example files say /var/lib/planetGen and
# /usr/local/planetgen/venv, and the installed ones say where the checkout and the virtual environment (PLANETGEN_VENV, as
# install.sh and update.sh read it) really are, and the state folder deploy-paths.py gives. A path with a space or a special
# character is refused. A unit already installed that differs from the rendered one (an older copy of the example, or another
# location) is replaced and the replacement is printed; edit the example, or a drop-in under /etc/systemd/system/<unit>.d/,
# for a setting of your own. install.sh and update.sh run the same function on every deploy, so they put the timers in place on
# their own (a deploy enables only a unit it installs for the first time, and never the unattended update timer, which it installs
# and leaves off with a line saying how to turn it on); this script is for choosing the databases and the options, and it enables
# everything it was asked for, the update timer too. The gunicorn unit (examples/systemd/planetgen-gunicorn.service) is not
# installed by either: copy it by hand, changing its two paths (see its own header).
#
# What it does:
#   1. Installs planetgen-orbits@.service/.timer (always), planetgen-update.service/.timer (unless --skip-update-timer),
#      planetgen-update.path (unless --skip-update-button), the tiles units (unless --skip-tiles-timer) and, with
#      --with-smoke-timer, the smoke units, into /etc/systemd/system/.
#   2. Installs /etc/planetgen/maintenance.env (mode 600, since it holds a
#      DB password) from maintenance.env.example, but only if it doesn't
#      already exist -- never overwrites credentials a previous run or a
#      human already set up here.
#   3. `systemctl daemon-reload`, then `enable --now` planetgen-update.timer
#      (unless skipped) and each requested database's orbit timer instance.
#
# On macOS, the same schedule as launchd daemons instead:
# /Library/LaunchDaemons/org.planetgen.orbits.<database>.plist per
# database (examples/macos/org.planetgen.orbits.planetgen.plist) and
# org.planetgen.update.plist (unless --skip-update-timer), with the
# checkout's path put in. No credentials file there: the orbit update
# reads config.json. An existing plist is left alone.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

install_update_timer=1
install_update_button=1
install_smoke_timer=0
install_tiles_timer=1
args=()
for arg in "$@"; do
    if [[ "$arg" == "--skip-update-timer" ]]; then
        install_update_timer=0
    elif [[ "$arg" == "--skip-update-button" ]]; then
        install_update_button=0
    elif [[ "$arg" == "--skip-tiles-timer" ]]; then
        install_tiles_timer=0
    elif [[ "$arg" == "--with-smoke-timer" ]]; then
        install_smoke_timer=1
    else
        args+=("$arg")
    fi
done

if [[ $EUID -ne 0 ]]; then
    echo "error: must be run as root, e.g.:" >&2
    echo "  sudo $0 $*" >&2
    exit 1
fi

if [[ "$(uname -s)" != Darwin ]] && ! command -v systemctl >/dev/null 2>&1; then
    echo "error: systemctl not found -- this script is for systemd-based" >&2
    echo "  systems (Ubuntu/Debian's default). Use docs/reference/database-schema.md's" >&2
    echo "  plain crontab example instead if this host has no systemd." >&2
    exit 1
fi

databases=(${args[@]+"${args[@]}"})
if [[ ${#databases[@]} -eq 0 ]]; then
    databases=("${PLANETGEN_MYSQL_DATABASE:-planetgen}")
fi

if [[ "$(uname -s)" == Darwin ]]; then
    REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
    MACOS_DIR="$REPO_DIR/examples/macos"
    PYTHON="${PYTHON:-python3}"
    SCRIPT_DIR="$REPO_DIR"  # what deploy-common.sh's helpers expect
    # shellcheck source=../../scripts/deploy-common.sh
    source "$REPO_DIR/scripts/deploy-common.sh"
    if [[ "$install_smoke_timer" -eq 1 ]]; then
        echo "  --with-smoke-timer is for systemd only; run 'python -m planetgen.cli.smoke' by hand on macOS."
    fi
    mkdir -p "$MACOS_LOG_DIR"
    chown _www:_www "$MACOS_LOG_DIR"
    echo "== Installing launchd daemons =="
    if [[ "$install_update_timer" -eq 1 ]]; then
        install_launchd_plist "$MACOS_DIR/org.planetgen.update.plist" org.planetgen.update.plist
        echo "  org.planetgen.update runs update.sh at 03:00 on the 1st of each month"
    else
        echo "  --skip-update-timer given -- not installing org.planetgen.update"
    fi
    for db in "${databases[@]}"; do
        install_launchd_plist "$MACOS_DIR/org.planetgen.orbits.planetgen.plist" "org.planetgen.orbits.$db.plist" \
            -e "s|org\.planetgen\.orbits\.planetgen|org.planetgen.orbits.$db|" \
            -e "s|<string>planetgen</string>|<string>$db</string>|"
        echo "  org.planetgen.orbits.$db updates orbits at 03:30 on the 1st of each month"
    done
    echo
    echo "Done. Check them with: sudo launchctl print system/org.planetgen.orbits.${databases[0]}"
    echo "Run one now with:      sudo launchctl kickstart system/org.planetgen.orbits.${databases[0]}"
    echo "Logs: $MACOS_LOG_DIR"
    exit 0
fi

REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
PYTHON="${PYTHON:-python3}"
[[ -x "${PLANETGEN_VENV:-/usr/local/planetgen/venv}/bin/python" ]] && PYTHON="${PLANETGEN_VENV:-/usr/local/planetgen/venv}/bin/python"
SCRIPT_DIR="$REPO_DIR"  # what deploy-common.sh's helpers expect
# shellcheck source=../../scripts/deploy-common.sh
source "$REPO_DIR/scripts/deploy-common.sh"

options=(--enable-all --credentials)
[[ "$install_update_timer" -eq 1 ]] || options+=(--skip-update-timer)
[[ "$install_update_button" -eq 1 ]] || options+=(--skip-update-button)
[[ "$install_tiles_timer" -eq 1 ]] || options+=(--skip-tiles-timer)
[[ "$install_smoke_timer" -eq 0 ]] || options+=(--with-smoke-timer)
echo "== Installing systemd units for $REPO_DIR (virtual environment $VENV) =="
install_maintenance_units "${options[@]}" "${databases[@]}" || { echo "error: the units were not all installed (see above)." >&2; exit 1; }

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
