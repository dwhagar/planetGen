#!/usr/bin/env bash
#
# examples/apache/setup-debug-log.sh
#
# Prepares planetGen's debug log (see src/planetgen/util/log.py), its
# always-on activity log (src/planetgen/admin/activity_log.py), and the
# rotation for both. Safe to run again. Called by install.sh, and by
# update.sh when there's nothing new to install.
#
#   1. Sets up the debug log ("log_file" in config.json, or
#      PLANETGEN_LOG_FILE; default /var/log/planetgen.log) whether or not
#      debug is on, so turning it on later needs nothing else: its folder
#      is made if it's missing and the file is created, owned by Apache's
#      worker user and group, mode 0660 (log-locations.py). The web
#      interface writes to it as that user, root (a systemd timer, sudo)
#      can always write to it, and a login user who runs the generator
#      from a shell must be in the Apache group to append to it
#      (`sudo usermod -aG <apache group> <user>`, then log in again; on
#      macOS `sudo dseditgroup -o edit -a <user> -t user _www`).
#      Other users can neither read nor write it -- a world-writable log
#      would let any local user forge or flood its entries.
#   2. Checks that Apache's user, and its group, can really write each
#      log (every folder above it must let them through). Anything that
#      can't be set up -- no rights, a read-only or missing drive, a user
#      that doesn't exist yet -- only warns, with the exact commands that
#      fix it or how to point "log_file"/"log_dir" elsewhere (OPS.5): a
#      log never stops the install or update.
#   3. Installs /etc/logrotate.d/planetgen, so a debug log someone forgot to
#      turn off can't fill the disk: rotated daily, or as soon as it passes
#      100 MB, keeping 7 compressed copies. Where /etc/cron.hourly exists it
#      also adds a job that checks the size every hour, since a galaxy run
#      with debug on can write gigabytes between daily logrotate runs.
#      On macOS, which has newsyslog instead of logrotate, installs
#      /etc/newsyslog.d/planetgen.conf: rotated past 100 MB, 7 compressed
#      copies kept (newsyslog runs hourly by itself).
#   4. Creates the activity log's folder ("log_dir" in config.json, default
#      /var/log/planetgen on Linux, /Library/Logs/planetgen on macOS),
#      owned by root and Apache's group, mode 2770 (the setgid bit keeps
#      new files in that group), and the log file itself, mode 0660 owned
#      by Apache's user and group: the same "Apache and its group write,
#      nobody else reads" rule as the debug log (log-locations.py again,
#      with the same checks and warnings). Installs its rotation:
#      /etc/logrotate.d/planetgen-log (daily or past 100 MB, 30 kept) on
#      Linux, /etc/newsyslog.d/planetgen-log.conf (daily or past 100 MB,
#      30 kept) on macOS. While that file exists the program leaves
#      rotation to the system (logpaths.log_rotation_mode).
#
# Rotation files that can't be written only warn too. Exits 0 whatever
# happens to the logs; only running it without root is an error.
#
# Runs on Linux and macOS. Windows has no rotation for it; install.ps1
# and update.ps1 check the log folders themselves
# (Test-LogLocations in scripts/deploy-common.ps1).
#
# Usage:
#   sudo examples/apache/setup-debug-log.sh

set -euo pipefail

APACHE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$APACHE_DIR/../.." && pwd)"

if [[ $EUID -ne 0 ]]; then
    echo "error: must be run as root (writes to /var/log and /etc/logrotate.d), e.g.:" >&2
    echo "  sudo $0" >&2
    exit 1
fi

PYTHON="${PYTHON:-$(command -v python3 || command -v python || true)}"
if [[ -z "$PYTHON" ]]; then
    echo "error: no python3/python found on PATH." >&2
    exit 1
fi

# Prints "<debug on: 1|0> <log file path>", read the same way the program
# itself reads them (planetgen/util/logpaths.py). logpaths is loaded
# straight from its file, not through the planetgen package, so this
# works before the package's own dependencies are installed. This runs as
# root, so -I keeps the current directory, user site-packages and PYTHON*
# variables off sys.path (logpaths.py itself lives outside src/html, in
# the root-owned part of the checkout).
# (Through a temp file: bash 3.2, macOS's, can't parse a heredoc inside
# a process substitution.)
SETTINGS="$(mktemp)"
if ! (cd / && "$PYTHON" -I - "$REPO_DIR" > "$SETTINGS") <<'PY'
import importlib.util
import os
import sys

path = os.path.join(sys.argv[1], "src", "planetgen", "util", "logpaths.py")
spec = importlib.util.spec_from_file_location("planetgen_logpaths", path)
logpaths = importlib.util.module_from_spec(spec)
spec.loader.exec_module(logpaths)
config = logpaths.read_config_file()
print(1 if logpaths.debug_enabled(config) else 0, logpaths.log_file_path(config))
print(logpaths.activity_log_path(config))
PY
then
    rm -f "$SETTINGS"
    echo "warning: couldn't work out where the logs go (see above), so they and their rotation" >&2
    echo "  were not set up. Check config.json's \"log_file\" and \"log_dir\" (each a path, or leave" >&2
    echo "  them out for the defaults), then run sudo ./update.sh again." >&2
    exit 0
fi
{ read -r _debug_on LOG_FILE; read -r ACTIVITY_LOG; } < "$SETTINGS"
rm -f "$SETTINGS"
ACTIVITY_DIR="$(dirname "$ACTIVITY_LOG")"

# shellcheck source=examples/apache/apache-identity.sh
source "$APACHE_DIR/apache-identity.sh"
if ! read -r APACHE_USER APACHE_GROUP < <(detect_apache_group) || [[ -z "${APACHE_GROUP:-}" ]]; then
    echo "warning: couldn't work out Apache's user and group (see above), so the logs were not set up." >&2
    echo "  Install the web server (sudo apt install apache2), then run sudo ./update.sh again." >&2
    exit 0
fi

# Both log files, their folders, owners and modes, and whether Apache's
# user and group can really write them; it only ever warns.
"$PYTHON" -I "$APACHE_DIR/log-locations.py" "$REPO_DIR" "$APACHE_USER" "$APACHE_GROUP" || true

# Writes stdin to a rotation config file, mode 0644, or warns (a
# read-only /etc, say) and carries on. Returns 1 when it couldn't.
write_rotation_conf() {
    local target="$1" text
    text="$(cat)"
    if mkdir -p "$(dirname "$target")" 2>/dev/null \
        && printf '%s\n' "$text" > "$target" 2>/dev/null \
        && chmod 0644 "$target" 2>/dev/null; then
        return 0
    fi
    echo "warning: couldn't write $target, so the logs it covers won't be rotated." >&2
    echo "  Fix it by making $(dirname "$target") writable for root, then run sudo ./update.sh again." >&2
    return 1
}

# Rotation is set up whether or not debug is on right now, so turning it on
# later is covered.
if [[ "$(uname -s)" == Darwin ]]; then
    # The activity log: daily at midnight ($D0) or past 100 MB, 30 kept,
    # bzip2 (J), no process to signal (N): the program notices the move
    # and reopens the file itself (WatchedFileHandler).
    if write_rotation_conf /etc/newsyslog.d/planetgen-log.conf <<EOF
# planetGen activity log -- written by examples/apache/setup-debug-log.sh
# (install.sh/update.sh). Edits here are overwritten on the next run.
# logfilename            [owner:group]              mode count size(KB) when  flags
$ACTIVITY_LOG  $APACHE_USER:$APACHE_GROUP  660  30  102400  \$D0  JN
EOF
    then
        echo "Activity log rotation: /etc/newsyslog.d/planetgen-log.conf (daily or past 100 MB, 30 kept)"
    fi

    # newsyslog: owner:group, mode, copies kept, size in KB, no time
    # rotation, Z = gzip. It creates the new file with that owner and mode.
    if write_rotation_conf /etc/newsyslog.d/planetgen.conf <<EOF
# planetGen debug log -- written by examples/apache/setup-debug-log.sh
# (install.sh/update.sh). Edits here are overwritten on the next run.
# logfilename            [owner:group]              mode count size(KB) when  flags
$LOG_FILE  $APACHE_USER:$APACHE_GROUP  660  7  102400  *  Z
EOF
    then
        echo "Log rotation: /etc/newsyslog.d/planetgen.conf (past 100 MB, 7 kept)"
    fi
    exit 0
fi

# `missingok` makes it a no-op while there's no file.
# `create` recreates the file with the same owner and mode after each
# rotation; the program reopens it on its own (WatchedFileHandler).
if write_rotation_conf /etc/logrotate.d/planetgen <<EOF
# planetGen debug log -- written by examples/apache/setup-debug-log.sh
# (install.sh/update.sh). Edits here are overwritten on the next run.
$LOG_FILE {
    daily
    maxsize 100M
    rotate 7
    missingok
    notifempty
    compress
    delaycompress
    create 0660 $APACHE_USER $APACHE_GROUP
    su root root
}
EOF
then
    echo "Log rotation: /etc/logrotate.d/planetgen (daily or past 100 MB, 7 kept)"
fi

# The activity log: kept longer (30 copies) since it is the record of
# who signed in and what changed. `su` because the folder is writable by
# Apache's group, which logrotate otherwise refuses.
if write_rotation_conf /etc/logrotate.d/planetgen-log <<EOF
# planetGen activity log -- written by examples/apache/setup-debug-log.sh
# (install.sh/update.sh). Edits here are overwritten on the next run.
$ACTIVITY_DIR/*.log {
    daily
    maxsize 100M
    rotate 30
    missingok
    notifempty
    compress
    delaycompress
    dateext
    # Seconds in the name: the hourly size check can rotate twice a day.
    dateformat -%Y%m%d-%s
    su root $APACHE_GROUP
    create 0660 $APACHE_USER $APACHE_GROUP
}
EOF
then
    echo "Activity log rotation: /etc/logrotate.d/planetgen-log (daily or past 100 MB, 30 kept)"
fi

if [[ -d /etc/cron.hourly ]]; then
    if write_rotation_conf /etc/cron.hourly/planetgen-logrotate <<'EOF'
#!/bin/sh
# planetGen: rotate the debug and activity logs as soon as they pass their
# maxsize instead of waiting for the daily logrotate run. Written by
# examples/apache/setup-debug-log.sh.
command -v logrotate >/dev/null 2>&1 || exit 0
logrotate /etc/logrotate.d/planetgen
exec logrotate /etc/logrotate.d/planetgen-log
EOF
    then
        chmod 0755 /etc/cron.hourly/planetgen-logrotate
        echo "Hourly size check: /etc/cron.hourly/planetgen-logrotate"
    fi
fi

if ! command -v logrotate >/dev/null 2>&1 && command -v apt-get >/dev/null 2>&1; then
    echo "Installing logrotate (not found)..."
    apt-get install -y logrotate >/dev/null || true
fi
if ! command -v logrotate >/dev/null 2>&1; then
    echo "warning: logrotate isn't installed, so the debug log won't be rotated." >&2
    echo "  Install it with: sudo apt install logrotate" >&2
fi
