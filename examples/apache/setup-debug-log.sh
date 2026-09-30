#!/usr/bin/env bash
#
# examples/apache/setup-debug-log.sh
#
# Prepares planetGen's debug log (see src/stellarObjects/log.py) and its
# log rotation. Safe to run again. Called by install.sh, and by update.sh
# when there's nothing new to install.
#
#   1. When debug is on (config.json's "debug": true, or PLANETGEN_DEBUG in
#      this shell), creates the log file -- "log_file" in config.json,
#      default /var/log/planetgen.log. An existing file is kept either way.
#   2. Makes the file owned by Apache's worker user and group, mode 0660:
#      the web interface writes to it as that user, root (a systemd timer,
#      sudo) can always write to it, and a login user who runs the
#      generator from a shell must be in the Apache group to append to it
#      (`sudo usermod -aG <apache group> <user>`, then log in again; on
#      macOS `sudo dseditgroup -o edit -a <user> -t user _www`).
#      Other users can neither read nor write it -- a world-writable log
#      would let any local user forge or flood its entries.
#   3. Installs /etc/logrotate.d/planetgen, so a debug log someone forgot to
#      turn off can't fill the disk: rotated daily, or as soon as it passes
#      100 MB, keeping 7 compressed copies. Where /etc/cron.hourly exists it
#      also adds a job that checks the size every hour, since a galaxy run
#      with debug on can write gigabytes between daily logrotate runs.
#      On macOS, which has newsyslog instead of logrotate, installs
#      /etc/newsyslog.d/planetgen.conf: rotated past 100 MB, 7 compressed
#      copies kept (newsyslog runs hourly by itself).
#
# Runs on Linux and macOS. Windows has no rotation for it; install.ps1
# only creates the logs folder.
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
# itself reads them (stellarObjects/appconfig.py). appconfig is loaded
# straight from its file, not through the stellarObjects package, so this
# works before the package's own dependencies are installed. This runs as
# root, so -I keeps the current directory, user site-packages and PYTHON*
# variables off sys.path (appconfig.py itself lives outside src/html, in
# the root-owned part of the checkout).
read -r DEBUG_ON LOG_FILE < <(cd / && "$PYTHON" -I - "$REPO_DIR" <<'PY'
import importlib.util
import os
import sys

path = os.path.join(sys.argv[1], "src", "stellarObjects", "appconfig.py")
spec = importlib.util.spec_from_file_location("planetgen_appconfig", path)
appconfig = importlib.util.module_from_spec(spec)
spec.loader.exec_module(appconfig)
config = appconfig.load_config()
print(1 if appconfig.debug_enabled(config) else 0, appconfig.log_file_path(config))
PY
)

# shellcheck source=examples/apache/apache-identity.sh
source "$APACHE_DIR/apache-identity.sh"
read -r APACHE_USER APACHE_GROUP < <(detect_apache_group)

if [[ "$DEBUG_ON" == "1" && ! -e "$LOG_FILE" ]]; then
    mkdir -p "$(dirname "$LOG_FILE")"
    touch "$LOG_FILE"
    echo "Debug is on: created $LOG_FILE"
fi

if [[ -e "$LOG_FILE" ]]; then
    chown "$APACHE_USER:$APACHE_GROUP" "$LOG_FILE"
    chmod 0660 "$LOG_FILE"
    echo "Debug log: $LOG_FILE (owned by $APACHE_USER:$APACHE_GROUP, mode 0660;"
    echo "  add CLI users who run the generator to the $APACHE_GROUP group so they can append)"
elif [[ "$DEBUG_ON" != "1" ]]; then
    echo "Debug is off: no debug log to create (set \"debug\": true in config.json and re-run to create $LOG_FILE)."
fi

# Rotation is set up whether or not debug is on right now, so turning it on
# later is covered.
if [[ "$(uname -s)" == Darwin ]]; then
    # newsyslog: owner:group, mode, copies kept, size in KB, no time
    # rotation, Z = gzip. It creates the new file with that owner and mode.
    mkdir -p /etc/newsyslog.d
    cat > /etc/newsyslog.d/planetgen.conf <<EOF
# planetGen debug log -- written by examples/apache/setup-debug-log.sh
# (install.sh/update.sh). Edits here are overwritten on the next run.
# logfilename            [owner:group]              mode count size(KB) when  flags
$LOG_FILE  $APACHE_USER:$APACHE_GROUP  660  7  102400  *  Z
EOF
    chmod 0644 /etc/newsyslog.d/planetgen.conf
    echo "Log rotation: /etc/newsyslog.d/planetgen.conf (past 100 MB, 7 kept)"
    exit 0
fi

# `missingok` makes it a no-op while there's no file.
# `create` recreates the file with the same owner and mode after each
# rotation; the program reopens it on its own (WatchedFileHandler).
cat > /etc/logrotate.d/planetgen <<EOF
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
chmod 0644 /etc/logrotate.d/planetgen
echo "Log rotation: /etc/logrotate.d/planetgen (daily or past 100 MB, 7 kept)"

if [[ -d /etc/cron.hourly ]]; then
    cat > /etc/cron.hourly/planetgen-logrotate <<'EOF'
#!/bin/sh
# planetGen: rotate the debug log as soon as it passes its maxsize instead of
# waiting for the daily logrotate run. Written by
# examples/apache/setup-debug-log.sh.
command -v logrotate >/dev/null 2>&1 || exit 0
exec logrotate /etc/logrotate.d/planetgen
EOF
    chmod 0755 /etc/cron.hourly/planetgen-logrotate
    echo "Hourly size check: /etc/cron.hourly/planetgen-logrotate"
fi

if ! command -v logrotate >/dev/null 2>&1 && command -v apt-get >/dev/null 2>&1; then
    echo "Installing logrotate (not found)..."
    apt-get install -y logrotate >/dev/null || true
fi
if ! command -v logrotate >/dev/null 2>&1; then
    echo "warning: logrotate isn't installed, so the debug log won't be rotated." >&2
    echo "  Install it with: sudo apt install logrotate" >&2
fi
