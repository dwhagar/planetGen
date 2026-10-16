#!/usr/bin/env bash
#
# examples/apache/set-permissions.sh
#
# Detects the user/group Apache2 actually runs as (on macOS, _www, which
# gunicorn runs as) and sets ownership and permissions on the deployed
# planetGen web directory accordingly. Runs on Linux and macOS (only
# chown, chmod and find flags both have).
# Deliberately bash rather than Python -- this is a one-shot
# root-privileged deployment step, not part of the portable application.
#
# Usage:
#   sudo examples/apache/set-permissions.sh [install-dir] [state-dir]
#
# Defaults match the layout documented in docs/html-interface.md and
# examples/apache/planetgen.conf.example:
#   install-dir defaults to /var/lib/planetGen (the checkout: src/, db/,
#               scripts/, config.json, .git, ... -- install.sh and
#               update.sh pass the checkout they run from)
#   state-dir   defaults to <install-dir>/state, the files the web app
#               writes at run time (the maintenance-mode flag, ...)
#
# What it does:
#   - The whole install tree, and the state directory (made when it is
#     missing), belong to Apache's detected user:group (www-data:www-data
#     on Debian/Ubuntu), so the web app can write what it needs to and
#     the files look the same whether install.sh or update.sh touched
#     them last. install.sh and update.sh both call this script, which is
#     what keeps a fresh install and an updated one identical
#     (test_deploy_scripts.py checks it).
#   - Directories: 750 (owner rwx, group r-x, others none).
#   - Regular files: 640 (owner rw, group r, others none); a file that
#     already has its owner execute bit, and every `*.py` and `*.sh`
#     file at any depth (no -maxdepth on purpose: an earlier version
#     only fixed top-level scripts, which was reported broken), gets 750,
#     so git sees no mode change and the scripts still run.
#   - Symbolic links are left alone (chown -h), never followed.
#   - What stays root's: the virtual environment (/usr/local/planetgen/
#     venv) and the NLTK corpus. install.sh and update.sh run their
#     Python as root; the web user must not be able to rewrite it.
#
# Apache itself is kept out of everything but the WSGI app and static/
# by the <Directory> rules in planetgen.conf.example.

set -euo pipefail

INSTALL_DIR="${1:-/var/lib/planetGen}"
STATE_DIR="${2:-$INSTALL_DIR/state}"

if [[ $EUID -ne 0 ]]; then
    echo "error: must be run as root (needs chown/chgrp), e.g.:" >&2
    echo "  sudo $0 $*" >&2
    exit 1
fi

if [[ ! -d "$INSTALL_DIR" ]]; then
    echo "error: install directory not found: $INSTALL_DIR" >&2
    exit 1
fi

if [[ ! -d "$INSTALL_DIR/src/html" ]]; then
    echo "error: $INSTALL_DIR doesn't look like a planetGen checkout (no src/html); nothing was changed." >&2
    exit 1
fi

# shellcheck source=examples/apache/apache-identity.sh
source "$(dirname "${BASH_SOURCE[0]}")/apache-identity.sh"

if ! read -r APACHE_USER APACHE_GROUP < <(detect_apache_group) || [[ -z "${APACHE_GROUP:-}" ]]; then
    echo "error: couldn't work out Apache's user and group (see above); nothing was changed." >&2
    exit 1
fi

echo "Detected Apache identity: user=$APACHE_USER group=$APACHE_GROUP"
echo "Applying permissions to:"
echo "  install: $INSTALL_DIR"
echo "  state:   $STATE_DIR"

apply_permissions() {
    local dir="$1"
    chown -R -h "$APACHE_USER:$APACHE_GROUP" "$dir"
    find "$dir" -type d -exec chmod 750 {} +
    # A file that is already owner-executable keeps that (750); every
    # other file is 640.
    find "$dir" -type f ! -perm -u+x -exec chmod 640 {} +
    find "$dir" -type f -perm -u+x -exec chmod 750 {} +
    find "$dir" -type f \( -name '*.py' -o -name '*.sh' \) -exec chmod 750 {} +
    local py_count
    py_count="$(find "$dir" -type f -name '*.py' | wc -l | tr -d ' ')"  # BSD wc pads
    echo "  $dir: owned by $APACHE_USER:$APACHE_GROUP, $py_count .py file(s) executable (owner+group)"
}

# The state directory first, so the walk below covers it when it is inside
# the install tree (the default).
mkdir -p "$STATE_DIR"
apply_permissions "$INSTALL_DIR"
if [[ "$(cd "$STATE_DIR" && pwd -P)/" != "$(cd "$INSTALL_DIR" && pwd -P)/"* ]]; then
    apply_permissions "$STATE_DIR"
fi

echo "Done."
