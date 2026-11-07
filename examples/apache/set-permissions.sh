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
# Defaults match the layout documented in docs/reference/html-interface.md and
# examples/apache/planetgen.conf.example:
#   install-dir defaults to /var/lib/planetGen (the checkout: src/, db/,
#               scripts/, config.json, .git, ... -- install.sh and
#               update.sh pass the checkout they run from)
#   state-dir   defaults to <install-dir>/state, the files the web app
#               writes at run time (the maintenance-mode flag, ...)
#
# What it does (SEC.37: no file the web user can write is ever run as root; install.sh and update.sh run git, python and the
# scripts in this tree as root):
#   - The code tree (everything in the install directory except the web user's own folders below: src/, scripts/,
#     update.sh, .git, the symbolic link config.json, ...) belongs to root, group = Apache's group. Directories 750 (the group
#     may read and enter, nobody may write), regular files 640, and every `*.py` and `*.sh` file at any depth (plus any file
#     already owner-executable) 750, so git sees no mode change and the scripts still run. The web user reads its code through
#     the group and cannot change a byte of it.
#   - The web user's own folders belong to Apache's detected user:group (www-data:www-data on Debian/Ubuntu): the state
#     directory (made when it is missing; 770, so an operator in the group can run `python3 -m planetgen.cli.maintenance`)
#     and the jobs directory the Generate page writes (750). Both are skipped by the code walk. Files in them are 640.
#   - config.json lives in the state directory (www-data, 640), where the settings page can replace it atomically, and
#     <install-dir>/config.json is a root-owned symbolic link to it. A real config.json found in the install directory is moved
#     there once (the old checkout-owned layout); if both exist nothing is moved and a warning says so.
#   - The backups folder inside the state directory (planetgen.backup) is 700 and its archives 600: the web user's own.
#   - Symbolic links are left alone (chown -h), never followed.
#   - What stays root's, as before: the virtual environment (/usr/local/planetgen/venv) and the NLTK corpus.
#   - The tile cache, the debug log and the other log folders are outside the tree and are set up by create-cache-dir.sh and
#     setup-debug-log.sh.
#
# Apache itself is kept out of everything but the WSGI app and static/
# by the <Directory> rules in planetgen.conf.example.

set -euo pipefail

# The lines this script prints are log records (OPS.98): a time, a severity and a component on every one.
LOG_COMPONENT="${LOG_COMPONENT:-permissions}"
export LOG_COMPONENT
# shellcheck source=scripts/logrecord.sh
source "$(dirname "${BASH_SOURCE[0]}")/../../scripts/logrecord.sh"

INSTALL_DIR="${1:-/var/lib/planetGen}"
STATE_DIR="${2:-$INSTALL_DIR/state}"

if [[ $EUID -ne 0 ]]; then
    log_error "must be run as root (needs chown/chgrp), e.g.:"
    log_error "  sudo $0 $*"
    exit 1
fi

if [[ ! -d "$INSTALL_DIR" ]]; then
    log_error "install directory not found: $INSTALL_DIR"
    exit 1
fi

if [[ ! -d "$INSTALL_DIR/src/html" ]]; then
    log_error "$INSTALL_DIR doesn't look like a planetGen checkout (no src/html); nothing was changed."
    if [[ -f "$INSTALL_DIR/wsgi.py" ]]; then
        log_error "  This takes the checkout (/var/lib/planetGen), not src/html. An update.sh that was already"
        log_error "  running when it pulled this change still passes src/html: run ./update.sh again."
    fi
    exit 1
fi

# shellcheck source=examples/apache/apache-identity.sh
source "$(dirname "${BASH_SOURCE[0]}")/apache-identity.sh"

if ! read -r APACHE_USER APACHE_GROUP < <(detect_apache_group) || [[ -z "${APACHE_GROUP:-}" ]]; then
    log_error "couldn't work out Apache's user and group (see above); nothing was changed."
    exit 1
fi

log_info "Detected Apache identity: user=$APACHE_USER group=$APACHE_GROUP"
log_info "Applying permissions to:"
log_info "  install: $INSTALL_DIR"
log_info "  state:   $STATE_DIR"

# The web user's own folders: state and the Generate page's jobs directory (deploy-paths.py gives it, read under -I).
apply_web_permissions() {
    local dir="$1" mode="$2"
    chown -R -h "$APACHE_USER:$APACHE_GROUP" "$dir"
    find "$dir" -type d -exec chmod "$mode" {} +
    find "$dir" -type f -exec chmod 640 {} +
    log_info "  $dir: owned by $APACHE_USER:$APACHE_GROUP, directories $mode, files 640"
}

# Everything else: root's, readable by the web user through its group. $WEB_DIRS (canonical paths inside the install tree)
# are skipped, whole. Bash 3.2: no empty-array expansion without the guard.
apply_code_permissions() {
    local dir="$1" prune=() path
    for path in ${WEB_DIRS[@]+"${WEB_DIRS[@]}"}; do
        if [[ ${#prune[@]} -gt 0 ]]; then prune+=(-o); fi
        prune+=(-path "$path")
    done
    local skip=()
    if [[ ${#prune[@]} -gt 0 ]]; then skip=(\( ${prune[@]+"${prune[@]}"} \) -prune -o); fi
    find "$dir" ${skip[@]+"${skip[@]}"} -exec chown -h "root:$APACHE_GROUP" {} +
    find "$dir" ${skip[@]+"${skip[@]}"} -type d -exec chmod 750 {} +
    # A file that is already owner-executable keeps that (750); every other file is 640.
    find "$dir" ${skip[@]+"${skip[@]}"} -type f ! -perm -u+x -exec chmod 640 {} +
    find "$dir" ${skip[@]+"${skip[@]}"} -type f -perm -u+x -exec chmod 750 {} +
    find "$dir" ${skip[@]+"${skip[@]}"} -type f \( -name '*.py' -o -name '*.sh' \) -exec chmod 750 {} +
    local py_count
    py_count="$(find "$dir" ${skip[@]+"${skip[@]}"} -type f -name '*.py' | wc -l | tr -d ' ')"  # BSD wc pads
    log_info "  $dir: owned by root:$APACHE_GROUP (the web user cannot write it), $py_count .py file(s) executable"
}

canonical() {
    (cd "$1" 2>/dev/null && pwd -P)
}

mkdir -p "$STATE_DIR"
INSTALL_P="$(canonical "$INSTALL_DIR")"
STATE_P="$(canonical "$STATE_DIR")"
WEB_DIRS=()
WEB_OUTSIDE=()
# Jobs directory (config.json or PLANETGEN_JOBS_DIR): the web user's, wherever it is. Not made here (create-cache-dir.sh does).
JOBS_P=""
if paths="$(cd / && "${PYTHON:-python3}" -I "$(dirname "${BASH_SOURCE[0]}")/deploy-paths.py" "$INSTALL_DIR" 2>/dev/null)"; then
    JOBS_DIR="$(printf '%s\n' "$paths" | sed -n 2p)"
    # SEC.44: config.json is the web user's file, so the folder it names is only given to the web user (chown -R) when
    # rootpaths.py allows it; otherwise it is left as it is, and inside the checkout it stays part of root's code tree.
    if [[ -n "$JOBS_DIR" ]] && ! why="$(cd / && "${PYTHON:-python3}" -I "$(dirname "${BASH_SOURCE[0]}")/deploy-paths.py" --check jobs "$JOBS_DIR" "$INSTALL_DIR" "$APACHE_USER" "$APACHE_GROUP" 2>&1)"; then
        log_warning "${why#refused: }"
        log_warning "  The jobs folder was left as it is. Name a folder under /var/lib, /var/cache, /var/log or /srv in config.json,"
        log_warning "  or set PLANETGEN_JOBS_DIR in root's shell, then run sudo ./update.sh again."
        JOBS_DIR=""
    fi
    [[ -n "$JOBS_DIR" && -d "$JOBS_DIR" ]] && JOBS_P="$(canonical "$JOBS_DIR")"
else
    log_warning "could not read the jobs directory from deploy-paths.py; the jobs folder is not touched (create-cache-dir.sh sets it up)."
fi
for dir in "$STATE_P" "$JOBS_P"; do
    [[ -n "$dir" ]] || continue
    if [[ "$dir/" == "$INSTALL_P/"* && "$dir" != "$INSTALL_P" ]]; then WEB_DIRS+=("$dir"); else WEB_OUTSIDE+=("$dir"); fi
done

# config.json belongs in the state directory, where the settings page can replace it; the checkout's own is a root-owned link.
CONFIG="$INSTALL_DIR/config.json"
WEB_CONFIG="$STATE_P/config.json"
if [[ -f "$CONFIG" && ! -L "$CONFIG" ]]; then
    if [[ -e "$WEB_CONFIG" ]]; then
        log_warning "both $CONFIG and $WEB_CONFIG exist; nothing was moved. Keep the one you want in $WEB_CONFIG and delete the other."
    else
        mv "$CONFIG" "$WEB_CONFIG"
        log_info "  moved config.json into $STATE_P (the settings page writes it there)"
    fi
fi
if [[ ! -e "$CONFIG" && ! -L "$CONFIG" && -e "$WEB_CONFIG" ]]; then
    ln -s "$WEB_CONFIG" "$CONFIG"
    log_info "  config.json in the checkout is now a link to $WEB_CONFIG"
fi

apply_code_permissions "$INSTALL_P"
for dir in ${WEB_DIRS[@]+"${WEB_DIRS[@]}"} ${WEB_OUTSIDE[@]+"${WEB_OUTSIDE[@]}"}; do
    if [[ "$dir" == "$STATE_P" ]]; then apply_web_permissions "$dir" 770; else apply_web_permissions "$dir" 750; fi
done
# The state directory: the one folder the group may also write (an operator in Apache's group runs
# `python3 -m planetgen.cli.maintenance`); never world-writable. create-cache-dir.sh sets the same mode.
chmod 770 "$STATE_P"
# The backup archives (planetgen.backup) hold a copy of config.json, and later the control database's password hashes:
# the web user's own, nobody else's. The folder is made by the first backup; this keeps it that way.
if [[ -d "$STATE_DIR/backups" ]]; then
    chmod 700 "$STATE_DIR/backups"
    find "$STATE_DIR/backups" -type f -exec chmod 600 {} +
fi

log_info "Done."
