#!/usr/bin/env bash
#
# examples/apache/set-permissions.sh
#
# Detects the user/group Apache2 actually runs as and sets ownership and
# permissions on the deployed planetGen web directory accordingly.
# Linux-only (reads /etc/apache2 and uses `ps`), deliberately bash rather
# than Python -- this is a one-shot root-privileged deployment step, not
# part of the portable application.
#
# Usage:
#   sudo examples/apache/set-permissions.sh [html-dir] [db-dir]
#
# Defaults match the layout documented in docs/html-interface.md and
# examples/apache/planetgen.conf.example:
#   html-dir defaults to /var/lib/planetGen/src/html
#   db-dir   defaults to /var/lib/planetGen/db, but the web interface has
#            read from a MySQL server (TODO.md Phase 5), not a local
#            file, since the SQLite-to-MySQL port -- this argument is
#            kept only for a deployment that still has an old `db/`
#            directory of pre-port `.db` files lying around (harmless to
#            chown, and skipped entirely if the directory doesn't exist).
#
# What it does:
#   - Recursively `chown`s both directories to Apache's detected
#     user:group (not just `chgrp` -- ownership as well as group
#     membership, so this doesn't depend on the deploying user's own
#     group memberships lining up).
#   - Directories: 750 (owner rwx, group r-x, others none) so Apache can
#     traverse and list them but other local users can't.
#   - Regular files: 640 (owner rw, group r, others none).
#   - Every `*.py` file anywhere under html-dir, at any subdirectory
#     depth (`html/*.py`, `html/lib/*.py`, ...): 750 (adds execute, since
#     Apache must be able to execute the CGI scripts, and `lib/`'s own
#     modules need at least read access to be importable). Reported with
#     a count at the end so a wrong `html-dir` path is obvious rather
#     than silently matching zero files.
#   - html/lib is included in the general file/directory pass like any
#     other subdirectory -- direct web access to it is denied at the
#     Apache config level (see examples/apache/planetgen.conf.example), not by
#     filesystem permissions, since Apache's own worker still needs to
#     read those modules to import them.

set -euo pipefail

HTML_DIR="${1:-/var/lib/planetGen/src/html}"
DB_DIR="${2:-/var/lib/planetGen/db}"

if [[ $EUID -ne 0 ]]; then
    echo "error: must be run as root (needs chown/chgrp), e.g.:" >&2
    echo "  sudo $0 $*" >&2
    exit 1
fi

if [[ ! -d "$HTML_DIR" ]]; then
    echo "error: html directory not found: $HTML_DIR" >&2
    exit 1
fi

# shellcheck source=examples/apache/apache-identity.sh
source "$(dirname "${BASH_SOURCE[0]}")/apache-identity.sh"

read -r APACHE_USER APACHE_GROUP < <(detect_apache_group)

echo "Detected Apache identity: user=$APACHE_USER group=$APACHE_GROUP"
echo "Applying permissions to:"
echo "  html: $HTML_DIR"
[[ -d "$DB_DIR" ]] && echo "  db:   $DB_DIR" || echo "  db:   $DB_DIR (not found -- skipping)"

apply_permissions() {
    local dir="$1"
    chown -R "$APACHE_USER:$APACHE_GROUP" "$dir"
    find "$dir" -type d -exec chmod 750 {} +
    find "$dir" -type f -exec chmod 640 {} +

    # No -maxdepth here on purpose: this must reach *.py files at any
    # subdirectory depth (html/lib/*.py included), not just directly
    # inside $dir -- a previous version of this script effectively only
    # fixed the top-level scripts, which is exactly what was reported
    # broken.
    find "$dir" -type f -name '*.py' -exec chmod 750 {} +
    local py_count
    py_count="$(find "$dir" -type f -name '*.py' | wc -l)"
    echo "  $dir: made $py_count .py file(s) executable (owner+group)"
}

apply_permissions "$HTML_DIR"
if [[ -d "$DB_DIR" ]]; then
    apply_permissions "$DB_DIR"
fi

echo "Done."
