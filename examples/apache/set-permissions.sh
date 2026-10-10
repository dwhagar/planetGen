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
#   - html-dir is code: recursively owned by root:<Apache group>, so
#     Apache's worker can read (and import) it but never change it. A web
#     process that could rewrite the code tree could plant Python that
#     install.sh/update.sh later run as root; nothing under src/html is
#     written at runtime (the tile cache and Generate jobs live in their
#     own directories, see create-cache-dir.sh).
#   - db-dir is runtime data (the legacy pre-MySQL `.db` files): owned by
#     Apache's detected user:group (ownership as well as group, so this
#     doesn't depend on the deploying user's own group memberships).
#   - Directories: 750 (owner rwx, group r-x, others none) so Apache can
#     traverse and list them but other local users can't.
#   - Regular files: 640 (owner rw, group r, others none).
#   - Every `*.py` file anywhere under html-dir, at any subdirectory
#     depth (`html/wsgi.py`, ...): 750 (owner+group
#     read, so mod_wsgi's daemon can import them; mod_wsgi only needs
#     read access, the execute bit is a leftover from the old CGI pages
#     and harmless). Reported with a count at the end so a wrong
#     `html-dir` path is obvious rather than silently matching zero files.
#   - config.json (at the repo root, two levels above html-dir), when it
#     exists: root:<Apache group>, mode 640 -- it holds the database
#     password and the session secret_key, so only root and Apache may
#     read it.

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

if ! read -r APACHE_USER APACHE_GROUP < <(detect_apache_group) || [[ -z "${APACHE_GROUP:-}" ]]; then
    echo "error: couldn't work out Apache's user and group (see above); nothing was changed." >&2
    exit 1
fi

echo "Detected Apache identity: user=$APACHE_USER group=$APACHE_GROUP"
echo "Applying permissions to:"
echo "  html: $HTML_DIR"
[[ -d "$DB_DIR" ]] && echo "  db:   $DB_DIR" || echo "  db:   $DB_DIR (not found -- skipping)"

apply_permissions() {
    local dir="$1" owner="$2"
    chown -R "$owner:$APACHE_GROUP" "$dir"
    find "$dir" -type d -exec chmod 750 {} +
    find "$dir" -type f -exec chmod 640 {} +

    # No -maxdepth here on purpose: this must reach *.py files at any
    # subdirectory depth, not just directly
    # inside $dir -- a previous version of this script effectively only
    # fixed the top-level scripts, which is exactly what was reported
    # broken.
    find "$dir" -type f -name '*.py' -exec chmod 750 {} +
    local py_count
    py_count="$(find "$dir" -type f -name '*.py' | wc -l | tr -d ' ')"  # BSD wc pads
    echo "  $dir: owned by $owner:$APACHE_GROUP, made $py_count .py file(s) executable (owner+group)"
}

# The code tree: root-owned, read-only for Apache.
apply_permissions "$HTML_DIR" root
# Runtime data: Apache-owned.
if [[ -d "$DB_DIR" ]]; then
    apply_permissions "$DB_DIR" "$APACHE_USER"
fi

# config.json holds the DB password and secret_key.
CONFIG_JSON="$(cd "$HTML_DIR/../.." && pwd)/config.json"
if [[ -f "$CONFIG_JSON" ]]; then
    chown root:"$APACHE_GROUP" "$CONFIG_JSON"
    chmod 640 "$CONFIG_JSON"
    echo "  $CONFIG_JSON: owned by root:$APACHE_GROUP, mode 640"
fi

echo "Done."
