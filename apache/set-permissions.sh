#!/usr/bin/env bash
#
# apache/set-permissions.sh
#
# Detects the user/group Apache2 actually runs as and sets ownership and
# permissions on the deployed planetGen web directory (and the database
# directory it reads from) accordingly. Linux-only (reads /etc/apache2 and
# uses `ps`), deliberately bash rather than Python -- this is a one-shot
# root-privileged deployment step, not part of the portable application.
#
# Usage:
#   sudo apache/set-permissions.sh [html-dir] [db-dir]
#
# Defaults match the layout documented in html/README.md and
# apache/planetgen.conf.example:
#   html-dir defaults to /var/lib/planetGen/html
#   db-dir   defaults to /var/lib/planetGen/db
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
#     Apache config level (see apache/planetgen.conf.example), not by
#     filesystem permissions, since Apache's own worker still needs to
#     read those modules to import them.

set -euo pipefail

HTML_DIR="${1:-/var/lib/planetGen/html}"
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

detect_apache_group() {
    # Debian/Ubuntu's apache2 package always defines APACHE_RUN_USER/
    # APACHE_RUN_GROUP in /etc/apache2/envvars (default www-data) -- this
    # is the authoritative source, since it's what apachectl itself reads
    # to decide who the worker processes run as.
    if [[ -r /etc/apache2/envvars ]]; then
        local user group
        # Run in a subshell so envvars' own `export`s and any unrelated
        # variables it sets don't leak into this script.
        user="$(source /etc/apache2/envvars >/dev/null 2>&1; echo "$APACHE_RUN_USER")"
        group="$(source /etc/apache2/envvars >/dev/null 2>&1; echo "$APACHE_RUN_GROUP")"
        if [[ -n "$group" ]]; then
            echo "$user" "$group"
            return 0
        fi
    fi

    # Fallback: ask a running apache2 for its worker identity directly --
    # the master process runs as root, so pick the first non-root owner
    # among its children. Covers non-Debian layouts (no envvars file)
    # where apache2 is already running.
    local line
    line="$(ps -eo user:32,group:32,comm | awk '$3 == "apache2" && $1 != "root" {print $1, $2; exit}')"
    if [[ -n "$line" ]]; then
        echo "$line"
        return 0
    fi

    # Last resort: neither the config file nor a running process was
    # found -- most likely because apache2 hasn't been started/enabled
    # yet (e.g. running as part of install.sh, before the vhost exists).
    # www-data:www-data is the standard Debian/Ubuntu apache2 identity,
    # so assume it rather than hard-failing the whole install.
    echo "warning: could not detect Apache's user/group (no readable" >&2
    echo "  /etc/apache2/envvars and no running apache2 process found)." >&2
    echo "  Assuming the Debian/Ubuntu default: www-data:www-data" >&2
    echo "www-data www-data"
}

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
