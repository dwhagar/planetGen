#!/usr/bin/env bash
#
# examples/apache/create-cache-dir.sh
#
# Creates the web interface's on-disk Galaxy Map tile cache (see
# src/html/lib/tilecache.py), and the admin Generate page's jobs directory
# (src/html/web/jobs.py), and gives them to Apache's worker user, so the
# web interface can write to them. Safe to run again: an existing directory is
# only re-owned. Called by install.sh, and by update.sh when there's
# nothing new to install.
#
# Usage:
#   sudo examples/apache/create-cache-dir.sh [cache-dir]
#
# Without an argument the directory is the one the web interface will use:
# `PLANETGEN_TILE_CACHE_DIR` if this shell has it, else config.json's
# `tile_cache.dir`, else /var/cache/planetgen/tiles. Nothing is created
# when `tile_cache.max_mb` is 0 (the disk cache is off). A directory set
# only through a `SetEnv PLANETGEN_TILE_CACHE_DIR` in the Apache vhost
# isn't visible here, so pass it as the argument.

set -euo pipefail

APACHE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$APACHE_DIR/../.." && pwd)"

if [[ $EUID -ne 0 ]]; then
    echo "error: must be run as root (needs chown), e.g.:" >&2
    echo "  sudo $0 $*" >&2
    exit 1
fi

CACHE_DIR="${1:-}"
if [[ -z "$CACHE_DIR" ]]; then
    PYTHON="$(command -v python3 || command -v python || true)"
    if [[ -z "$PYTHON" ]]; then
        echo "error: no python3/python found on PATH." >&2
        exit 1
    fi
    CACHE_DIR="$("$PYTHON" - "$REPO_DIR" <<'PY'
import os
import sys

repo = sys.argv[1]
sys.path[:0] = [os.path.join(repo, "src", "html", "lib"), os.path.join(repo, "src")]
import tilecache

print(tilecache.configured_cache_dir() or "")
PY
)"
fi

# shellcheck source=examples/apache/apache-identity.sh
source "$APACHE_DIR/apache-identity.sh"
read -r APACHE_USER APACHE_GROUP < <(detect_apache_group)

# Parents (e.g. /var/cache/planetgen) stay root-owned and world-traversable;
# only the cache directory itself belongs to Apache.
if [[ -z "$CACHE_DIR" ]]; then
    echo "Tile cache is off (tile_cache.max_mb is 0) -- nothing to create."
else
    mkdir -p "$CACHE_DIR"
    chown -R "$APACHE_USER:$APACHE_GROUP" "$CACHE_DIR"
    chmod 750 "$CACHE_DIR"
    echo "Tile cache: $CACHE_DIR (owned by $APACHE_USER:$APACHE_GROUP)"
fi

# The admin Generate page's background jobs (src/html/web/jobs.py):
# PLANETGEN_JOBS_DIR, else config.json's jobs.dir, else
# /var/lib/planetgen/jobs. Only when no tile cache directory was given as
# the argument, since that argument names the tile cache alone.
if [[ -z "${1:-}" ]]; then
    PYTHON="$(command -v python3 || command -v python || true)"
    JOBS_DIR="$("$PYTHON" - "$REPO_DIR" <<'PY'
import os
import sys

repo = sys.argv[1]
sys.path[:0] = [os.path.join(repo, "src", "html"), os.path.join(repo, "src")]
from stellarObjects.appconfig import load_config

print(os.environ.get("PLANETGEN_JOBS_DIR") or load_config()["jobs"].get("dir") or "/var/lib/planetgen/jobs")
PY
)"
    mkdir -p "$JOBS_DIR"
    chown "$APACHE_USER:$APACHE_GROUP" "$JOBS_DIR"
    chmod 750 "$JOBS_DIR"
    echo "Generate jobs: $JOBS_DIR (owned by $APACHE_USER:$APACHE_GROUP)"
fi
