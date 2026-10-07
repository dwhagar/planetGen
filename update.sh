#!/usr/bin/env bash
#
# update.sh
#
# Pulls the latest planetGen changes from git, then makes sure everything
# the site needs is present and usable -- without reinstalling anything
# that already is. Each step looks first and changes only what is
# missing, so running this on a server that's already current changes
# nothing (which matters for a scheduled caller, see
# `examples/maintenance/`, that may run it monthly for years):
#
#   1. Pulls (a `git reset --hard` to origin's branch tip, see below) and
#      puts back the executable bit on the repo's shell scripts and
#      `src/html/`'s Python files. A pull rewrites any changed file with
#      whatever mode is tracked in the repo, which has dropped that bit
#      before (see `docs/TODO.md`'s "Deployment bugs found in
#      production").
#   2. Checks every Python library the site needs by importing it with
#      the system Python (the one mod_wsgi and the CLI run under) and
#      comparing its version with setup.py's floor:
#      `scripts/install-python-deps.sh --check`. Only a library that is
#      missing, too old or broken gets installed, the same way install.sh
#      would on this host (apt first and then system-wide pip on an
#      externally managed Python, pip on an ordinary one). planetGen itself is never
#      reinstalled: the web app, the maintenance scripts and the
#      `planetgen` wrapper all run the checkout's code directly.
#   3. The NLTK 'words' corpus: fetched only if it's missing.
#   4. `generate.py check-math`: the math check (TEST.68). A failure only
#      warns (repeated at the end); bulk generation refuses to start until
#      it passes.
#   5. `planetgen.cli.migrate`: brings the database up to the current schema
#      (a no-op when it already is), with a progress bar. When a migration
#      is pending, it first asks (y/N, 30 seconds, default N) whether to
#      delete the galaxy data instead of migrating it; see
#      migrate_or_reset_db in scripts/deploy-common.sh. The update never
#      runs the population pass (OPS.7); run `generate.py population` by
#      hand when wanted.
#   6. Apache's headers, deflate and wsgi modules: enabled only if not
#      already (mod_wsgi installed first if it's missing). On macOS, the
#      gunicorn launchd daemon instead: installed only if it's missing.
#   7. Ownership/permissions for Apache (`examples/apache/set-permissions.sh`),
#      since a pull leaves new and changed files owned by root.
#   8. The tile cache and Generate jobs directories
#      (`examples/apache/create-cache-dir.sh`), and the debug and activity
#      logs (`examples/apache/setup-debug-log.sh`): each one's folder,
#      file, owner and mode, wherever config.json or the environment puts
#      it, and whether Apache's user can write it. A log it can't set up
#      only warns, with the commands that fix it (OPS.5).
#   9. Imports the web app as Apache's user, so anything still unusable
#      fails here instead of as a 500.
#
# Steps 2, 3 and 6 share their code with install.sh (scripts/), so the
# two can't disagree about what a working server needs. `sudo
# ./install.sh` is still there for a full reinstall.
#
# Usage:
#   sudo ./update.sh
#   sudo PYTHON=/usr/bin/python3.12 ./update.sh   (a Python other than python3)
#
# Forces the checkout to match origin's branch tip even if there are
# uncommitted local changes to tracked files (a `git reset --hard` after
# fetching, rather than a `git pull`, which would otherwise refuse or
# conflict) -- a deployment server is meant to always run exactly what's
# on that branch, not whatever got hand-edited on it since. This never
# runs `git clean`, so untracked files -- most importantly `config.json`,
# which is gitignored precisely so it survives this -- are left alone.
#
# Runs on Linux (Apache with mod_wsgi) and macOS (gunicorn under launchd,
# in the venv install-python-deps.sh makes; docs/deployment/macos.md),
# under macOS's bash 3.2 too. update.ps1 is the Windows counterpart and
# follows the same steps.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HTML_DIR="$SCRIPT_DIR/src/html"
DB_DIR="$SCRIPT_DIR/db"
NLTK_DATA_DIR="${PLANETGEN_NLTK_DATA_DIR:-/usr/local/share/nltk_data}"

if [[ $EUID -ne 0 ]]; then
    echo "error: must be run as root, e.g.:" >&2
    echo "  sudo $0" >&2
    exit 1
fi

cd "$SCRIPT_DIR"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "error: $SCRIPT_DIR is not a git checkout -- can't pull an update here." >&2
    exit 1
fi

# shellcheck source=scripts/deploy-common.sh
source "$SCRIPT_DIR/scripts/deploy-common.sh"

PYTHON="${PYTHON:-$(default_python)}"
if [[ -z "$PYTHON" ]]; then
    echo "error: no python3/python found on PATH." >&2
    exit 1
fi

echo "== 1/9: Pulling the latest changes =="

dirty="$(git status --porcelain)"
if [[ -n "$dirty" ]]; then
    echo "warning: uncommitted local changes in $SCRIPT_DIR -- these will be overwritten:" >&2
    echo "$dirty" >&2
fi

branch="$(git rev-parse --abbrev-ref HEAD)"
if [[ "$branch" == "HEAD" ]]; then
    echo "error: repository is in a detached HEAD state -- check out a branch first." >&2
    exit 1
fi

before="$(git rev-parse HEAD)"
git fetch origin "$branch"
# `reset --hard` rather than `pull`: forces every tracked file to match
# origin's branch tip regardless of local commits or uncommitted edits,
# instead of failing on diverged history or a dirty working tree. Only
# touches tracked files -- unlike `git clean`, it never removes untracked
# files, so config.json (gitignored) and anything else local survives.
git reset --hard "origin/$branch"
after="$(git rev-parse HEAD)"

if [[ "$before" == "$after" ]]; then
    echo "Already up to date ($before)."
else
    echo "Updated $before..$after:"
    git log --oneline "$before..$after"
fi
# The scripts below are run directly, so this comes before any of them.
find "$SCRIPT_DIR" -name '*.sh' -exec chmod +x {} +
find "$HTML_DIR" -name '*.py' -exec chmod +x {} +

# Again, now that the pull may have changed it.
# shellcheck source=scripts/deploy-common.sh
source "$SCRIPT_DIR/scripts/deploy-common.sh"

echo
echo "== 2/9: Checking the Python libraries =="
PYTHON="$PYTHON" bash "$SCRIPT_DIR/scripts/install-python-deps.sh" --check
use_site_python

echo
echo "== 3/9: Checking the NLTK 'words' corpus =="
ensure_nltk_words "$NLTK_DATA_DIR"

echo
echo "== 4/9: Checking the generator's math =="
check_math

echo
echo "== 5/9: Migrating the configured MySQL database to the current schema =="
migrate_or_reset_db

echo
echo "== 6/9: Checking Apache's modules (macOS: the gunicorn daemon) =="
APACHE_NEEDS_RESTART=0
if is_macos; then
    ensure_gunicorn_daemon
else
    ensure_apache_modules
fi

echo
echo "== 7/9: Setting directory ownership/permissions for Apache =="
PYTHON="$PYTHON" "$SCRIPT_DIR/examples/apache/set-permissions.sh" "$HTML_DIR" "$DB_DIR"

echo
echo "== 8/9: Checking the cache, jobs and debug log locations =="
PYTHON="$PYTHON" "$SCRIPT_DIR/examples/apache/create-cache-dir.sh"
# The logs never stop an update (OPS.5): setup-debug-log.sh warns, with
# the commands that fix it, for anything it can't set up itself.
PYTHON="$PYTHON" "$SCRIPT_DIR/examples/apache/setup-debug-log.sh" \
    || echo "warning: the log setup didn't finish (see above); the update carries on." >&2

echo
echo "== 9/9: Checking that the web app imports =="
check_app_imports

echo
echo "Checking the Redis server (config.json's redis.url):"
ensure_redis

echo
if [[ "$MATH_CHECK_FAILED" == 1 ]]; then
    echo "warning: the math check failed (step 4): bulk generation refuses to start until it passes." >&2
fi
echo "The update doesn't run the population pass (species, civilizations, territories); when wanted:"
echo "  $PYTHON $SCRIPT_DIR/generate.py population"
if is_macos; then
    if [[ "$before" != "$after" ]]; then
        echo "Done. Reload gunicorn so the site runs the new code (a running Generate job is left alone):"
        echo "  sudo launchctl kill SIGHUP system/$GUNICORN_LABEL"
    else
        echo "Done. Nothing new was pulled."
    fi
elif (( APACHE_NEEDS_RESTART )); then
    echo "Done. An Apache module was just enabled: restart Apache to load it:"
    echo "  sudo systemctl restart apache2"
elif [[ "$before" != "$after" ]]; then
    echo "Done. Reload Apache so the site runs the new code:"
    echo "  sudo systemctl reload apache2"
else
    echo "Done. Nothing new was pulled."
fi
