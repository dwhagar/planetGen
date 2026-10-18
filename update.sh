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
#      the venv's python (/usr/local/planetgen/venv, OPS.45: the one
#      mod_wsgi, gunicorn and the CLI run under) and comparing its version
#      with setup.py's floor: `scripts/install-python-deps.sh --check`. The
#      venv is made, and the locked libraries installed, only when it is
#      missing or a library is missing, too old or broken. planetGen itself
#      is never reinstalled: the web app, the maintenance scripts and the
#      `planetgen` wrapper all run the checkout's code directly.
#   3. The NLTK 'words' corpus: fetched only if it's missing.
#   4. `planetgen check-math`: the math check (TEST.68). A failure only
#      warns (repeated at the end); bulk generation refuses to start until
#      it passes.
#   5. `planetgen.cli.migrate`: brings the database up to the current schema
#      (a no-op when it already is), with a progress bar. When a migration
#      is pending, it first asks (y/N, 30 seconds, default N) whether to
#      delete the galaxy data instead of migrating it; see
#      migrate_or_reset_db in scripts/deploy-common.sh. The update never
#      runs the population pass (OPS.7); run `planetgen population` by
#      hand when wanted.
#      Then it records the version key (OPS.13): one row per galaxy in the
#      control database, the last 10 kept; `planetgen versions` lists them.
#   6. Apache's headers, deflate and wsgi modules: enabled only if not
#      already (mod_wsgi installed first if it's missing). On macOS, the
#      gunicorn launchd daemon instead: installed only if it's missing.
#   7. Ownership/permissions (`examples/apache/set-permissions.sh`): the
#      whole checkout and its state directory belong to Apache's user
#      (www-data:www-data), since a pull leaves new and changed files
#      owned by root. The virtual environment stays root's.
#   8. The tile cache and Generate jobs directories
#      (`examples/apache/create-cache-dir.sh`), and the debug and activity
#      logs (`examples/apache/setup-debug-log.sh`): each one's folder,
#      file, owner and mode, wherever config.json or the environment puts
#      it, and whether Apache's user can write it. A log it can't set up
#      only warns, with the commands that fix it (OPS.5).
#   9. Imports the web app as Apache's user, so anything still unusable
#      fails here instead of as a 500.
#
# At the end, on Linux, it reloads Apache itself when Apache is running
# (restarts it when step 6 just enabled a module) and says so; if Apache
# isn't running it prints the command instead (OPS.8).
#
# After the steps it starts building the site's caches in the background
# (`planetgen warm`, MAP.134, PERF.104 and PERF.127): the Galaxy Map's opening
# view, the lists, the sector, system and phenomenon pages and the search, so
# the first visit after the update is fast. Each part is its own command, run
# in order by the wrapper, and the log has a line for each; the update does not
# wait for it.
#
# Steps 2, 3 and 6 share their code with install.sh (scripts/), so the
# two can't disagree about what a working server needs. `sudo
# ./install.sh` is still there for a full reinstall.
#
# Before it changes anything (OPS.41, OPS.53) it takes the update guard
# (`planetgen.cli.update_guard`): it lists every running or queued piece of
# work (Generate jobs, work queue runs, queued jobs, the orbital update, the
# Galaxy Map warm-up) and stops, changing nothing, unless you agree to cancel
# them (a question on the terminal; without one, e.g. the maintenance timer,
# add --cancel-running). Then it closes the site (maintenance mode, cause
# "update") and nothing else can start until the update has finished. A good
# update opens the site again at the end; a failed one leaves it closed, says
# which step failed and how to open it (`python3 -m planetgen.cli.maintenance
# off`), because the database may be half migrated.
#
# Usage:
#   sudo ./update.sh [--cancel-running]
#   sudo PYTHON=/usr/bin/python3.12 ./update.sh   (make the venv from another Python)
#
# Forces the checkout to match origin's branch tip even if there are
# uncommitted local changes to tracked files (a `git reset --hard` after
# fetching, rather than a `git pull`, which would otherwise refuse or
# conflict) -- a deployment server is meant to always run exactly what's
# on that branch, not whatever got hand-edited on it since. This never
# runs `git clean`, so untracked files -- most importantly `config.json`,
# which is gitignored precisely so it survives this -- are left alone.
#
# Runs on Linux (Apache with mod_wsgi) and macOS (gunicorn under launchd;
# docs/deployment/macos.md), both from the venv install-python-deps.sh
# makes, under macOS's bash 3.2 too.

set -euo pipefail

CANCEL_RUNNING=()
for arg in "$@"; do
    case "$arg" in
        --cancel-running) CANCEL_RUNNING=(--cancel-running) ;;
        *) echo "error: unknown argument '$arg' (the only option is --cancel-running)." >&2; exit 2 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HTML_DIR="$SCRIPT_DIR/src/html"
NLTK_DATA_DIR="${PLANETGEN_NLTK_DATA_DIR:-/usr/local/share/nltk_data}"

if [[ $EUID -ne 0 ]]; then
    echo "error: must be run as root, e.g.:" >&2
    echo "  sudo $0" >&2
    exit 1
fi

cd "$SCRIPT_DIR"

# The checkout belongs to Apache's user (examples/apache/set-permissions.sh),
# not to root, and git refuses a repository owned by someone else
# ("detected dubious ownership"). Trust this one checkout for this run only,
# through the environment rather than any git config file.
export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=safe.directory GIT_CONFIG_VALUE_0="$SCRIPT_DIR"

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

# OPS.41, OPS.53: the guard comes before the first change. The checkout may
# not have it yet (the first update that brings it), then it is taken right
# after the pull; without a virtual environment yet, after step 2.
UPDATE_STEP="starting"
export PLANETGEN_UPDATE_PID=$$
trap finish_update_guard EXIT
ensure_update_guard

UPDATE_STEP="step 1/9, pulling the latest changes"
echo "== 1/9: Pulling the latest changes =="

# Mode changes are not local edits: set-permissions.sh makes scripts and
# *.py files executable, which git (with core.fileMode on) would list.
dirty="$(git -c core.fileMode=false status --porcelain)"
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
ensure_update_guard
# The scripts below are run directly, so this comes before any of them.
find "$SCRIPT_DIR" -name '*.sh' -exec chmod +x {} +
find "$HTML_DIR" -name '*.py' -exec chmod +x {} +

# Again, now that the pull may have changed it.
# shellcheck source=scripts/deploy-common.sh
source "$SCRIPT_DIR/scripts/deploy-common.sh"

echo
UPDATE_STEP="step 2/9, checking the Python libraries"
echo "== 2/9: Checking the Python libraries =="
PYTHON="$PYTHON" bash "$SCRIPT_DIR/scripts/install-python-deps.sh" --check
use_site_python
ensure_update_guard required

echo
UPDATE_STEP="step 3/9, checking the NLTK words corpus"
echo "== 3/9: Checking the NLTK 'words' corpus =="
ensure_nltk_words "$NLTK_DATA_DIR"

echo
UPDATE_STEP="step 4/9, checking the generator's math"
echo "== 4/9: Checking the generator's math =="
check_math

echo
UPDATE_STEP="step 5/9, migrating the database"
echo "== 5/9: Migrating the configured MySQL database to the current schema =="
migrate_or_reset_db
record_version_key

echo
UPDATE_STEP="step 6/9, checking Apache's modules"
echo "== 6/9: Checking Apache's modules (macOS: the gunicorn daemon) =="
APACHE_NEEDS_RESTART=0
if is_macos; then
    ensure_gunicorn_daemon
else
    ensure_apache_modules
fi

echo
UPDATE_STEP="step 7/9, setting ownership and permissions"
echo "== 7/9: Setting ownership/permissions: the checkout belongs to Apache's user =="
PYTHON="$PYTHON" "$SCRIPT_DIR/examples/apache/set-permissions.sh" "$SCRIPT_DIR"

echo
UPDATE_STEP="step 8/9, checking the cache, jobs and log locations"
echo "== 8/9: Checking the cache, jobs and debug log locations =="
PYTHON="$PYTHON" "$SCRIPT_DIR/examples/apache/create-cache-dir.sh"
# The logs never stop an update (OPS.5): setup-debug-log.sh warns, with
# the commands that fix it, for anything it can't set up itself.
PYTHON="$PYTHON" "$SCRIPT_DIR/examples/apache/setup-debug-log.sh" \
    || echo "warning: the log setup didn't finish (see above); the update carries on." >&2

echo
UPDATE_STEP="step 9/9, checking that the web app imports"
echo "== 9/9: Checking that the web app imports =="
check_app_imports

echo
echo "Checking the Redis server (config.json's redis.url):"
ensure_redis
check_redis_settings

echo
warm_galaxy_map

echo
note_changed_example_configs "$before" "$after"
if [[ "$MATH_CHECK_FAILED" == 1 ]]; then
    echo "warning: the math check failed (step 4): bulk generation refuses to start until it passes." >&2
fi
echo "The update doesn't run the population pass (species, civilizations, territories); when wanted:"
echo "  planetgen population"
if is_macos; then
    if [[ "$before" != "$after" ]]; then
        echo "Done. Reload gunicorn so the site runs the new code (a running Generate job is left alone):"
        echo "  sudo launchctl kill SIGHUP system/$GUNICORN_LABEL"
    else
        echo "Done. Nothing new was pulled."
    fi
elif (( APACHE_NEEDS_RESTART )); then
    if reload_apache_if_running; then
        echo "Done. An Apache module was just enabled, so Apache was restarted to load it."
    else
        echo "Done. An Apache module was just enabled: restart Apache to load it:"
        echo "  sudo systemctl restart apache2"
    fi
elif [[ "$before" != "$after" ]]; then
    if reload_apache_if_running; then
        echo "Done. The site is running the new code."
    else
        echo "Done. Reload Apache so the site runs the new code:"
        echo "  sudo systemctl reload apache2"
    fi
else
    echo "Done. Nothing new was pulled."
fi
