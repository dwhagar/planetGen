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
#   3. The NLTK 'words' corpus: fetched only if it's missing. Then the static
#      data sets (GEN.548): `planetgen.cli.datasets fetch` downloads each file
#      once, checks its sha256 and keeps it in /var/cache/planetgen/datasets
#      (PLANETGEN_DATA_DIR), outside the checkout. A host that is down only
#      warns; bulk generation refuses to start until the next update has
#      fetched it (GEN.590). Every file's https URL is printed with what
#      happened to it (fetched, kept, failed and why), so the output says
#      where each piece comes from (GEN.600). Then the data of the References
#      pages (`scripts/build_reference_pages.py`, DOC.88), built from the
#      reference list and the stand-in registry; a failure only warns.
#   4. `planetgen check-math`: the math check (TEST.68). A failure only
#      warns (repeated at the end). Generation does not run it (PERF.170); the
#      admin page can run it any time.
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
# At the end it reloads what serves the site when that is running (OPS.8, OPS.29): Apache
# (after its own configtest; restarted only when step 6 just installed mod_wsgi), the
# planetgen-gunicorn unit on Linux, or the launchd gunicorn daemon on macOS, then checks
# that /api/health answers, at config.json's api_base_url and nowhere else (OPS.86: no address is guessed; with none configured the
# update stops before it changes the database). A failed step ends the update non-zero and prints the command
# to run by hand; a reload that worked but whose /api/health probe never got a good answer only warns and the site opens; nothing running prints the command and ends zero. A reload aborts a
# mod_wsgi request still running after about 4 seconds; touching src/html/wsgi.py does not.
#
# Once the code, the database, the permissions and Redis are in place it starts building
# the site's caches (`planetgen warm`, MAP.134, PERF.104 and PERF.127) in two stages, each in the
# foreground with a progress bar: stage 1 builds only the opening phase (the cluster field, the object index, the
# opening views, the first page of each list), what the reload and the health check need, then it waits 5 seconds
# (UPDATE_SETTLE_SECONDS), reloads, waits 10 seconds (HEALTH_SETTLE_SECONDS) and checks /api/health; stage 2 is the
# full warm-up (the sector, system and phenomenon pages, the search and the rest), which the update waits for, showing
# its bar, so the reload cannot interrupt it.
#
# Steps 2, 3 and 6 share their code with install.sh (scripts/), so the
# two can't disagree about what a working server needs. `sudo
# ./install.sh` is still there for a full reinstall.
#
# Before it changes anything (OPS.41, OPS.53) it takes the update guard
# (`planetgen.cli.update_guard`, DB.29): it lists every running or queued piece
# of work (Generate jobs, runs in the run record, queued jobs, the orbital
# update, the Galaxy Map warm-up) and stops, changing nothing, unless you
# agree to cancel them: --yes, or a y typed at the question on the terminal.
# Nobody to ask is a no: a run started by systemd or cron, one with no
# terminal, and 30 seconds without an answer all stop it, with exit status 1,
# unless --yes was given (the maintenance timer adds it when it should cancel
# the work). Then it closes the site (maintenance mode, cause
# "update") and nothing else can start until the update has finished. A good
# update opens the site again at the end; a failed one leaves it closed, says
# which step failed and how to open it (`python3 -m planetgen.cli.maintenance
# off`), because the database may be half migrated.
#
# A job that is stuck and will not stop (the guard waits for it and refuses)
# is cleared with --reset-jobs: before the guard it stops everything it lists,
# kills the queue workers and empties the jobs table, the work lease and the
# Redis queues (`planetgen.cli.reset_jobs`, which can also be run alone). It
# asks nothing and cancels what is running, so it is the word to use when you
# know nothing needs to finish.
#
# --kill-all (PERF.186) is the full stop: before the guard it runs `planetgen queue stop-all --no-hold`, the same
# routine as the Stop everything button on the Queue page. It cancels every run and the Generate page's job, empties
# every queue, and ends every planetgen worker and command process on the host with the processes they started,
# including the ones nothing lists any more, and implies --yes. It is only run when you give it. It does
# not leave the stop-all hold on (the update's own lock keeps new work out, and the site works after it).
#
# Usage:
#   sudo ./update.sh [--yes] [--reset-jobs] [--kill-all]
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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/logrecord.sh
source "$SCRIPT_DIR/scripts/logrecord.sh"

ANSWER_YES=()
RESET_JOBS=0
KILL_ALL=0
for arg in "$@"; do
    case "$arg" in
        --yes) ANSWER_YES=(--yes) ;;
        --reset-jobs) RESET_JOBS=1 ;;
        --kill-all) KILL_ALL=1; ANSWER_YES=(--yes) ;;
        *) log_error "unknown argument '$arg' (the options are --yes, --reset-jobs and --kill-all)."; exit 2 ;;
    esac
done

HTML_DIR="$SCRIPT_DIR/src/html"
NLTK_DATA_DIR="${PLANETGEN_NLTK_DATA_DIR:-/usr/local/share/nltk_data}"

if [[ $EUID -ne 0 ]]; then
    log_error "must be run as root, e.g.:"
    log_error "  sudo $0"
    exit 1
fi

cd "$SCRIPT_DIR"

# The checkout belongs to Apache's user (examples/apache/set-permissions.sh),
# not to root, and git refuses a repository owned by someone else
# ("detected dubious ownership"). Trust this one checkout for this run only,
# through the environment rather than any git config file.
export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=safe.directory GIT_CONFIG_VALUE_0="$SCRIPT_DIR"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    log_error "$SCRIPT_DIR is not a git checkout -- can't pull an update here."
    exit 1
fi

# shellcheck source=scripts/deploy-common.sh
source "$SCRIPT_DIR/scripts/deploy-common.sh"

PYTHON="${PYTHON:-$(default_python)}"
if [[ -z "$PYTHON" ]]; then
    log_error "no python3/python found on PATH."
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
log_stage 1 9 "Pulling the latest changes"

# Mode changes are not local edits: set-permissions.sh makes scripts and
# *.py files executable, which git (with core.fileMode on) would list.
# The state/ folder (the maintenance switch, the update lock) is the update's own, never a local edit (OPS.76).
dirty="$(git -c core.fileMode=false status --porcelain -- . ':(exclude)state')"
if [[ -n "$dirty" ]]; then
    log_warning "uncommitted local changes in $SCRIPT_DIR -- these will be overwritten:"
    log_warning "$dirty"
fi

branch="$(git rev-parse --abbrev-ref HEAD)"
if [[ "$branch" == "HEAD" ]]; then
    log_error "repository is in a detached HEAD state -- check out a branch first."
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
    log_info "Already up to date ($before)."
else
    log_info "Updated $before..$after:"
    git log --oneline "$before..$after"
fi
ensure_update_guard
# The scripts below are run directly, so this comes before any of them.
find "$SCRIPT_DIR" -name '*.sh' -exec chmod +x {} +
find "$HTML_DIR" -name '*.py' -exec chmod +x {} +

# Again, now that the pull may have changed it.
# shellcheck source=scripts/deploy-common.sh
source "$SCRIPT_DIR/scripts/deploy-common.sh"

warn_if_no_swap

echo
UPDATE_STEP="step 2/9, checking the Python libraries"
log_stage 2 9 "Checking the Python libraries"
PYTHON="$PYTHON" bash "$SCRIPT_DIR/scripts/install-python-deps.sh" --check
use_site_python
ensure_update_guard required
# OPS.86: the address the update asks whether the site answers is config.json's api_base_url; missing, it stops here
# (before the database is touched) and says so, instead of falling back to a loopback address.
resolve_api_health || { log_error "set api_base_url in config.json, then run the update again."; exit 1; }
log_info "The API address (api_base_url in config.json) is ${API_HEALTH_URL%/health}; the update asks only that address."

echo
UPDATE_STEP="step 3/9, checking the NLTK words corpus and the downloaded data"
log_stage 3 9 "Checking the NLTK 'words' corpus and the downloaded data sets"
ensure_nltk_words "$NLTK_DATA_DIR"
ensure_datasets
build_reference_pages

echo
UPDATE_STEP="step 4/9, checking the generator's math"
log_stage 4 9 "Checking the generator's math"
check_math

echo
UPDATE_STEP="step 5/9, migrating the database"
log_stage 5 9 "Migrating the configured MySQL database to the current schema"
migrate_or_reset_db
record_version_key
warm_shape_means

echo
UPDATE_STEP="step 6/9, checking Apache's modules"
log_stage 6 9 "Checking Apache's modules (macOS: the gunicorn daemon)"
APACHE_NEEDS_RESTART=0
if is_macos; then
    ensure_gunicorn_daemon
else
    ensure_apache_modules
fi

echo
UPDATE_STEP="step 7/9, setting ownership and permissions"
log_stage 7 9 "Setting ownership/permissions: the checkout belongs to Apache's user"
PYTHON="$PYTHON" "$SCRIPT_DIR/examples/apache/set-permissions.sh" "$SCRIPT_DIR"

echo
UPDATE_STEP="step 8/9, checking the cache, jobs and log locations"
log_stage 8 9 "Checking the cache, jobs and debug log locations"
PYTHON="$PYTHON" "$SCRIPT_DIR/examples/apache/create-cache-dir.sh"
# The logs never stop an update (OPS.5): setup-debug-log.sh warns, with
# the commands that fix it, for anything it can't set up itself.
PYTHON="$PYTHON" "$SCRIPT_DIR/examples/apache/setup-debug-log.sh" \
    || log_warning "the log setup didn't finish (see above); the update carries on."

# Set the secret key if none is (never printed); as Apache's user so the web app can read it.
set_secret_key

echo
UPDATE_STEP="step 9/9, checking that the web app imports"
log_stage 9 9 "Checking that the web app imports"
check_app_imports

echo
log_info "Checking the Redis server (config.json's redis.url):"
ensure_redis
check_redis_settings

echo
note_changed_example_configs "$before" "$after"
if [[ "$MATH_CHECK_FAILED" == 1 ]]; then
    log_warning "the math check failed (step 4): charted data may be wrong until it is fixed."
fi
if [[ "$DATASETS_MISSING" == 1 ]]; then
    log_warning "some data sets could not be downloaded (step 3): bulk generation refuses to start until they are."
fi
if [[ "${REFERENCE_PAGES_FAILED:-0}" == 1 ]]; then
    log_warning "the References pages could not be built (step 3): /references may show old numbers."
fi
log_info "The update doesn't run the population pass (species, civilizations, territories); when wanted:"
log_info "  planetgen population && planetgen territories"
if [[ "$before" == "$after" ]] && (( ! APACHE_NEEDS_RESTART )); then
    log_notice "Done. Nothing new was pulled."
else
    UPDATE_STEP="building what the reload and the health check need (cache warm-up, stage 1)"
    run_warm_stage "stage 1 of 2, what the reload and the health check need" opening
    echo "The update steps are finished: waiting ${UPDATE_SETTLE_SECONDS:-5} seconds before reloading the site."
    sleep "${UPDATE_SETTLE_SECONDS:-5}"
    UPDATE_STEP="reloading what serves the site and checking that it answers /api/health"
    reload_status=0
    reload_web_if_running || reload_status=$?
    if (( reload_status == 0 )); then
        if (( APACHE_NEEDS_RESTART )); then
            log_notice "Done. mod_wsgi was just installed, so Apache was restarted to load it."
        else
            log_notice "Done. The site is running the new code."
        fi
    elif (( reload_status == 1 )); then
        log_error "Done, but the site was not reloaded cleanly (see the warning above); the old code may still be serving."
        print_reload_commands >&2
        exit 1
    elif (( reload_status == 3 )); then
        # The migration finished (step 5) and the reload itself succeeded: only the probe got no good answer, which
        # can be the address's fault (another virtual host, no https listener) as much as the app's. The site opens.
        log_warning "the site was reloaded but /api/health did not answer well (see the warning above); the database is migrated and current, so the site is opened. Check it by hand."
    else
        log_notice "Done. Nothing that serves the site is running, so nothing was reloaded. When it is, run:"
        print_reload_commands
    fi
fi

# Stage 2 (Boss, 11-01 EDT): the full warm-up, after the reload and the health check, so nothing can interrupt it; the
# update stays up, showing its bar, until it ends. The update is over as far as the site goes (the guard is released and
# the site opens) before it starts, so the site is not closed for the length of the warm-up.
if (( UPDATE_GUARD_HELD )); then
    update_guard end --pid "$$" && UPDATE_GUARD_HELD=0
fi
echo
UPDATE_STEP="building the site's caches (cache warm-up, stage 2)"
run_warm_stage "stage 2 of 2, the full warm-up"
