#!/usr/bin/env bash
#
# update.sh
#
# Pulls the latest planetGen changes from git, then:
#   - If the pull actually brought new commits, re-runs `install.sh` so
#     everything it covers (the Python package, the NLTK corpus, Apache's
#     CGI module, permissions) stays correct afterward -- a `git pull` on
#     its own isn't enough: pulling a changed file rewrites it with
#     whatever mode is tracked in the repo (non-executable, historically
#     -- see `docs/TODO.md`'s "Deployment bugs found in production"
#     section), silently undoing any executable bit a previous
#     `install.sh`/`set-permissions.sh` run had fixed.
#   - Otherwise (already up to date), skips that -- there is nothing new
#     to reinstall, re-fetch, or re-`chmod`, so redoing all of it anyway
#     (a `pip install --force-reinstall`, Apache module/permission
#     churn, ...) on every single run would just be wasted work, which
#     matters for a scheduled/unattended caller (see
#     `examples/maintenance/`) that may run this monthly for years
#     without ever actually finding a new commit. Either way, `src/migrateDb.py`
#     still runs directly (see below) -- database migrations are cheap
#     and idempotent (a no-op once the schema is already current, per its
#     own docstring), and a database can need migrating even when this
#     checkout's code didn't just change (e.g. this is the first time
#     update.sh has run against it since install.sh's own last-run
#     migration).
#
# Usage:
#   sudo ./update.sh
#
# Refuses to run over uncommitted local changes (checks `git status`
# first) rather than stashing or discarding them for you -- if you've
# hand-edited something on this deployment, resolve that yourself first
# (commit, stash, or discard it deliberately) and re-run.
#
# Linux only -- same scope as install.sh/examples/apache/set-permissions.sh.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

echo "== 1/2: Pulling the latest changes =="

dirty="$(git status --porcelain)"
if [[ -n "$dirty" ]]; then
    echo "error: uncommitted local changes in $SCRIPT_DIR -- refusing to pull over them:" >&2
    echo "$dirty" >&2
    echo "Commit, stash, or discard these yourself, then re-run $0." >&2
    exit 1
fi

branch="$(git rev-parse --abbrev-ref HEAD)"
if [[ "$branch" == "HEAD" ]]; then
    echo "error: repository is in a detached HEAD state -- check out a branch first." >&2
    exit 1
fi

before="$(git rev-parse HEAD)"
git fetch origin "$branch"
# --ff-only rather than a plain `pull`: a deployment server should never
# end up with a surprise merge commit. If history has diverged, fail
# loudly here rather than silently merging or, worse, needing a manual
# conflict resolution on a live server.
git pull --ff-only origin "$branch"
after="$(git rev-parse HEAD)"

if [[ "$before" == "$after" ]]; then
    echo "Already up to date ($before)."
else
    echo "Updated $before..$after:"
    git log --oneline "$before..$after"
fi

echo
if [[ "$before" == "$after" ]]; then
    echo "== 2/2: No new commits -- checking the database schema only =="
    # Nothing changed on disk, so nothing needs reinstalling/re-chmod-ing:
    # skip straight to the one step that isn't conditional on the code
    # having changed. Uses the package/interpreter already installed from
    # a previous install.sh run -- there's nothing to reinstall it from
    # here, since this branch is exactly "no new commit landed".
    PYTHON="$(command -v python3 || command -v python || true)"
    if [[ -z "$PYTHON" ]]; then
        echo "error: no python3/python found on PATH." >&2
        exit 1
    fi
    "$PYTHON" "$SCRIPT_DIR/src/migrateDb.py"
else
    echo "== 2/2: Re-running install.sh to keep permissions (and everything else it covers) correct =="
    # A pull rewrites any changed file with whatever mode is tracked in the
    # repo -- including install.sh (and this script) itself -- so a prior
    # run's executable-bit fix doesn't survive a pull that touched them.
    # Fixing that here, before invoking install.sh, matters because
    # install.sh is run directly below ("$SCRIPT_DIR/install.sh", not
    # `bash install.sh`): if the pull just dropped its executable bit,
    # install.sh's own step 4 (which re-chmods every *.sh in the repo) never
    # gets a chance to run at all -- the shell refuses to exec it first with
    # "Permission denied", exactly as install.sh's own step 4 fix already
    # had to for src/html/*.py. install.sh's own step 2 covers the database
    # migration in this branch, so it isn't run a second time here.
    find "$SCRIPT_DIR" -name '*.sh' -exec chmod +x {} +
    "$SCRIPT_DIR/install.sh"
fi
