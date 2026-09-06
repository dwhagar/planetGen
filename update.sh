#!/usr/bin/env bash
#
# update.sh
#
# Pulls the latest planetGen changes from git and re-runs `install.sh` so
# permissions (and everything else `install.sh` covers -- the Python
# package, the NLTK corpus, Apache's CGI module) stay correct afterward.
# This exists because a `git pull` on its own isn't enough: pulling a
# changed file rewrites it with whatever mode is tracked in the repo
# (non-executable, historically -- see `TODO.md`'s "Deployment bugs
# found in production" section), silently undoing any executable bit a
# previous `install.sh`/`set-permissions.sh` run had fixed.
#
# Usage:
#   sudo ./update.sh
#
# Refuses to run over uncommitted local changes (checks `git status`
# first) rather than stashing or discarding them for you -- if you've
# hand-edited something on this deployment, resolve that yourself first
# (commit, stash, or discard it deliberately) and re-run.
#
# Linux only -- same scope as install.sh/apache/set-permissions.sh.

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
echo "== 2/2: Re-running install.sh to keep permissions (and everything else it covers) correct =="
"$SCRIPT_DIR/install.sh"
