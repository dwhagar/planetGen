#!/usr/bin/env bash
#
# scripts/lock-requirements.sh
#
# Rewrites requirements.lock: every library planetGen needs at run time
# (setup.py's install_requires plus its 'api' extra) and everything those
# need in turn, each pinned to one version and to the sha256 hashes of
# the files PyPI publishes for it. scripts/install-python-deps.sh installs
# only those exact files whenever it uses pip (`--require-hashes`), so a
# release that appears on PyPI later, or a file swapped under the same
# version, never reaches a server by surprise. apt-installed libraries
# aren't covered: apt checks its own signatures.
#
# The lock is universal: one file for Linux and macOS and every
# Python from setup.py's python_requires (3.9) up. Where the newest
# release of a library needs a newer Python, the lock holds one pin per
# Python range, each with an environment marker.
#
# Run it after changing a requirement in setup.py (a test fails until the
# lock satisfies setup.py's floors), or with --upgrade to move every pin
# to the newest release. Without --upgrade, pins already in the lock are
# kept wherever they still satisfy setup.py.
#
# Needs uv (https://docs.astral.sh/uv/): `pipx install uv` or
# `python3 -m pip install --user uv`. Not needed on a server.
#
# Usage:
#   scripts/lock-requirements.sh [--upgrade]

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

UPGRADE=()
case "${1:-}" in
    "") ;;
    --upgrade) UPGRADE=(--upgrade) ;;
    *)
        echo "error: unknown argument '$1' (expected nothing or --upgrade)." >&2
        exit 1 ;;
esac

if ! command -v uv >/dev/null 2>&1; then
    echo "error: uv not found; install it with 'pipx install uv' or 'python3 -m pip install --user uv'." >&2
    exit 1
fi

cd "$ROOT"
uv pip compile setup.py --extra api \
    --universal --python-version 3.9 --generate-hashes \
    --custom-compile-command "scripts/lock-requirements.sh" \
    ${UPGRADE[@]+"${UPGRADE[@]}"} \
    --output-file requirements.lock
uv pip compile setup.py --extra api --extra server \
    --universal --python-version 3.9 --generate-hashes \
    --custom-compile-command "scripts/lock-requirements.sh" \
    ${UPGRADE[@]+"${UPGRADE[@]}"} \
    --output-file requirements-server.lock
