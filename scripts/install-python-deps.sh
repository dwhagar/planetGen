#!/usr/bin/env bash
#
# scripts/install-python-deps.sh
#
# Gives planetGen its own Python: a virtual environment (PLANETGEN_VENV,
# default /usr/local/planetgen/venv) made from the system's python3 (on
# macOS, Homebrew's), with every library planetGen needs, and the app
# server (gunicorn), installed into it from requirements-server.lock, the
# exact versions scripts/lock-requirements.sh pinned, checked against the
# lock's sha256 hashes (--require-hashes), so nothing newer or altered on
# PyPI gets in unnoticed. Nothing is installed into the system Python, and
# the system Python's own packages (apt's python3-* libraries) are never
# used or touched: the venv does not see them (OPS.45).
#
# Everything runs from the venv's python:
#   - the web app: examples/apache/planetgen.conf.example gives
#     WSGIDaemonProcess `python-home=/usr/local/planetgen/venv`, and the
#     gunicorn units and plists name the venv's gunicorn;
#   - the command-line tools and background jobs: the `planetgen` command
#     (/usr/local/bin/planetgen links to the venv's) and
#     `/usr/local/planetgen/venv/bin/python -m planetgen.cli.NAME`;
#   - install.sh and update.sh themselves, once this script has run.
#
# planetGen itself is installed as an editable package of this checkout
# (`pip install -e`, install_checkout below), so the web app
# (src/html/wsgi.py) and every command-line tool import the checkout's
# code from anywhere, and a `git pull` takes effect with no reinstall.
#
# mod_wsgi embeds the Python it was built against: the venv must be made
# from that same version (deploy-common.sh's check_mod_wsgi_python warns
# when it is not). Debian and Ubuntu need the python3-venv package to
# make one; it is installed here when `python3 -m venv` fails without it.
#
# Usage (as root):
#   scripts/install-python-deps.sh           install.sh: full install
#   scripts/install-python-deps.sh --check   update.sh: install nothing that
#                                            is already there
#
# --check imports every requirement with the venv's python and compares
# its version with the floor, and only when one is missing, below its
# floor or fails to import does it install from the lock again. Both
# actions print one line per requirement and exit non-zero if anything is
# still unusable. Both also make the editable install of this checkout and
# the `planetgen` command if either is missing or points elsewhere.
#
# Runs under macOS's bash 3.2 as well as Linux's bash: no associative
# arrays or mapfile, and every array that can be empty is expanded as
# ${a[@]+"${a[@]}"} (bash before 4.4 calls an empty "${a[@]}" unbound
# under set -u).
#
# Environment:
#   PYTHON          the interpreter the venv is made from (default: python3 on
#                   PATH; on macOS, Homebrew's python3). Python 3.12 or newer.
#   PLANETGEN_VENV  the venv's path (default /usr/local/planetgen/venv)

set -euo pipefail

ACTION=install
case "${1:-}" in
    "") ;;
    --check) ACTION=check ;;
    *)
        echo "error: unknown argument '$1' (expected nothing or --check)." >&2
        exit 1 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -z "${PYTHON:-}" && "$(uname -s)" == Darwin ]]; then
    # Homebrew's python3 (Apple Silicon, then Intel), not the Command Line
    # Tools' /usr/bin/python3. Homebrew itself won't run as root, so it
    # is never called here.
    for candidate in /opt/homebrew/bin/python3 /usr/local/bin/python3; do
        if [[ -x "$candidate" ]]; then PYTHON="$candidate"; break; fi
    done
fi
PYTHON="${PYTHON:-$(command -v python3 || command -v python || true)}"
VENV="${PLANETGEN_VENV:-/usr/local/planetgen/venv}"
PLANETGEN_COMMAND=/usr/local/bin/planetgen
# Exact versions and file hashes for everything pip installs
# (scripts/lock-requirements.sh writes both): the libraries plus the app
# server (gunicorn).
SERVER_LOCK="$SCRIPT_DIR/requirements-server.lock"

# Every runtime requirement install.sh needs: setup.py's install_requires
# plus its 'api' extra. Keep in step with setup.py
# (src/tests/test_install_python_deps.py checks this).
# Each is imported by its pip name with "-" turned into "_"
# (flask-limiter -> flask_limiter), except where
# scripts/probe_requirements.py's IMPORT_NAMES says otherwise
# (scikit-image -> skimage).
REQUIREMENTS=(
    "nltk>=3.9.1"
    "pymysql>=1.1.1"
    "cryptography>=42.0.0"
    "werkzeug>=3.0.0"
    "rich>=13.7.0"
    "flask>=3.0.3"
    "flask-limiter>=3.7.0"
    "redis>=5.0.0"
    "rq>=1.16.0"
    "pyotp>=2.9.0"
    "segno>=1.6.0"
    "markdown>=3.6"
    "cachetools>=5.3.0"
    "sqlalchemy>=2.0.30"
    "alembic>=1.13.0"
    "pydantic>=2.7.0"
    "numpy>=1.26.0"
    "scipy>=1.13.0"
    "scikit-image>=0.22.0"
    "numba>=0.68.0"
)

if [[ -z "$PYTHON" ]]; then
    echo "error: no python3/python found on PATH." >&2
    exit 1
fi

SPECS=()
for line in "${REQUIREMENTS[@]}"; do
    SPECS+=("$line")
done

# Prints "<state> <spec> <version or error> <directory>" for each
# requirement (from the arguments) as the venv's python sees it: state is
# ok, missing (not installed), old (below its floor) or broken (installed
# but the import raised). Imports each one for real, since an installed
# distribution whose own dependencies are missing is not usable either.
probe() {
    "$VENV/bin/python" "$SCRIPT_DIR/scripts/probe_requirements.py" "$@"
}

# Prints one line per requirement from the probe lines on stdin. Returns
# 1 if any requirement is unusable.
report() {
    local state spec detail where failed=0
    while read -r state spec detail where; do
        if [[ "$state" != ok ]]; then
            printf '  %-10s %s (%s: %s)\n' failed "$spec" "$state" "$detail"
            failed=1
        else
            printf '  %-10s %-24s %-8s (pip, in %s)\n' ok "$spec" "$detail" "$VENV"
        fi
    done
    return "$failed"
}

# Probes every requirement (and gunicorn) and reports; fails if any is
# unusable or the probe itself didn't run.
final_report() {
    local after
    after="$(probe "${SPECS[@]}" "gunicorn>=0")" || true
    if (( $(grep -c . <<< "$after") != ${#SPECS[@]} + 1 )); then
        echo "error: could not check the Python libraries with $VENV/bin/python." >&2
        return 1
    fi
    if ! report <<< "$after"; then
        echo "error: some Python libraries are still unusable (see above)." >&2
        return 1
    fi
}

# Installs this checkout into the venv as an editable package (`pip
# install -e`, no dependencies: those come from the lock), so
# `planetgen` imports from anywhere, the web app and every command-line
# tool run the checkout's code, and a `git pull` takes effect without a
# reinstall. pip writes the `planetgen` console script into the venv's
# bin; /usr/local/bin/planetgen links to it so the command is on PATH.
# Skipped when the package already imports from this checkout and the
# command is in place, so update.sh reinstalls nothing.
install_checkout() {
    local scripts
    scripts="$("$VENV/bin/python" -c 'import sysconfig; print(sysconfig.get_path("scripts"))')"
    if ! (cd / && "$VENV/bin/python" - "$SCRIPT_DIR/src" <<'EOF2'
import os, sys
import planetgen
sys.exit(0 if os.path.dirname(os.path.dirname(os.path.realpath(planetgen.__file__))) == os.path.realpath(sys.argv[1]) else 1)
EOF2
    ) 2>/dev/null || [[ ! -x "$scripts/planetgen" ]]; then
        echo "Installing planetGen from $SCRIPT_DIR as an editable package."
        "$VENV/bin/python" -m pip install --quiet --no-deps -e "$SCRIPT_DIR"
    fi
    if [[ "$(readlink "$PLANETGEN_COMMAND" 2>/dev/null || true)" != "$scripts/planetgen" ]]; then
        mkdir -p "$(dirname "$PLANETGEN_COMMAND")"
        ln -sf "$scripts/planetgen" "$PLANETGEN_COMMAND"
        echo "Linked $PLANETGEN_COMMAND to $scripts/planetgen."
    fi
}

# Makes the venv from $PYTHON, installing Debian's python3-venv when
# `-m venv` can't (no ensurepip) and apt is there to do it.
make_venv() {
    mkdir -p "$(dirname "$VENV")"
    if "$PYTHON" -m venv "$VENV" 2>/dev/null; then return 0; fi
    rm -rf "$VENV"
    if command -v apt-get >/dev/null 2>&1; then
        echo "python3 -m venv failed: installing python3-venv."
        apt-get update -qq || true
        DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3-venv || true
    fi
    "$PYTHON" -m venv "$VENV"
}

# Creates the venv if it is not there and installs the locked libraries
# and gunicorn into it. pip leaves alone whatever is already installed at
# the locked version.
install_venv() {
    # setup.py's python_requires.
    "$PYTHON" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 12) else 1)' || {
        echo "error: $PYTHON is older than Python 3.12. Install a newer python3 (macOS: brew install python@3.12)," >&2
        echo "  and give it to this script as PYTHON=..." >&2
        return 1
    }
    if [[ ! -x "$VENV/bin/python" ]]; then
        echo "Creating the virtual environment $VENV with $PYTHON."
        make_venv
        "$VENV/bin/python" -m pip install --quiet --upgrade pip
    fi
    echo "Installing the locked libraries and gunicorn into $VENV."
    "$VENV/bin/python" -m pip install --quiet --require-hashes -r "$SERVER_LOCK"
}

venv_ready() {
    [[ -x "$VENV/bin/python" ]] || return 1
    local status
    status="$(probe "${SPECS[@]}" "gunicorn>=0")" || return 1
    ! grep -qv '^ok ' <<< "$status"
}

if [[ "$ACTION" == check ]]; then
    echo "Checking the Python libraries in $VENV."
    if ! venv_ready; then
        echo "Missing, too old or not importable in $VENV: installing from requirements-server.lock."
        install_venv
    fi
    install_checkout
    if ! final_report; then
        echo "  Running sudo ./install.sh does a full reinstall." >&2
        exit 1
    fi
else
    install_venv
    install_checkout
    echo "Python install path: venv ($VENV)."
    final_report
fi
