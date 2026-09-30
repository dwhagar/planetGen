#!/usr/bin/env bash
#
# scripts/install-python-deps.sh
#
# Puts every Python library planetGen needs into one dedicated virtual
# environment, /opt/planetgen/venv (PLANETGEN_VENV_DIR), the usual way to
# deploy a Python service on a Linux host. The system Python is never
# changed: nothing is pip-installed into /usr/lib or /usr/local, and
# `--break-system-packages` is never used, so apt/dpkg stay in charge of
# everything they installed (PEP 668). apt still provides the OS side:
# the base interpreter, python3-venv, and mod_wsgi (see
# scripts/deploy-common.sh).
#
# The venv is built from the system python3 (PYTHON), because mod_wsgi
# (libapache2-mod-wsgi-py3) embeds that same Python and runs the venv via
# `python-home=/opt/planetgen/venv` on WSGIDaemonProcess (see
# examples/apache/planetgen.conf.example). It's isolated (no
# --system-site-packages): everything the site imports comes from here.
#
# planetGen itself isn't installed into the venv: every entry point adds
# the checkout's src/ to sys.path itself (generate.py,
# src/html/lib/apiclient.py, src/html/wsgi.py, and src/migrateDb.py
# through its own sys.path[0]), so the code that runs is always the code
# that was just pulled. /usr/local/bin/planetgen is a small wrapper that
# runs the checkout's generate.py with the venv's Python.
#
# Usage (as root):
#   scripts/install-python-deps.sh           install.sh: create the venv if
#                                            needed, install every library
#                                            at its latest version
#   scripts/install-python-deps.sh --check   update.sh: install nothing
#                                            that is already there
#
# --check imports every requirement with the venv's Python and compares
# its version with setup.py's floor. Only a requirement that is missing,
# below its floor or fails to import gets installed. The venv is rebuilt
# only when it's missing, its Python no longer runs, or its Python
# version no longer matches the system one (a distribution upgrade),
# since mod_wsgi would then load it with the wrong Python. Both actions
# print one line per requirement and exit non-zero if anything is still
# unusable.
#
# Environment:
#   PYTHON              base interpreter the venv is built from (default:
#                       python3 on PATH -- the one mod_wsgi is built for)
#   PLANETGEN_VENV_DIR  the venv (default /opt/planetgen/venv)

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
PYTHON="${PYTHON:-$(command -v python3 || command -v python || true)}"
VENV_DIR="${PLANETGEN_VENV_DIR:-/opt/planetgen/venv}"
VENV_PYTHON="$VENV_DIR/bin/python"
WRAPPER=/usr/local/bin/planetgen

# Every runtime requirement: setup.py's install_requires plus its 'api'
# extra. Keep in step with setup.py (src/tests/test_install_python_deps.py
# checks this). Each is imported by its pip name with "-" turned into "_"
# (flask-limiter -> flask_limiter).
REQUIREMENTS=(
    "nltk>=3.9.1"
    "pymysql>=1.1.1"
    "dbutils>=3.1.0"
    "werkzeug>=3.0.0"
    "rich>=13.7.0"
    "flask>=3.0.3"
    "flask-limiter>=3.7.0"
)

if [[ -z "$PYTHON" ]]; then
    echo "error: no python3/python found on PATH." >&2
    exit 1
fi

py_version() {
    "$1" -c 'import sys; print("%d.%d" % sys.version_info[:2])' 2>/dev/null
}

# Creates the venv when it's missing, its Python doesn't run, its Python
# version isn't the system one's any more, or it was made with
# --system-site-packages (the old fallback venv). Otherwise leaves it
# exactly as it is.
ensure_venv() {
    local want have
    want="$(py_version "$PYTHON")"
    have="$(py_version "$VENV_PYTHON" || true)"
    local shared=0
    grep -qi '^include-system-site-packages *= *true' "$VENV_DIR/pyvenv.cfg" 2>/dev/null && shared=1
    if [[ -n "$have" && "$have" == "$want" ]] && (( ! shared )); then
        echo "Using the venv at $VENV_DIR (Python $have)."
        return 0
    fi
    if (( shared )); then
        # PR #99's venv saw apt's packages and held only what apt lacked.
        echo "The venv at $VENV_DIR shares the system's packages: rebuilding it on its own."
    elif [[ -n "$have" ]]; then
        echo "The venv at $VENV_DIR is Python $have but $PYTHON is $want: rebuilding it."
    elif [[ -e "$VENV_DIR" ]]; then
        echo "The venv at $VENV_DIR doesn't run: rebuilding it."
    else
        echo "Creating the venv at $VENV_DIR (Python $want)."
    fi
    # python3-venv provides ensurepip, which `-m venv` needs for pip.
    if ! "$PYTHON" -c "import ensurepip" >/dev/null 2>&1 && command -v apt-get >/dev/null 2>&1; then
        DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3-venv
    fi
    mkdir -p "$(dirname "$VENV_DIR")"
    "$PYTHON" -m venv --clear "$VENV_DIR"
    "$VENV_PYTHON" -m pip install --quiet --upgrade pip
    VENV_REBUILT=1
}

# Prints "<state> <spec> <version or error>" for each requirement (from
# the arguments) as the venv's Python sees it: state is ok, missing (not
# installed), old (below its floor) or broken (installed but the import
# raised). Imports each one for real, since an installed distribution
# whose own dependencies are missing is not usable either.
probe() {
    "$VENV_PYTHON" - "$@" <<'EOF'
import importlib
import re
import sys
from importlib import metadata


def key(version):
    parts = []
    for piece in version.split("."):
        match = re.match(r"\d+", piece)
        if not match:
            break
        parts.append(int(match.group()))
        if match.group() != piece:
            break
    return tuple(parts)


for spec in sys.argv[1:]:
    name, floor = spec.split(">=")
    try:
        version = metadata.version(name)
    except metadata.PackageNotFoundError:
        print("missing", spec, "-")
        continue
    if key(version) < key(floor):
        print("old", spec, version)
        continue
    try:
        importlib.import_module(name.replace("-", "_"))
    except Exception as exc:  # any import failure makes it unusable
        print("broken", spec, f"{type(exc).__name__}: {exc}".replace("\n", " "))
        continue
    print("ok", spec, version)
EOF
}

# Prints one line per requirement from the probe lines in $2 (a file).
# Given the probe lines from before anything was installed in $1, labels
# each one present/installed/upgraded/repaired rather than ok. Returns 1
# if any requirement is still unusable.
report() {
    local before_file="$1" after_file="$2" state spec detail old label failed=0 i=0
    local -a before=()
    [[ -n "$before_file" ]] && mapfile -t before < "$before_file"
    while read -r state spec detail; do
        if [[ "$state" != ok ]]; then
            printf '  %-10s %s (%s: %s)\n' failed "$spec" "$state" "$detail"
            failed=1
        else
            label=ok
            if (( ${#before[@]} )); then
                read -r old _ <<< "${before[$i]}"
                case "$old" in
                    ok) label=present ;; missing) label=installed ;;
                    old) label=upgraded ;; *) label=repaired ;;
                esac
            fi
            printf '  %-10s %-24s %s\n' "$label" "$spec" "$detail"
        fi
        i=$((i + 1))
    done < "$after_file"
    return "$failed"
}

# Removes what earlier versions of this script put into the system
# Python: PR #99's planetgen-venv.pth, which put the venv on every system
# Python process's sys.path. The venv is reached through python-home and
# the wrapper now, so the system Python goes back to being the
# distribution's own.
remove_legacy_pth() {
    local pth
    for pth in /usr/local/lib/python3*/dist-packages/planetgen-venv.pth; do
        if [[ -f "$pth" ]]; then
            rm -f "$pth"
            echo "Removed $pth (the system Python no longer loads the venv)."
        fi
    done
}

# Stands in for the `planetgen` console script pip would make, but runs
# the checkout's generate.py with the venv's Python. Rewritten only when
# it differs.
write_wrapper() {
    local want
    want="$(cat <<EOF
#!/bin/sh
# Written by planetGen's scripts/install-python-deps.sh.
exec "$VENV_PYTHON" "$SCRIPT_DIR/generate.py" "\$@"
EOF
)"
    if [[ "$(cat "$WRAPPER" 2>/dev/null || true)" != "$want" ]]; then
        printf '%s\n' "$want" > "$WRAPPER"
        echo "Wrote $WRAPPER (runs $SCRIPT_DIR/generate.py with $VENV_PYTHON)."
    fi
    chmod 755 "$WRAPPER"
}

VENV_REBUILT=0
remove_legacy_pth
ensure_venv

before_file="$(mktemp)"
after_file="$(mktemp)"
trap 'rm -f "$before_file" "$after_file"' EXIT

if [[ "$ACTION" == install ]]; then
    "$VENV_PYTHON" -m pip install --upgrade "${REQUIREMENTS[@]}"
else
    echo "Checking the Python libraries in $VENV_DIR."
    probe "${REQUIREMENTS[@]}" > "$before_file"
    mapfile -t need < <(awk '$1 != "ok" {print $2}' "$before_file")
    if (( ${#need[@]} )); then
        echo "Missing, too old or not importable: ${need[*]}"
        # pip's default strategy upgrades only what these specs need.
        "$VENV_PYTHON" -m pip install "${need[@]}" || true
    fi
fi

# The venv is read by Apache's user and anyone running the CLI.
chmod -R a+rX "$VENV_DIR"
write_wrapper

probe "${REQUIREMENTS[@]}" > "$after_file" || true
if (( $(grep -c . "$after_file") != ${#REQUIREMENTS[@]} )); then
    echo "error: could not check the Python libraries with $VENV_PYTHON." >&2
    exit 1
fi
labels=""
[[ "$ACTION" == check ]] && labels="$before_file"
if ! report "$labels" "$after_file"; then
    echo "error: some Python libraries are still unusable (see above)." >&2
    exit 1
fi
if (( VENV_REBUILT )); then
    echo "The venv was (re)built: restart Apache so mod_wsgi loads it."
fi
