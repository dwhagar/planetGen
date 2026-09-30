#!/usr/bin/env bash
#
# scripts/install-python-deps.sh
#
# Step 1 of install.sh: makes planetGen and the libraries it needs
# importable by the system Python, the one Apache's CGI scripts
# (`#!/usr/bin/env python3`) and mod_wsgi run under. Picks one of two paths
# and prints which one it took:
#
#   unmanaged  The interpreter lets pip install into it. Same pip install
#              install.sh has always done (see the comment on that branch
#              below for why each flag is there).
#
#   managed    The interpreter is "externally managed" (PEP 668): the
#              distribution ships an EXTERNALLY-MANAGED file in its stdlib
#              directory, as Debian 12+/Ubuntu 23.04+ do, and pip refuses
#              to install into it ("error: externally-managed-environment").
#              Here the libraries come from the distribution's own packages
#              instead (python3-flask, python3-nltk, ... via apt). Any
#              library apt has no package for, or only one older than the
#              floor setup.py asks for, is pip-installed into a venv
#              instead (PLANETGEN_VENV_DIR, default /opt/planetgen/venv),
#              which a .pth file then puts at the front of the system
#              Python's sys.path. planetGen itself isn't installed into
#              site-packages at all on this path: every entry point already
#              adds the checkout's src/ to sys.path itself (generate.py,
#              src/html/lib/apiclient.py, src/html/wsgi.py, and
#              src/migrateDb.py through its own sys.path[0]), so a small
#              /usr/local/bin/planetgen wrapper around generate.py stands in
#              for the console script pip would have made.
#
# Usage (as root; install.sh calls it this way):
#   scripts/install-python-deps.sh
#
# Environment:
#   PYTHON                 interpreter to install for (default: python3 on PATH)
#   PLANETGEN_PYTHON_MODE  auto (default), managed or unmanaged
#   PLANETGEN_VENV_DIR     fallback venv location (default /opt/planetgen/venv)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON:-$(command -v python3 || command -v python || true)}"
MODE="${PLANETGEN_PYTHON_MODE:-auto}"
VENV_DIR="${PLANETGEN_VENV_DIR:-/opt/planetgen/venv}"
PTH_NAME="planetgen-venv.pth"
WRAPPER=/usr/local/bin/planetgen

# Every runtime requirement install.sh needs: setup.py's install_requires
# plus its 'api' extra, as "<pip requirement> <apt package>". Keep in step
# with setup.py (src/tests/test_install_python_deps.py checks this).
REQUIREMENTS=(
    "nltk>=3.9.1 python3-nltk"
    "pymysql>=1.1.1 python3-pymysql"
    "dbutils>=3.1.0 python3-dbutils"
    "werkzeug>=3.0.0 python3-werkzeug"
    "rich>=13.7.0 python3-rich"
    "flask>=3.0.3 python3-flask"
    "flask-limiter>=3.7.0 python3-flask-limiter"
)

if [[ -z "$PYTHON" ]]; then
    echo "error: no python3/python found on PATH." >&2
    exit 1
fi

# The system site-packages directory a .pth file goes in, e.g.
# /usr/local/lib/python3.13/dist-packages on Debian/Ubuntu.
site_dir() {
    "$PYTHON" -c "import site; print(site.getsitepackages()[0])"
}

is_managed() {
    local stdlib
    stdlib="$("$PYTHON" -c "import sysconfig; print(sysconfig.get_path('stdlib'))")"
    [[ -f "$stdlib/EXTERNALLY-MANAGED" ]]
}

# Prints each requirement (from the arguments) the system Python doesn't
# satisfy: not installed, or older than its floor. Ignores the fallback
# venv (argument 1, may be empty), so this reports what the distribution
# itself provides.
unsatisfied() {
    "$PYTHON" - "$@" <<'EOF'
import re
import sys
from importlib import metadata

ignore = sys.argv[1]
path = [p for p in sys.path if not ignore or p != ignore]


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


for spec in sys.argv[2:]:
    name, floor = spec.split(">=")
    dists = list(metadata.distributions(name=name, path=path))
    if not dists and "-" in name:
        dists = list(metadata.distributions(name=name.replace("-", "_"), path=path))
    if not dists or key(dists[0].version) < key(floor):
        print(spec)
EOF
}

install_unmanaged() {
    echo "Python at $PYTHON is not externally managed: installing with pip."
    # `pip install .` (a proper, build-isolated PEP 517 install), NOT the
    # legacy `python3 setup.py install` this used to run. setuptools itself
    # now prints "Please avoid running setup.py directly" for that direct
    # invocation, and it's not just a style complaint: that legacy code path
    # is where two separate production incidents happened back to back (see
    # docs/TODO.md's "Deployment bugs found in production"). Both had the same
    # root cause -- setuptools' own vendoring shim (`extern`) prefers a
    # *real*, already-installed copy of a dependency it vendors
    # (`importlib_metadata`, then `packaging`) over its own newer bundled
    # copy whenever a real one is importable, so this system's old
    # apt-provided copies of each one in turn silently shadowed the working
    # vendored copy and crashed on a missing/changed API
    # (`importlib_metadata.EntryPoints`, then
    # `packaging.version.canonicalize_version`'s `strip_trailing_zero`
    # kwarg) -- and chasing each one individually with another `pip install
    # --upgrade <whatever's shadowed this time>` only fixes the specific
    # dependency that happened to break today, not the next one. `pip
    # install .`'s build isolation builds this package in a throwaway
    # environment that can't see this system's site-packages at all (only
    # the stdlib and pip's own freshly fetched build dependencies), so the
    # shadowing can't happen there regardless of which dependency it would
    # have hit -- avoiding this whole class of bug instead of patching it
    # dependency-by-dependency. This also means the global `setuptools`
    # system install no longer needs to be upgraded at all for this step,
    # which is one less thing on this box's system-wide Python environment
    # for this script to touch.
    #
    # --force-reinstall (not a plain `pip install .`): unlike the old `setup.py
    # install`, which unconditionally redid the install every run, plain `pip
    # install .` skips reinstalling when pip thinks the same version is
    # already installed -- true on every run between version bumps in
    # `stellarObjects/_version.py`. Since this script's whole point (via
    # update.sh) is redeploying whatever was just `git pull`-ed regardless of
    # whether the version string changed, skipping would silently leave the
    # previous run's installed copy in place, shadowing the freshly pulled
    # source the same way this whole section is otherwise about avoiding.
    #
    # The `api` extra (Flask/Flask-Limiter, see setup.py's `extras_require`)
    # is included here, not left to a separate manual `pip install .[api]`
    # some other doc might mention: every page under src/html/ is a thin
    # HTTP client over GET /api/... now (see html/lib/apiclient.py's own
    # docstring), so the web interface this script exists to deploy simply
    # doesn't work without it -- confirmed in production as
    # "ModuleNotFoundError: No module named 'flask'" from mod_wsgi once the
    # vhost's own handler-conflict and sys.path bugs (see wsgi.py) were fixed
    # and this became the next thing standing between a fresh install and a
    # working /api/search. A CLI-only use of this package (just `sectorgen`/
    # `systemgen`, no web interface ever deployed) wouldn't need it, but
    # nothing reaches this script without wanting the web interface.
    #
    # --ignore-installed: pulling in Flask this way surfaced a second,
    # unrelated production failure -- Flask 3.x needs blinker>=1.9.0, but
    # Ubuntu 22.04 ships blinker 1.4 pre-installed the old `distutils`
    # way (no `RECORD` file, so pip can't tell which files are its to
    # remove). `--force-reinstall`/`--upgrade` both still need to *upgrade*
    # it, which means uninstalling that old copy first, which fails with
    # "Cannot uninstall blinker 1.4 ... distutils installed project" and
    # aborts the whole install before it ever reaches flask/flask-limiter/
    # planetGen itself (see pip's own install order in its output -- it
    # aborts alphabetically-ish partway through, well before the packages
    # that actually matter here). `--ignore-installed` sidesteps the
    # uninstall step entirely: pip just installs its own copy into
    # /usr/local's site-packages, which already comes before apt's
    # /usr/lib/python3/dist-packages on sys.path, so the newer pip-managed
    # copy shadows the old system one without ever touching it -- the
    # standard workaround for this well-known Debian/Ubuntu packaging class
    # of error, not specific to blinker (a future dependency bump could hit
    # the same wall with some other apt-provided package).
    "$PYTHON" -m pip install --upgrade pip
    "$PYTHON" -m pip install --upgrade --force-reinstall --ignore-installed "${SCRIPT_DIR}[api]"
    echo "Python install path: pip (unmanaged interpreter)."
}

install_managed() {
    echo "Python at $PYTHON is externally managed (PEP 668): using distribution packages."

    local specs=() packages=() line
    for line in "${REQUIREMENTS[@]}"; do
        specs+=("${line%% *}")
        packages+=("${line##* }")
    done

    local have_apt=0
    if command -v apt-get >/dev/null 2>&1; then
        have_apt=1
        apt-get update -qq
        # Only packages this distribution actually has: one it doesn't
        # would make apt-get fail the whole batch, and is exactly what the
        # venv below is for.
        local available=() missing=() pkg candidate
        for pkg in "${packages[@]}"; do
            candidate="$(apt-cache policy "$pkg" 2>/dev/null | awk '/Candidate:/ {print $2}')"
            if [[ -n "$candidate" && "$candidate" != "(none)" ]]; then
                available+=("$pkg")
            else
                missing+=("$pkg")
            fi
        done
        if (( ${#available[@]} )); then
            DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends "${available[@]}"
        fi
        if (( ${#missing[@]} )); then
            echo "No distribution package for: ${missing[*]}"
        fi
    else
        echo "warning: apt-get not found -- can't install distribution packages;" >&2
        echo "  anything missing goes into the venv at $VENV_DIR instead." >&2
    fi

    local pth venv_site=""
    pth="$(site_dir)/$PTH_NAME"
    if [[ -f "$VENV_DIR/bin/python" ]]; then
        venv_site="$("$VENV_DIR/bin/python" -c "import sysconfig; print(sysconfig.get_path('purelib'))" 2>/dev/null || true)"
    fi

    local need=()
    mapfile -t need < <(unsatisfied "$venv_site" "${specs[@]}")

    if (( ${#need[@]} == 0 )); then
        # Everything comes from the distribution now; a venv left over from
        # an earlier run would only shadow it.
        if [[ -f "$pth" ]]; then
            rm -f "$pth"
            echo "Removed $pth (the distribution now provides everything)."
        fi
        echo "Python install path: distribution packages (apt)."
    else
        echo "Not provided by the distribution (missing or too old): ${need[*]}"
        if (( have_apt )) && ! "$PYTHON" -c "import ensurepip" >/dev/null 2>&1; then
            DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3-venv
        fi
        # --clear rebuilds it every run (like the unmanaged path's
        # --force-reinstall), which also keeps it working after a
        # distribution upgrade changes the Python version under it.
        # --system-site-packages lets pip see what apt already installed,
        # so only the requirements listed above (and whatever they need
        # that apt doesn't have) go in here.
        mkdir -p "$(dirname "$VENV_DIR")"
        "$PYTHON" -m venv --clear --system-site-packages "$VENV_DIR"
        "$VENV_DIR/bin/python" -m pip install --upgrade pip
        "$VENV_DIR/bin/python" -m pip install --upgrade "${need[@]}"
        venv_site="$("$VENV_DIR/bin/python" -c "import sysconfig; print(sysconfig.get_path('purelib'))")"
        chmod -R a+rX "$VENV_DIR"
        # An `import` line rather than a bare path, which site.py would
        # append after the distribution's own dist-packages: a package apt
        # has but at a version below our floor must lose to the newer copy
        # here, not win over it.
        mkdir -p "$(dirname "$pth")"
        printf 'import sys; p = %s; p in sys.path or sys.path.insert(0, p)\n' \
            "$("$PYTHON" -c 'import sys; print(repr(sys.argv[1]))' "$venv_site")" > "$pth"
        echo "Wrote $pth, putting $venv_site first on the system Python's sys.path."
        echo "Python install path: distribution packages (apt) plus a venv at $VENV_DIR for: ${need[*]}"
    fi

    # Stands in for the `planetgen` console script pip would have written.
    cat > "$WRAPPER" <<EOF
#!/bin/sh
# Written by planetGen's scripts/install-python-deps.sh (managed Python).
exec "$PYTHON" "$SCRIPT_DIR/generate.py" "\$@"
EOF
    chmod 755 "$WRAPPER"
}

case "$MODE" in
    auto)
        if is_managed; then install_managed; else install_unmanaged; fi ;;
    managed)
        install_managed ;;
    unmanaged)
        install_unmanaged ;;
    *)
        echo "error: PLANETGEN_PYTHON_MODE must be auto, managed or unmanaged (got '$MODE')." >&2
        exit 1 ;;
esac
