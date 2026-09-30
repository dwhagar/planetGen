#!/usr/bin/env bash
#
# scripts/install-python-deps.sh
#
# Step 1 of install.sh: makes planetGen and the libraries it needs
# importable by the system Python, the one mod_wsgi (and the CLI tools)
# run under. Picks one of two paths
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
# Usage (as root):
#   scripts/install-python-deps.sh           install.sh: full install
#   scripts/install-python-deps.sh --check   update.sh: install nothing that
#                                            is already there
#
# --check imports every requirement with the interpreter the site runs
# under (so the venv .pth, if any, is in effect) and compares its version
# with the floor. Only a requirement that is missing, below its floor or
# fails to import gets installed, the same way the full install would
# have on this host: apt first and then the venv on a managed Python, pip
# on an unmanaged one. Nothing already present is reinstalled or
# rebuilt. It prints one line per requirement (present, installed, upgraded,
# repaired or failed) and exits non-zero if anything is still unusable. It also puts
# back the /usr/local/bin/planetgen wrapper if it is missing or stale.
#
# Environment:
#   PYTHON                 interpreter to install for (default: python3 on PATH)
#   PLANETGEN_PYTHON_MODE  auto (default), managed or unmanaged
#   PLANETGEN_VENV_DIR     fallback venv location (default /opt/planetgen/venv)

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
MODE="${PLANETGEN_PYTHON_MODE:-auto}"
VENV_DIR="${PLANETGEN_VENV_DIR:-/opt/planetgen/venv}"
PTH_NAME="planetgen-venv.pth"
WRAPPER=/usr/local/bin/planetgen

# Every runtime requirement install.sh needs: setup.py's install_requires
# plus its 'api' extra, as "<pip requirement> <apt package>". Keep in step
# with setup.py (src/tests/test_install_python_deps.py checks this).
# --check imports each one by its pip name with "-" turned into "_"
# (flask-limiter -> flask_limiter).
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

# Prints "<state> <spec> <detail>" for each requirement (from the
# arguments) as the system Python sees it, venv .pth and all: state is ok,
# missing (not installed), old (below its floor; detail is the version
# found) or broken (installed but the import raised; detail is the
# error). Imports each one for real, since an installed distribution
# whose own dependencies are missing is not usable either.
probe() {
    "$PYTHON" - "$@" <<'EOF'
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
        try:
            version = metadata.version(name.replace("-", "_"))
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

# Sets VENV_SITE to the fallback venv's site-packages, (re)creating the
# venv first if needed: with "clear" always (the full install rebuilds
# it every run), otherwise only when it's missing or its interpreter no
# longer runs (a distribution upgrade removed the Python it was made
# from). Then pip-installs the given requirements into it and writes the
# .pth that puts it first on the system Python's sys.path.
venv_install() {
    local how="$1"; shift
    if [[ "$how" == clear ]] || ! "$VENV_DIR/bin/python" -c "" >/dev/null 2>&1; then
        if command -v apt-get >/dev/null 2>&1 && ! "$PYTHON" -c "import ensurepip" >/dev/null 2>&1; then
            DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3-venv
        fi
        # --clear also keeps it working after a distribution upgrade
        # changes the Python version under it. --system-site-packages
        # lets pip see what apt already installed, so only the given
        # requirements (and whatever they need that apt doesn't have) go
        # in here.
        mkdir -p "$(dirname "$VENV_DIR")"
        "$PYTHON" -m venv --clear --system-site-packages "$VENV_DIR"
        "$VENV_DIR/bin/python" -m pip install --upgrade pip
    fi
    "$VENV_DIR/bin/python" -m pip install --upgrade "$@"
    VENV_SITE="$("$VENV_DIR/bin/python" -c "import sysconfig; print(sysconfig.get_path('purelib'))")"
    chmod -R a+rX "$VENV_DIR"
    # An `import` line rather than a bare path, which site.py would
    # append after the distribution's own dist-packages: a package apt
    # has but at a version below our floor must lose to the newer copy
    # here, not win over it.
    local pth
    pth="$(site_dir)/$PTH_NAME"
    mkdir -p "$(dirname "$pth")"
    printf 'import sys; p = %s; p in sys.path or sys.path.insert(0, p)\n' \
        "$("$PYTHON" -c 'import sys; print(repr(sys.argv[1]))' "$VENV_SITE")" > "$pth"
    echo "Wrote $pth, putting $VENV_SITE first on the system Python's sys.path."
}

# Stands in for the `planetgen` console script pip would make, but runs
# this checkout's generate.py, so the CLI always runs the code that was
# just pulled (like the web app, which imports from the checkout too)
# without anything being reinstalled. Rewritten only when it differs.
write_wrapper() {
    local want
    want="$(cat <<EOF
#!/bin/sh
# Written by planetGen's scripts/install-python-deps.sh.
exec "$PYTHON" "$SCRIPT_DIR/generate.py" "\$@"
EOF
)"
    if [[ "$(cat "$WRAPPER" 2>/dev/null || true)" != "$want" ]]; then
        printf '%s\n' "$want" > "$WRAPPER"
        echo "Wrote $WRAPPER (runs $SCRIPT_DIR/generate.py)."
    fi
    chmod 755 "$WRAPPER"
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
    # --force-reinstall (not a plain `pip install .`): plain `pip install .`
    # skips reinstalling when pip thinks the same version is already
    # installed -- true on every run between version bumps in
    # `stellarObjects/_version.py` -- and install.sh is the full reinstall.
    # (update.sh never comes here: it runs --check, which installs only
    # what's missing. Nothing needs the pip-installed copy of planetGen to
    # be current anyway, since every entry point imports from the checkout
    # and write_wrapper below replaces pip's console script.)
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
    # Replaces pip's own console script, which would run the copy pip
    # just installed and go stale after the next update.sh.
    write_wrapper
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
        venv_install clear "${need[@]}"
        echo "Python install path: distribution packages (apt) plus a venv at $VENV_DIR for: ${need[*]}"
    fi

    write_wrapper
}

# update.sh's path: install nothing that is already usable. See the
# header for what it does; $1 is "managed" or "unmanaged".
check_requirements() {
    local mode="$1" line
    local specs=() packages=()
    declare -A package_of=()
    for line in "${REQUIREMENTS[@]}"; do
        specs+=("${line%% *}")
        packages+=("${line##* }")
        package_of["${line%% *}"]="${line##* }"
    done

    echo "Checking the Python libraries with $PYTHON ($mode Python)."
    local before=() need=() state spec detail
    mapfile -t before < <(probe "${specs[@]}")
    for line in "${before[@]}"; do
        read -r state spec detail <<< "$line"
        [[ "$state" == ok ]] || need+=("$spec")
    done

    if (( ${#need[@]} )); then
        echo "Missing, too old or not importable: ${need[*]}"
        if [[ "$mode" == unmanaged ]]; then
            # A plain install upgrades only what the specs need. Falls back
            # to install.sh's --ignore-installed for the case it exists for:
            # an old distutils-installed dependency pip can't uninstall.
            "$PYTHON" -m pip install --upgrade "${need[@]}" \
                || "$PYTHON" -m pip install --upgrade --ignore-installed "${need[@]}" \
                || true
        else
            local wanted=() pkg candidate
            if command -v apt-get >/dev/null 2>&1; then
                apt-get update -qq || true
                for spec in "${need[@]}"; do
                    pkg="${package_of[$spec]}"
                    candidate="$(apt-cache policy "$pkg" 2>/dev/null | awk '/Candidate:/ {print $2}')"
                    [[ -n "$candidate" && "$candidate" != "(none)" ]] && wanted+=("$pkg")
                done
                if (( ${#wanted[@]} )); then
                    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends "${wanted[@]}" || true
                fi
            fi
            # Whatever apt couldn't satisfy goes into the venv, which is
            # kept (not rebuilt) if it still works.
            local still=()
            mapfile -t still < <(probe "${need[@]}" | awk '$1 != "ok" {print $2}')
            if (( ${#still[@]} )); then
                venv_install keep "${still[@]}" || true
            fi
        fi
    fi

    # Report against the state before, so each line says what happened.
    local after=() failed=0 old_state
    mapfile -t after < <(probe "${specs[@]}")
    if (( ${#after[@]} != ${#specs[@]} )); then
        echo "error: could not check the Python libraries with $PYTHON." >&2
        return 1
    fi
    local i
    for i in "${!after[@]}"; do
        read -r state spec detail <<< "${after[$i]}"
        read -r old_state _ _ <<< "${before[$i]}"
        if [[ "$state" != ok ]]; then
            printf '  %-10s %s (%s: %s)\n' failed "$spec" "$state" "$detail"
            failed=1
        elif [[ "$old_state" == ok ]]; then
            printf '  %-10s %s %s\n' present "$spec" "$detail"
        elif [[ "$old_state" == missing ]]; then
            printf '  %-10s %s %s\n' installed "$spec" "$detail"
        elif [[ "$old_state" == old ]]; then
            printf '  %-10s %s %s\n' upgraded "$spec" "$detail"
        else
            printf '  %-10s %s %s\n' repaired "$spec" "$detail"
        fi
    done

    write_wrapper
    if (( failed )); then
        echo "error: some Python libraries are still unusable (see above)." >&2
        echo "  Running sudo ./install.sh does a full reinstall." >&2
        return 1
    fi
}

case "$MODE" in
    auto)
        if is_managed; then MODE=managed; else MODE=unmanaged; fi ;;
    managed|unmanaged) ;;
    *)
        echo "error: PLANETGEN_PYTHON_MODE must be auto, managed or unmanaged (got '$MODE')." >&2
        exit 1 ;;
esac

if [[ "$ACTION" == check ]]; then
    check_requirements "$MODE"
elif [[ "$MODE" == managed ]]; then
    install_managed
else
    install_unmanaged
fi
