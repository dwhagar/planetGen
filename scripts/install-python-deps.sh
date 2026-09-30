#!/usr/bin/env bash
# TODO(installers #50): macOS ships bash 3.2, where declare -A and mapfile
# below fail; rewrite them (or require Homebrew bash) and use Homebrew or a
# venv instead of apt on macOS. See docs/TODO.md item 50.
#
# scripts/install-python-deps.sh
#
# Makes planetGen's libraries importable by the system Python, the one
# mod_wsgi (and the CLI tools) run under. Everything goes into that
# Python's own system-wide site-packages; there is no virtual
# environment. Picks one of two paths and prints which one it took:
#
#   unmanaged  The interpreter lets pip install into it. Same pip install
#              install.sh has always done (see the comment on that branch
#              below for why each flag is there).
#
#   managed    The interpreter is "externally managed" (PEP 668): the
#              distribution ships an EXTERNALLY-MANAGED file in its stdlib
#              directory, as Debian 12+/Ubuntu 23.04+ do. Every library the
#              distribution packages at or above setup.py's floor comes from
#              apt (python3-flask, python3-nltk, ...). Only one apt has no
#              package for, or ships below the floor, is pip-installed
#              system-wide (--break-system-packages) into
#              /usr/local/lib/python3.X/dist-packages, which comes before
#              apt's /usr/lib/python3/dist-packages on sys.path, so it
#              wins without apt's copy being touched. The report says
#              which libraries came from pip and why.
#
#              Earlier versions put that fallback in a venv
#              (/opt/planetgen/venv, or PLANETGEN_VENV_DIR) with a
#              planetgen-venv.pth file; both are removed here when found.
#
# Whatever pip installs on either path comes from requirements.lock: the
# exact versions scripts/lock-requirements.sh pinned, checked against the
# lock's sha256 hashes (--require-hashes), so nothing newer or altered on
# PyPI gets in unnoticed. apt packages are apt's to verify.
#
# planetGen itself is never installed into site-packages: every entry
# point adds the checkout's src/ to sys.path itself (generate.py,
# src/html/lib/apiclient.py, src/html/wsgi.py, and src/migrateDb.py
# through its own sys.path[0]), and /usr/local/bin/planetgen is a small
# wrapper around the checkout's generate.py. (The unmanaged path still
# pip-installs the package, as it always has, but nothing runs that copy.)
#
# Usage (as root):
#   scripts/install-python-deps.sh           install.sh: full install
#   scripts/install-python-deps.sh --check   update.sh: install nothing that
#                                            is already there
#
# --check imports every requirement with the system Python and compares
# its version with the floor. Only a requirement that is missing, below
# its floor or fails to import gets installed, the same way the full
# install would have on this host (apt first, then system-wide pip on a
# managed Python; pip on an unmanaged one). Nothing already present is
# reinstalled. Both actions print one line per requirement with where it
# came from, and exit non-zero if anything is still unusable. Both also
# put back the /usr/local/bin/planetgen wrapper if it's missing or stale.
#
# Environment:
#   PYTHON                 interpreter to install for (default: python3 on PATH)
#   PLANETGEN_PYTHON_MODE  auto (default), managed or unmanaged
#   PLANETGEN_VENV_DIR     old fallback venv to remove (default /opt/planetgen/venv)

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
LEGACY_VENV_DIR="${PLANETGEN_VENV_DIR:-/opt/planetgen/venv}"
LEGACY_PTH_NAME="planetgen-venv.pth"
WRAPPER=/usr/local/bin/planetgen
# Exact versions and file hashes for everything pip installs
# (scripts/lock-requirements.sh writes it).
LOCK="$SCRIPT_DIR/requirements.lock"
LOCK_PINS="$SCRIPT_DIR/scripts/lock_pins.py"

# Every runtime requirement install.sh needs: setup.py's install_requires
# plus its 'api' extra, as "<pip requirement> <apt package>". Keep in step
# with setup.py (src/tests/test_install_python_deps.py checks this).
# Each is imported by its pip name with "-" turned into "_"
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

SPECS=()
declare -A PACKAGE_OF=()
for line in "${REQUIREMENTS[@]}"; do
    SPECS+=("${line%% *}")
    PACKAGE_OF["${line%% *}"]="${line##* }"
done

is_managed() {
    local stdlib
    stdlib="$("$PYTHON" -c "import sysconfig; print(sysconfig.get_path('stdlib'))")"
    [[ -f "$stdlib/EXTERNALLY-MANAGED" ]]
}

# Prints "<state> <spec> <version or error> <directory>" for each
# requirement (from the arguments) as the system Python sees it: state is
# ok, missing (not installed), old (below its floor) or broken (installed
# but the import raised). Imports each one for real, since an installed
# distribution whose own dependencies are missing is not usable either.
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
    dist = None
    for candidate in (name, name.replace("-", "_")):
        try:
            dist = metadata.distribution(candidate)
            break
        except metadata.PackageNotFoundError:
            pass
    if dist is None:
        print("missing", spec, "-", "-")
        continue
    where = str(dist.locate_file("")).rstrip("/") or "-"
    if key(dist.version) < key(floor):
        print("old", spec, dist.version, where)
        continue
    try:
        importlib.import_module(name.replace("-", "_"))
    except Exception as exc:  # any import failure makes it unusable
        print("broken", spec, f"{type(exc).__name__}:{exc}".replace(" ", "_").replace("\n", "_"), where)
        continue
    print("ok", spec, dist.version, where)
EOF
}

# The version apt would install for a package, or nothing if it has none.
apt_candidate() {
    command -v apt-cache >/dev/null 2>&1 || return 0
    local candidate
    candidate="$(apt-cache policy "$1" 2>/dev/null | awk '/Candidate:/ {print $2}')"
    [[ "$candidate" == "(none)" ]] || printf '%s' "$candidate"
}

# apt-installs the distribution package of each given requirement that
# apt has at all (one it doesn't would make apt-get fail the whole
# batch). Packages already installed are left alone by apt-get.
apt_install() {
    command -v apt-get >/dev/null 2>&1 || {
        echo "warning: apt-get not found -- installing with pip instead." >&2
        return 0
    }
    apt-get update -qq || true
    local spec pkg wanted=() lacking=()
    for spec in "$@"; do
        pkg="${PACKAGE_OF[$spec]}"
        if [[ -n "$(apt_candidate "$pkg")" ]]; then wanted+=("$pkg"); else lacking+=("$pkg"); fi
    done
    if (( ${#wanted[@]} )); then
        DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends "${wanted[@]}" || true
    fi
    if (( ${#lacking[@]} )); then
        echo "No distribution package for: ${lacking[*]}"
    fi
}

# pip-installs the given requirements into the system Python's own
# site-packages (/usr/local/lib/python3.X/dist-packages on Debian/Ubuntu),
# never removing anything. A plain `pip install --upgrade` would
# uninstall the older copy first, and for a library apt installed that
# deletes apt's own files (it does for egg-info packages such as
# python3-pymysql), leaving dpkg's view of the system broken. So pip
# resolves first (--dry-run --report), and then exactly the distributions
# it picked -- the requirements plus any dependency they need newer than
# what's installed -- go into /usr/local with --ignore-installed
# --no-deps. /usr/local comes before /usr/lib/python3/dist-packages on
# sys.path, so each shadows apt's copy without touching it, and
# everything else apt provides stays apt's. A managed Python also needs
# --break-system-packages for pip to write there at all.
#
# scripts/lock_pins.py does the resolving, holding everything pip adds to
# the version in requirements.lock, and the install uses --require-hashes,
# so pip only accepts the exact files the lock names.
pip_install_system() {
    if ! "$PYTHON" -m pip --version >/dev/null 2>&1 && command -v apt-get >/dev/null 2>&1; then
        DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3-pip || true
    fi
    local flags=() hashed status=0 pins=()
    [[ "$MODE" == managed ]] && flags+=(--break-system-packages)
    hashed="$(mktemp)"
    "$PYTHON" -I "$LOCK_PINS" resolve "$LOCK" "$hashed" ${flags[@]+"${flags[@]}"} -- "$@" || status=$?
    if (( status != 3 )); then
        (( status == 0 )) && mapfile -t pins < <(awk '{print $1}' "$hashed")
        if (( ${#pins[@]} )); then
            # An older copy pip itself installed (outside apt's
            # directory) is pip's to replace, and would otherwise leave its
            # metadata behind next to the new one.
            local owned=()
            mapfile -t owned < <("$PYTHON" - "${pins[@]%%==*}" <<'EOF'
import sys
from importlib import metadata

for name in sys.argv[1:]:
    try:
        where = str(metadata.distribution(name).locate_file(""))
    except metadata.PackageNotFoundError:
        continue
    if not where.startswith("/usr/lib/python3/"):
        print(name)
EOF
            )
            if (( ${#owned[@]} )); then
                "$PYTHON" -m pip uninstall -y "${flags[@]}" "${owned[@]}" || true
            fi
            echo "Installing with pip into the system Python: ${pins[*]}"
            "$PYTHON" -m pip install "${flags[@]}" --ignore-installed --no-deps --require-hashes -r "$hashed" || true
        fi
    else
        # pip older than 22.2 has no --report: the same --ignore-installed
        # the unmanaged install uses, dependencies and all, held to the
        # locked versions (though not their hashes, which would need
        # every dependency listed up front).
        "$PYTHON" -I "$LOCK_PINS" constraints "$LOCK" > "$hashed"
        "$PYTHON" -m pip install "${flags[@]}" --ignore-installed -c "$hashed" "$@" || true
    fi
    rm -f "$hashed"
}

# Earlier versions put the fallback libraries in a venv that a .pth file
# put first on sys.path. Removed before anything is checked, so the
# check sees the system Python on its own and fills any gap system-wide.
remove_legacy_venv() {
    local pth
    for pth in "$("$PYTHON" -c "import site; print(site.getsitepackages()[0])")/$LEGACY_PTH_NAME" \
               /usr/local/lib/python3*/dist-packages/"$LEGACY_PTH_NAME"; do
        if [[ -f "$pth" ]]; then
            rm -f "$pth"
            echo "Removed $pth (the old venv fallback; everything is system-wide now)."
        fi
    done
    if [[ -f "$LEGACY_VENV_DIR/pyvenv.cfg" ]]; then
        rm -rf "$LEGACY_VENV_DIR"
        echo "Removed the old fallback venv at $LEGACY_VENV_DIR."
        rmdir "$(dirname "$LEGACY_VENV_DIR")" 2>/dev/null || true
    fi
}

# Where a library came from, for the report: apt, or pip and why apt
# couldn't provide it.
source_of() {
    local spec="$1" where="$2" pkg candidate
    case "$where" in
        /usr/lib/python3/dist-packages|/usr/lib/python3/dist-packages/*) echo "apt"; return ;;
        -) echo "-"; return ;;
    esac
    if [[ "$MODE" != managed ]]; then
        echo "pip"; return
    fi
    pkg="${PACKAGE_OF[$spec]}"
    candidate="$(apt_candidate "$pkg")"
    if [[ -z "$candidate" ]]; then
        echo "pip, apt has no $pkg"
    elif dpkg --compare-versions "$candidate" lt "${spec#*>=}" 2>/dev/null; then
        echo "pip, apt's $pkg is $candidate, below ${spec#*>=}"
    else
        echo "pip, as a newer dependency of another library than apt's $candidate"
    fi
}

# Prints one line per requirement from the probe lines on stdin; with
# the probe lines from before any install in $1 (a file), labels each as
# present/installed/upgraded/repaired instead of ok. Returns 1 if any
# requirement is still unusable.
report() {
    local before_file="${1:-}" state spec detail where old failed=0 label i=0
    local -a before=()
    [[ -n "$before_file" ]] && mapfile -t before < "$before_file"
    while read -r state spec detail where; do
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
            printf '  %-10s %-24s %-8s (%s)\n' "$label" "$spec" "$detail" "$(source_of "$spec" "$where")"
        fi
        i=$((i + 1))
    done
    return "$failed"
}

# Probes every requirement and reports; fails if any is unusable or the
# probe itself didn't run.
final_report() {
    local after
    after="$(probe "${SPECS[@]}")" || true
    if (( $(grep -c . <<< "$after") != ${#SPECS[@]} )); then
        echo "error: could not check the Python libraries with $PYTHON." >&2
        return 1
    fi
    if ! report "${1:-}" <<< "$after"; then
        echo "error: some Python libraries are still unusable (see above)." >&2
        return 1
    fi
}

# Just the specs from probe lines on stdin that aren't ok.
not_ok() {
    awk '$1 != "ok" {print $2}'
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
    #
    # The libraries come from requirements.lock with --require-hashes, so
    # pip installs exactly the locked files, and planetGen itself after
    # them with --no-deps, since a local directory has no hash to check.
    "$PYTHON" -m pip install --upgrade pip
    "$PYTHON" -m pip install --ignore-installed --require-hashes -r "$LOCK"
    "$PYTHON" -m pip install --upgrade --force-reinstall --ignore-installed --no-deps "$SCRIPT_DIR"
    # Replaces pip's own console script, which would run the copy pip
    # just installed and go stale after the next update.sh.
    write_wrapper
    echo "Python install path: pip (unmanaged interpreter)."
    final_report
}

install_managed() {
    echo "Python at $PYTHON is externally managed (PEP 668): using distribution packages."
    remove_legacy_venv
    apt_install "${SPECS[@]}"
    local need=()
    mapfile -t need < <(probe "${SPECS[@]}" | not_ok)
    if (( ${#need[@]} )); then
        echo "Not provided well enough by apt (missing, too old or broken): ${need[*]}"
        echo "Installing those system-wide with pip."
        pip_install_system "${need[@]}"
    fi
    write_wrapper
    final_report
}

# update.sh's path: install nothing that is already usable.
check_requirements() {
    echo "Checking the Python libraries with $PYTHON ($MODE Python)."
    [[ "$MODE" == managed ]] && remove_legacy_venv
    local before_file need=()
    before_file="$(mktemp)"
    probe "${SPECS[@]}" > "$before_file"
    mapfile -t need < <(not_ok < "$before_file")

    if (( ${#need[@]} )); then
        echo "Missing, too old or not importable: ${need[*]}"
        if [[ "$MODE" == managed ]]; then
            apt_install "${need[@]}"
            mapfile -t need < <(probe "${need[@]}" | not_ok)
            if (( ${#need[@]} )); then
                echo "apt can't provide: ${need[*]} -- installing those system-wide with pip."
            fi
        fi
        if (( ${#need[@]} )); then
            pip_install_system "${need[@]}"
        fi
    fi

    write_wrapper
    if ! final_report "$before_file"; then
        rm -f "$before_file"
        echo "  Running sudo ./install.sh does a full reinstall." >&2
        return 1
    fi
    rm -f "$before_file"
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
    check_requirements
elif [[ "$MODE" == managed ]]; then
    install_managed
else
    install_unmanaged
fi
