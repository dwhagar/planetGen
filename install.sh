#!/usr/bin/env bash
#
# install.sh
#
# One-shot Linux installer for deploying planetGen's web interface on an
# Apache2 server. `setup.py` stays scoped to the Python side only (the
# `stellarObjects` package plus the `sectorgen`/`systemgen` console
# scripts, installable on any OS); everything Linux/Apache-specific lives
# here instead:
#
#   1. Installs the Python package via `setup.py install`.
#   2. Pre-fetches the NLTK `words` corpus into a shared, world-readable
#      location (not a per-user home directory) so it works under any
#      user that later imports `stellarObjects` -- a login shell running
#      `sectorgen`/`systemgen`, or Apache's own locked-down `www-data`
#      running the `html/` CGI scripts. `stellarObjects/names.py` checks
#      `nltk.data.find()` before ever calling `download()`, so once this
#      step has populated a directory nltk's default search path already
#      covers, nothing later attempts a download of its own. See
#      `TODO.md`'s "Deployment bugs found in production" section for the
#      incident (`PermissionError: [Errno 13] ... '/var/www/nltk_data'`)
#      this fixes.
#   3. Makes the `html/` CGI scripts executable, independent of whatever
#      executable bit git happened to preserve on checkout (also see
#      `TODO.md` -- a `core.fileMode=false` git config on the authoring
#      machine silently dropped this once already, and nothing about a
#      git checkout should be trusted to carry it reliably).
#   4. Enables Apache's CGI module (`a2enmod cgid`).
#   5. Runs `apache/set-permissions.sh` to set ownership/permissions on
#      the deployed `html/`/`db/` directories for Apache's worker
#      user/group.
#   6. Prints the one remaining manual step: copying and enabling the
#      example virtual host config. This script never touches Apache's
#      site configuration itself -- ServerName, TLS, and logging are
#      site-specific decisions for a human to make, not something to
#      silently create or overwrite.
#
# Usage:
#   sudo ./install.sh
#
# Linux only (apt/a2enmod/systemd conventions) -- same scope as
# apache/set-permissions.sh, which this script calls.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HTML_DIR="$SCRIPT_DIR/html"
DB_DIR="$SCRIPT_DIR/db"
NLTK_DATA_DIR="${PLANETGEN_NLTK_DATA_DIR:-/usr/local/share/nltk_data}"

if [[ $EUID -ne 0 ]]; then
    echo "error: must be run as root, e.g.:" >&2
    echo "  sudo $0" >&2
    exit 1
fi

PYTHON="$(command -v python3 || command -v python || true)"
if [[ -z "$PYTHON" ]]; then
    echo "error: no python3/python found on PATH." >&2
    exit 1
fi

echo "== 1/5: Installing the Python package =="
# Also upgrades importlib_metadata, not just setuptools. Recent setuptools
# releases vendor their own newer importlib_metadata internally, but
# setuptools' own `extern` shim prefers a *real*, already-installed
# importlib_metadata package over that bundled copy when one is present on
# the path. On older distributions (e.g. Ubuntu 20.04's apt-provided
# python3-importlib-metadata, ~1.5.0) that old real package shadows the
# newer vendored one and is missing APIs (`EntryPoints`) current setuptools
# calls, crashing every `setup.py install` with `AttributeError: module
# 'importlib_metadata' has no attribute 'EntryPoints'`. Upgrading
# importlib_metadata alongside setuptools installs a compatible version
# ahead of the stale system one on the import path.
"$PYTHON" -m pip install --upgrade setuptools importlib_metadata
(cd "$SCRIPT_DIR" && "$PYTHON" setup.py install)

echo
echo "== 2/5: Fetching the NLTK 'words' corpus into $NLTK_DATA_DIR =="
mkdir -p "$NLTK_DATA_DIR"
"$PYTHON" -m nltk.downloader -d "$NLTK_DATA_DIR" words
chmod -R a+rX "$NLTK_DATA_DIR"

echo
echo "== 3/5: Making the CGI scripts and shell scripts executable =="
# No -maxdepth: every *.py under html/, at any subdirectory depth
# (html/lib/*.py included), needs this -- a previous version of this
# line was restricted to the top level only, which silently left
# html/lib/*.py non-executable/unreadable-as-intended after every
# install. apache/set-permissions.sh (below) re-does this same walk
# anyway with the correct final ownership, but doing it correctly here
# too means a plain `sudo ./install.sh` is never the reason this is wrong.
find "$HTML_DIR" -name '*.py' -exec chmod +x {} +
# Every *.sh anywhere in the repo (this script, update.sh,
# apache/set-permissions.sh, and any future one), not a hardcoded list --
# git checkouts made from a `core.fileMode=false` machine silently drop
# the executable bit on ANY file type, not just html/'s .py scripts (see
# TODO.md's "Deployment bugs found in production"), so a script added
# later doesn't need this list updated to be covered.
find "$SCRIPT_DIR" -name '*.sh' -exec chmod +x {} +

echo
echo "== 4/5: Enabling Apache's CGI module =="
if command -v a2enmod >/dev/null 2>&1; then
    a2enmod cgid
else
    echo "warning: a2enmod not found -- is apache2 installed?" >&2
    echo "  Try: sudo apt install apache2" >&2
fi

echo
echo "== 5/5: Setting directory ownership/permissions for Apache =="
"$SCRIPT_DIR/apache/set-permissions.sh" "$HTML_DIR" "$DB_DIR"

if [[ ! -f /etc/apache2/sites-available/planetgen.conf ]]; then
    cat <<EOF

------------------------------------------------------------------------
Install steps complete. One manual step remains -- create the Apache2
site from the example config (never done automatically, since
ServerName/TLS/logging are your call):

  1. Copy the example and edit it (at minimum, set ServerName):

       sudo cp "$SCRIPT_DIR/apache/planetgen.conf.example" /etc/apache2/sites-available/planetgen.conf
       sudo \${EDITOR:-nano} /etc/apache2/sites-available/planetgen.conf

  2. Enable the site and reload Apache:

       sudo a2ensite planetgen
       sudo systemctl reload apache2

See apache/README.md and html/README.md for more detail.
------------------------------------------------------------------------
EOF
else
    echo
    echo "Done. (/etc/apache2/sites-available/planetgen.conf already exists --"
    echo "not touching it; reload Apache yourself if this update needs it:"
    echo "  sudo systemctl reload apache2)"
fi
