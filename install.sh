#!/usr/bin/env bash
#
# install.sh
#
# One-shot installer for deploying planetGen's web interface: on Linux
# behind an Apache2 server with mod_wsgi (docs/deployment/apache.md), and
# on macOS under gunicorn and launchd with Homebrew's nginx in front
# (docs/deployment/macos.md). install.ps1 is the Windows counterpart;
# the three follow the same steps, so a change to one belongs in all. `setup.py` stays scoped to the Python side only (the
# `stellarObjects` package plus the `sectorgen`/`systemgen` console
# scripts, installable on any OS); everything Linux/Apache-specific lives
# here instead:
#
#   1. Runs `scripts/install-python-deps.sh` to install the Python package
#      and its libraries: a build-isolated `pip install` on an ordinary
#      Python, or distribution (apt) packages on an externally managed one
#      (PEP 668, e.g. Ubuntu 24.04+), with system-wide pip only for
#      libraries the distribution lacks or ships too old. No venv.
#   2. Runs `planetgen.cli.migrate` (with a progress bar; when a migration is
#      pending it first asks, y/N with a 30-second timeout defaulting to
#      N, whether to delete the galaxy data instead) against the
#      configured MySQL database
#      ($PLANETGEN_MYSQL_* in this shell's environment, or the vhost's
#      `SetEnv` directives once deployed), bringing it up to the current
#      schema (`planetgen/db/schema.sql`) if it isn't already. A no-op
#      for a database that's already current. Needs step 1 done first,
#      since it imports `stellarObjects`.
#   3. Pre-fetches (unless it's already there) the NLTK `words` corpus into a shared, world-readable
#      location (not a per-user home directory) so it works under any
#      user that later imports `stellarObjects` -- a login shell running
#      `sectorgen`/`systemgen`, or Apache's own locked-down `www-data`
#      running the `src/html/` web app. `planetgen/names/wordlists.py` checks
#      `nltk.data.find()` before ever calling `download()`, so once this
#      step has populated a directory nltk's default search path already
#      covers, nothing later attempts a download of its own. See
#      `docs/TODO.md`'s "Deployment bugs found in production" section for the
#      incident (`PermissionError: [Errno 13] ... '/var/www/nltk_data'`)
#      this fixes. Then, unless --skip-database, it offers (y/N, 30
#      seconds, default N; skipped with no terminal) to run the
#      population pass, `generate.py population`; POPULATION=1 runs it
#      without asking (offer_population_pass in scripts/deploy-common.sh).
#   4. Makes the repo's shell scripts (and `src/html/`'s Python files) executable, independent of whatever
#      executable bit git happened to preserve on checkout (also see
#      `docs/TODO.md` -- a `core.fileMode=false` git config on the authoring
#      machine silently dropped this once already, and nothing about a
#      git checkout should be trusted to carry it reliably).
#   5. Enables Apache's headers, deflate and wsgi modules, installing
#      mod_wsgi (libapache2-mod-wsgi-py3) first if apt can and it's
#      missing. On macOS, installs and starts the gunicorn launchd daemon
#      (examples/macos/org.planetgen.gunicorn.plist) instead. Steps 3 and
#      5 come from scripts/deploy-common.sh, which update.sh shares.
#   6. Runs `examples/apache/set-permissions.sh` to set ownership/permissions on
#      the deployed `src/html/`/`db/` directories for Apache's worker
#      user/group, and config.json to root:<apache group>, mode 640.
#   7. Runs `examples/apache/create-cache-dir.sh` to create the Galaxy Map's
#      on-disk tile cache (`tile_cache.dir` in config.json, default
#      /var/cache/planetgen/tiles) owned by Apache's worker user.
#   8. Runs `examples/apache/setup-debug-log.sh` to set up the debug log
#      (`log_file` in config.json, default /var/log/planetgen.log) whether
#      or not `debug` is on, and the activity log's folder (`log_dir`),
#      mode 0660 for Apache's user and group (CLI users must be in that
#      group to append), checks Apache's user can write both, and
#      installs their rotation. A log it can't set up only warns, with the
#      commands that fix it; it never stops the install (OPS.5).
#   9. Prints the one remaining manual step: copying and enabling the
#      example virtual host config. This script never touches Apache's
#      site configuration itself -- ServerName, TLS, and logging are
#      site-specific decisions for a human to make, not something to
#      silently create or overwrite.
#
# Usage:
#   sudo ./install.sh [--skip-database]
#   sudo POPULATION=1 ./install.sh   (also run the population pass)
#
# --skip-database leaves out step 2, for a host whose database isn't set
# up yet (and for CI); run `sudo ./update.sh` once it is.
#
# On macOS, step 1 makes a virtual environment from Homebrew's python3
# (brew install python@3.12 first; Homebrew itself won't run as root, so
# this never calls it), steps 6-8 use _www, and the debug log rotates
# with newsyslog. Runs under macOS's bash 3.2.

set -euo pipefail

SKIP_DATABASE=0
case "${1:-}" in
    "") ;;
    --skip-database) SKIP_DATABASE=1 ;;
    *)
        echo "error: unknown argument '$1' (expected nothing or --skip-database)." >&2
        exit 1 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HTML_DIR="$SCRIPT_DIR/src/html"
DB_DIR="$SCRIPT_DIR/db"
NLTK_DATA_DIR="${PLANETGEN_NLTK_DATA_DIR:-/usr/local/share/nltk_data}"

if [[ $EUID -ne 0 ]]; then
    echo "error: must be run as root, e.g.:" >&2
    echo "  sudo $0" >&2
    exit 1
fi

# shellcheck source=scripts/deploy-common.sh
source "$SCRIPT_DIR/scripts/deploy-common.sh"

PYTHON="${PYTHON:-$(default_python)}"
if [[ -z "$PYTHON" ]]; then
    echo "error: no python3/python found on PATH." >&2
    exit 1
fi

echo "== 1/8: Installing the Python package and its libraries =="
# pip on an ordinary Python; distribution packages (plus system-wide pip
# for anything the distribution lacks or ships too old) on an externally
# managed one (PEP 668, e.g. Ubuntu 24.04+). See that script for
# the details of each path; it prints which one it took.
PYTHON="$PYTHON" bash "$SCRIPT_DIR/scripts/install-python-deps.sh"
use_site_python

echo
echo "== 2/8: Migrating the configured MySQL database to the current schema =="
if (( SKIP_DATABASE )); then
    echo "Skipped (--skip-database). Run sudo ./update.sh once the database is set up."
else
    migrate_or_reset_db
fi

echo
echo "== 3/8: Fetching the NLTK 'words' corpus into $NLTK_DATA_DIR =="
# Skipped when it's already there (scripts/deploy-common.sh).
ensure_nltk_words "$NLTK_DATA_DIR"

if (( ! SKIP_DATABASE )); then
    # After the corpus: the population pass names species with it.
    offer_population_pass
fi

echo
echo "== 4/8: Making the web app's Python files and the shell scripts executable =="
# No -maxdepth: every *.py under src/html/, at any subdirectory depth
# (src/html/web/*.py included), needs this -- a previous version of this
# line was restricted to the top level only, which silently left
# src/html/web/*.py non-executable/unreadable-as-intended after every
# install. examples/apache/set-permissions.sh (below) re-does this same walk
# anyway with the correct final ownership, but doing it correctly here
# too means a plain `sudo ./install.sh` is never the reason this is wrong.
find "$HTML_DIR" -name '*.py' -exec chmod +x {} +
# Every *.sh anywhere in the repo (this script, update.sh,
# examples/apache/set-permissions.sh, and any future one), not a hardcoded list --
# git checkouts made from a `core.fileMode=false` machine silently drop
# the executable bit on ANY file type, not just src/html/'s .py scripts (see
# docs/TODO.md's "Deployment bugs found in production"), so a script added
# later doesn't need this list updated to be covered.
find "$SCRIPT_DIR" -name '*.sh' -exec chmod +x {} +

echo
echo "== 5/8: Enabling Apache's wsgi, headers and deflate modules (macOS: the gunicorn daemon) =="
# Installs mod_wsgi (libapache2-mod-wsgi-py3) if it's missing, and
# enables whichever of the three aren't already (scripts/deploy-common.sh).
APACHE_NEEDS_RESTART=0
if is_macos; then
    ensure_gunicorn_daemon
else
    ensure_apache_modules
fi

echo
echo "== 6/8: Setting directory ownership/permissions for Apache =="
# Also sets config.json (DB password, secret_key) to root:<apache group>, 640.
PYTHON="$PYTHON" "$SCRIPT_DIR/examples/apache/set-permissions.sh" "$HTML_DIR" "$DB_DIR"

echo
echo "== 7/8: Creating the Galaxy Map tile cache directory =="
PYTHON="$PYTHON" "$SCRIPT_DIR/examples/apache/create-cache-dir.sh"

echo
echo "== 8/8: Setting up the debug log and its rotation =="
# The logs never stop an install (OPS.5): setup-debug-log.sh warns, with
# the commands that fix it, for anything it can't set up itself.
PYTHON="$PYTHON" "$SCRIPT_DIR/examples/apache/setup-debug-log.sh" \
    || echo "warning: the log setup didn't finish (see above); the install carries on." >&2

echo
echo "Checking that the web app imports:"
check_app_imports

echo
echo "Checking the Redis server (config.json's redis.url):"
ensure_redis

if is_macos; then
    cat <<EOF

------------------------------------------------------------------------
Install steps complete. gunicorn serves the site on 127.0.0.1:8000
(check: curl -s http://127.0.0.1:8000/api/health). What remains is the
web server in front, never set up automatically:

  1. Copy the nginx example and edit server_name and the certificate paths:

       cp "$SCRIPT_DIR/examples/macos/planetgen-nginx.conf" "\$(brew --prefix)/etc/nginx/servers/planetgen.conf"

  2. Set "proxy_fix": {"x_for": 1, "x_proto": 1, "x_host": 0} in
     config.json, and "user _www _www;" at the top of nginx.conf, then:

       sudo nginx -t && sudo brew services start nginx

  3. Log in at https://<server_name>/login with the admin username and
     password printed once in step 2/8 above, and change both.

See docs/deployment/macos.md for the details.
------------------------------------------------------------------------
EOF
elif [[ ! -f /etc/apache2/sites-available/planetgen.conf ]]; then
    cat <<EOF

------------------------------------------------------------------------
Install steps complete. One manual step remains -- create the Apache2
site from the example config (never done automatically, since
ServerName/TLS/logging are your call):

  1. Copy the example and edit it (at minimum, set ServerName):

       sudo cp "$SCRIPT_DIR/examples/apache/planetgen.conf.example" /etc/apache2/sites-available/planetgen.conf
       sudo \${EDITOR:-nano} /etc/apache2/sites-available/planetgen.conf

  2. Enable the site and reload Apache:

       sudo a2ensite planetgen
       sudo systemctl reload apache2

  3. Log in at https://<ServerName>/login with the admin username and
     password printed once in step 2/8 above, and change both (the
     admin pages stay locked until you do). Lost it? See "Resetting the
     admin login" in docs/api.md.

See docs/deployment/apache.md and docs/html-interface.md for more detail.
------------------------------------------------------------------------
EOF
else
    echo
    echo "Done. (/etc/apache2/sites-available/planetgen.conf already exists --"
    echo "not touching it; reload Apache yourself if this update needs it:"
    echo "  sudo systemctl reload apache2)"

    # install.sh/update.sh deliberately never overwrite an existing site
    # file (ServerName/TLS/logging are the admin's own edits), so a site
    # file from before the CGI pages were removed keeps its CGI rules
    # until someone edits it. They do no harm, but they keep old
    # `/<name>.py` bookmarks from reaching the app's 301 to the page
    # that replaced them (Apache answers 404 itself instead).
    if grep -q 'ScriptAliasMatch' /etc/apache2/sites-available/planetgen.conf 2>/dev/null; then
        cat <<'EOF'

------------------------------------------------------------------------
NOTE: /etc/apache2/sites-available/planetgen.conf still has the old CGI
rules. Every page is served by the Flask app now; the CGI scripts are
gone. Remove from that file:

    ScriptAliasMatch "^/((?!wsgi\.py$)[a-z_]+\.py)$" ...
    Options +ExecCGI            (in <Directory .../src/html>; keep -Indexes)
    AddHandler cgi-script .py   (same block)
    <Files "wsgi.py"> SetHandler wsgi-script </Files>   (only once the
                                AddHandler line above is gone)

and any SetEnv PLANETGEN_* lines (they only ever reached the CGI
pages). Keep Alias /static/ and WSGIScriptAlias /. Compare with
examples/apache/planetgen.conf.example, then:

    sudo apache2ctl configtest && sudo systemctl reload apache2
------------------------------------------------------------------------
EOF
    fi
fi
