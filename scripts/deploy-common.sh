# scripts/deploy-common.sh
#
# Checks shared by install.sh and update.sh, sourced by both so the two
# can't drift apart. Each one looks first and only changes what is
# missing, so running it again on a working server does nothing.
# (Python libraries have their own script, install-python-deps.sh.)
# scripts/deploy-common.ps1 is the Windows counterpart for install.ps1
# and update.ps1: a change to a step here belongs there too.
#
# Two platforms: Linux with Debian's Apache and mod_wsgi
# (docs/deployment/apache.md), and macOS with gunicorn under launchd and
# Homebrew's nginx in front (docs/deployment/macos.md). Written for
# macOS's bash 3.2 as well: no mapfile or associative arrays, and empty
# arrays expanded as ${a[@]+"${a[@]}"}.
#
# Expects PYTHON and SCRIPT_DIR (the repo root) to be set, and to run as
# root.

MACOS_VENV="${PLANETGEN_VENV:-/usr/local/planetgen/venv}"
MACOS_LOG_DIR=/usr/local/planetgen/log
GUNICORN_LABEL=org.planetgen.gunicorn
GUNICORN_PLIST="/Library/LaunchDaemons/$GUNICORN_LABEL.plist"

is_macos() {
    [[ "$(uname -s)" == Darwin ]]
}

# The Python to set up the site with when PYTHON isn't given: on macOS
# Homebrew's python3 (Apple Silicon, then Intel) rather than the Command
# Line Tools' /usr/bin/python3; elsewhere python3 on the PATH.
default_python() {
    local candidate
    if is_macos; then
        for candidate in /opt/homebrew/bin/python3 /usr/local/bin/python3; do
            if [[ -x "$candidate" ]]; then
                echo "$candidate"
                return
            fi
        done
    fi
    command -v python3 || command -v python || true
}

# After install-python-deps.sh: on macOS the site runs from the venv it
# made, so every later step uses that Python.
use_site_python() {
    if is_macos; then
        PYTHON="$MACOS_VENV/bin/python"
    fi
}

# Copies an example launchd plist into /Library/LaunchDaemons with the
# checkout and venv paths put in, when it isn't there already (an
# existing one is the admin's; it's never overwritten), and loads it.
# Args: example file, installed name, then extra sed expressions.
install_launchd_plist() {
    local example="$1" target="/Library/LaunchDaemons/$2" label="${2%.plist}"
    shift 2
    if [[ -f "$target" ]]; then
        echo "launchd: $target present"
        return 0
    fi
    sed -e "s|/var/lib/planetGen|$SCRIPT_DIR|g" -e "s|/usr/local/planetgen/venv|$MACOS_VENV|g" "$@" \
        "$example" > "$target"
    chown root:wheel "$target"
    chmod 644 "$target"
    launchctl bootstrap system "$target" 2>/dev/null || launchctl load -w "$target" || true
    echo "launchd: installed and loaded $target ($label)"
}

# macOS's counterpart of the Apache module check: gunicorn under launchd
# (examples/macos/org.planetgen.gunicorn.plist), serving on 127.0.0.1:8000
# for nginx. Installed and started the first time only. Its log
# directory belongs to _www.
ensure_gunicorn_daemon() {
    # gunicorn's own output; a folder that can't be made only warns (OPS.5).
    if ! { mkdir -p "$MACOS_LOG_DIR" && chown _www:_www "$MACOS_LOG_DIR"; }; then
        echo "warning: couldn't set up gunicorn's log folder $MACOS_LOG_DIR; launchd can't start" >&2
        echo "  gunicorn until it exists. Fix it with:" >&2
        echo "    sudo mkdir -p $MACOS_LOG_DIR" >&2
        echo "    sudo chown _www:_www $MACOS_LOG_DIR" >&2
    fi
    install_launchd_plist "$SCRIPT_DIR/examples/macos/$GUNICORN_LABEL.plist" "$GUNICORN_LABEL.plist"
}

# The NLTK 'words' corpus, in a shared, world-readable directory (not a
# per-user home directory) so it works for every user that imports
# the generator: a login shell running the generator, or Apache's
# www-data running the web app. planetgen/names/wordlists.py checks
# nltk.data.find() before ever calling download(), so once this
# directory (on nltk's default search path) has it, nothing downloads
# again. See docs/TODO.md's "Deployment bugs found in production" for
# the PermissionError on /var/www/nltk_data this avoids.
ensure_nltk_words() {
    local dir="$1"
    mkdir -p "$dir"
    if "$PYTHON" -c "import nltk, sys; nltk.data.find('corpora/words', paths=[sys.argv[1]])" "$dir" >/dev/null 2>&1; then
        echo "NLTK 'words' corpus: present in $dir"
    else
        # A plain nltk.download() rather than `python -m nltk.downloader`,
        # which prints a spurious "found in sys.modules" RuntimeWarning
        # (nltk's __init__ already imported nltk.downloader).
        "$PYTHON" -c "import nltk, sys; sys.exit(0 if nltk.download('words', download_dir=sys.argv[1]) else 1)" "$dir"
        echo "NLTK 'words' corpus: installed in $dir"
    fi
    chmod -R a+rX "$dir"
}

# The Redis server at config.json's redis.url (OPS.21), which the work
# queue and the rate limits will use. On Linux with a local address it
# installs redis-server with apt and starts it when nothing answers; on
# macOS (Homebrew refuses to run as root) and for a remote address it
# only says what to run. Nothing uses Redis yet, so a server that doesn't
# answer warns rather than stopping the install or update.
ensure_redis() {
    local url host
    url="$("$PYTHON" -c "from planetgen.util.settings import get_settings; print(get_settings().redis.url)")" || {
        echo "warning: couldn't read redis.url from config.json; skipping the Redis check." >&2
        return 0
    }
    host="$("$PYTHON" -c "import sys; from urllib.parse import urlsplit; print(urlsplit(sys.argv[1]).hostname or '')" "$url")"
    if redis_answers "$url"; then
        echo "Redis: answering at $url"
        return 0
    fi
    if [[ "$host" == 127.0.0.1 || "$host" == localhost || "$host" == ::1 ]] && ! is_macos \
            && command -v apt-get >/dev/null 2>&1; then
        if ! command -v redis-server >/dev/null 2>&1; then
            echo "Redis: installing redis-server."
            DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends redis-server || true
        fi
        if command -v systemctl >/dev/null 2>&1; then
            systemctl enable --now redis-server >/dev/null 2>&1 || true
        fi
        if redis_answers "$url"; then
            echo "Redis: answering at $url"
            return 0
        fi
    fi
    echo "warning: no Redis server answers at $url (config.json's redis.url)." >&2
    if is_macos; then
        echo "  Install and start it as your own user: brew install redis && brew services start redis" >&2
    else
        echo "  Install and start it: sudo apt install redis-server && sudo systemctl enable --now redis-server" >&2
    fi
    echo "  Nothing needs it yet; the work queue will (docs/deployment/README.md)." >&2
}

# Whether a Redis server answers PING at the URL in $1, using the redis
# library install-python-deps.sh installed.
redis_answers() {
    "$PYTHON" -c "import sys, redis; redis.Redis.from_url(sys.argv[1], socket_connect_timeout=3).ping()" "$1" >/dev/null 2>&1
}

# Apache's modules: headers (static/'s Cache-Control/nosniff lines in
# examples/apache/planetgen.conf.example), deflate (its compression
# block) and wsgi (runs the Flask app: every page and the API). mod_wsgi
# comes from libapache2-mod-wsgi-py3, installed here when apt can.
# Warns rather than fails when Apache itself isn't installed yet. Sets
# APACHE_NEEDS_RESTART=1 when it enabled anything, so the caller can say
# Apache needs a restart.
ensure_apache_modules() {
    if ! command -v a2enmod >/dev/null 2>&1; then
        echo "warning: a2enmod not found -- is apache2 installed?" >&2
        echo "  Try: sudo apt install apache2" >&2
        return 0
    fi
    if [[ ! -e /etc/apache2/mods-available/wsgi.load ]] && command -v apt-get >/dev/null 2>&1; then
        echo "mod_wsgi is not installed: installing libapache2-mod-wsgi-py3."
        # The package enables the module itself, so a2query below finds
        # it already on; Apache still has to be restarted to load it.
        DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends libapache2-mod-wsgi-py3 \
            && APACHE_NEEDS_RESTART=1 || true
    fi
    local mod enabled=()
    for mod in headers deflate wsgi; do
        if a2query -q -m "$mod" 2>/dev/null; then
            echo "Apache module $mod: enabled"
        elif a2enmod -q "$mod" >/dev/null 2>&1; then
            echo "Apache module $mod: enabled now"
            enabled+=("$mod")
        else
            echo "warning: could not enable Apache module $mod." >&2
            if [[ "$mod" == wsgi ]]; then
                echo "  Install it first: sudo apt install libapache2-mod-wsgi-py3 && sudo a2enmod wsgi" >&2
            fi
        fi
    done
    if (( ${#enabled[@]} )); then
        APACHE_NEEDS_RESTART=1
    fi
    check_mod_wsgi_python
}

# Reloads Apache itself when it is running (OPS.8), or restarts it when a
# module was just enabled (APACHE_NEEDS_RESTART=1). Prints "reloaded" or
# "restarted" on success; returns 1 when it didn't act (no systemctl, Apache
# not running, or the command failed), so the caller prints the command for
# the admin to run. The caller must be root; update.sh always is.
reload_apache_if_running() {
    local verb=reload
    (( APACHE_NEEDS_RESTART )) && verb=restart
    command -v systemctl >/dev/null 2>&1 || return 1
    systemctl is-active --quiet apache2 2>/dev/null || return 1
    if systemctl "$verb" apache2 >/dev/null 2>&1; then
        echo "${verb}ed Apache."
        return 0
    fi
    echo "warning: systemctl $verb apache2 failed." >&2
    return 1
}

# mod_wsgi embeds the Python it was built against, not whatever `python3`
# is. The libraries (which live in that Python's own
# site-packages) are set up for $PYTHON, so the two must be the same
# version or Apache won't see them. Warns rather than fails: the fix is
# a package choice for the admin.
check_mod_wsgi_python() {
    local so=/usr/lib/apache2/modules/mod_wsgi.so built ours
    [[ -e "$so" ]] && command -v ldd >/dev/null 2>&1 || return 0
    built="$(ldd "$so" 2>/dev/null | grep -o 'libpython[0-9]*\.[0-9]*' | head -n1 | sed 's/libpython//')"
    ours="$("$PYTHON" -c 'import sys; print("%d.%d" % sys.version_info[:2])')"
    if [[ -z "$built" ]]; then
        return 0
    elif [[ "$built" == "$ours" ]]; then
        echo "mod_wsgi runs Python $built, the same as $PYTHON."
    else
        echo "warning: mod_wsgi is built for Python $built, but the libraries were" >&2
        echo "  checked for $PYTHON (Python $ours). Apache won't see them. Either install" >&2
        echo "  the libapache2-mod-wsgi-py3 that matches $PYTHON, or rerun with" >&2
        echo "  PYTHON=/usr/bin/python$built." >&2
    fi
}

# Builds the Galaxy Map's opening view into the tile cache (MAP.134), in the
# background as Apache's user (the cache belongs to it), so the update
# doesn't wait for it and the first visitor after it doesn't either. At
# most 30 minutes; the output goes to $WARM_MAP_LOG. Never fails the
# update: the first visit builds the view itself when this doesn't.
warm_galaxy_map() {
    local user="" runner=()
    WARM_MAP_LOG="${WARM_MAP_LOG:-/tmp/planetgen-warm-map.log}"
    if [[ -r "$SCRIPT_DIR/examples/apache/apache-identity.sh" ]]; then
        # shellcheck source=../examples/apache/apache-identity.sh
        source "$SCRIPT_DIR/examples/apache/apache-identity.sh"
        user="$(detect_apache_group 2>/dev/null | awk '{print $1}')"
    fi
    if [[ -n "$user" ]] && id "$user" >/dev/null 2>&1 && command -v runuser >/dev/null 2>&1; then
        runner=(runuser -u "$user" --)
    elif [[ -n "$user" ]] && id "$user" >/dev/null 2>&1 && is_macos; then
        runner=(sudo -u "$user" --)
    fi
    local limit=()
    if command -v timeout >/dev/null 2>&1; then
        limit=(timeout 1800)
    fi
    (cd / && nohup ${runner[@]+"${runner[@]}"} ${limit[@]+"${limit[@]}"} "$PYTHON" -m planetgen.cli.warm_map \
        >"$WARM_MAP_LOG" 2>&1 </dev/null &) || true
    echo "Building the Galaxy Map's opening view in the background (log: $WARM_MAP_LOG)."
}

# Imports the web app and the generator package with the interpreter the
# site runs under, as Apache's own user when that user exists, so a
# library that's installed but unreadable to www-data (or the checkout's
# own code failing to import) shows up here rather than as a 500.
check_app_imports() {
    local user="" runner=()
    if [[ -r "$SCRIPT_DIR/examples/apache/apache-identity.sh" ]]; then
        # shellcheck source=../examples/apache/apache-identity.sh
        source "$SCRIPT_DIR/examples/apache/apache-identity.sh"
        user="$(detect_apache_group 2>/dev/null | awk '{print $1}')"
    fi
    if [[ -n "$user" ]] && id "$user" >/dev/null 2>&1 && command -v runuser >/dev/null 2>&1; then
        runner=(runuser -u "$user" --)
    elif [[ -n "$user" ]] && id "$user" >/dev/null 2>&1 && is_macos; then
        runner=(sudo -u "$user" --)  # macOS has no runuser
    else
        user=root
    fi
    if (cd / && ${runner[@]+"${runner[@]}"} "$PYTHON" - <<'EOF'
import planetgen.generation.system  # noqa: F401
from planetgen.web.app import create_app  # noqa: F401
EOF
    ); then
        echo "The web app and the generator import cleanly as $user with $PYTHON."
    else
        echo "error: the web app does not import as $user with $PYTHON (see above)." >&2
        return 1
    fi
}

# Runs one of planetGen's command-line tools, `python3 -m planetgen.cli.NAME
# ARGS...` (the editable install makes planetgen import from anywhere).
run_cli() {
    "$PYTHON" -m "planetgen.cli.$1" "${@:2}"
}

# Brings the configured database up to the current schema. When a
# migration is pending, first asks whether to delete the galaxy data
# instead: y wipes every generated sector and system (planetgen.cli.reset, the
# same as the Generate page's Reset; admin logins are kept) and the empty
# database is then brought to the current schema. Anything else, no
# answer within 30 seconds, or no terminal to ask on (the maintenance
# timer) keeps the data and migrates it. Nothing is asked when the
# database is already current.
migrate_or_reset_db() {
    local status current target pending database answer=""
    if ! status="$(run_cli migrate --status)"; then
        return 1
    fi
    read -r current target pending database <<< "$status"
    if (( pending > 0 )); then
        echo "Database '$database' is at schema v$current; this version needs v$target ($pending migration step(s))."
        if [[ -t 0 ]]; then
            read -r -t 30 -p "Delete all galaxy data in '$database' instead of migrating it? [y/N] (default N in 30s): " answer \
                || { echo; answer=""; }
        else
            echo "(No terminal to ask on: keeping the data and migrating it.)"
        fi
        if [[ "$answer" =~ ^[Yy]([Ee][Ss])?$ ]]; then
            echo "Deleting the galaxy data in '$database' (admin logins are kept)."
            run_cli reset --yes
        else
            echo "Keeping the data and migrating it."
        fi
    fi
    run_cli migrate
}

# Records the version key this update runs under, one row per galaxy in the
# control database's history (OPS.13; the last 10 are kept), after the
# schemas are current. A failure only warns: the update carries on.
record_version_key() {
    run_cli version_history || echo "warning: could not record the version key (see above); the update carries on." >&2
}

# The math check (`planetgen check-math`, TEST.68): known answers from
# real astronomy, identities and sampler distributions. A failure only
# warns -- the update carries on and the site keeps serving -- but sets
# MATH_CHECK_FAILED=1, so offer_population_pass skips the population pass
# (which would refuse anyway) and the closing message repeats the warning.
MATH_CHECK_FAILED=0
check_math() {
    if run_cli generate check-math; then
        MATH_CHECK_FAILED=0
    else
        MATH_CHECK_FAILED=1
        echo "warning: the math check failed (above). Bulk generation refuses to start until it passes;" >&2
        echo "         run 'planetgen check-math -v' for every check." >&2
    fi
}

# Optionally runs the population pass (`planetgen population`: species,
# civilizations and territories, docs/design/population-and-politics.md)
# over the stored galaxy. Off by default (Boss, 2026-10-01): it runs only
# when POPULATION=1 (or yes) is set, or when someone answers y to the
# prompt (y/N, 30 seconds, default N). With no terminal it is skipped.
# `planetgen population` can always be run by hand later.
offer_population_pass() {
    local answer=""
    if [[ "${MATH_CHECK_FAILED:-0}" == 1 ]]; then
        echo "Skipping the population pass: the math check failed."
        return 0
    fi
    case "${POPULATION:-}" in
        1|[Yy]|[Yy][Ee][Ss]) answer="y" ;;
        *)
            if [[ -t 0 ]]; then
                read -r -t 30 -p "Run the population pass now (species, civilizations, territories)? [y/N] (default N in 30s): " answer \
                    || { echo; answer=""; }
            fi ;;
    esac
    if [[ "$answer" =~ ^[Yy]([Ee][Ss])?$ ]]; then
        echo "Running the population pass."
        run_cli generate population
    else
        echo "Skipping the population pass (run 'planetgen population' any time, or set POPULATION=1)."
    fi
}
