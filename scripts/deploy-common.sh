# scripts/deploy-common.sh
#
# Checks shared by install.sh and update.sh, sourced by both so the two
# can't drift apart. Each one looks first and only changes what is
# missing, so running it again on a working server does nothing.
# (Python libraries have their own script, install-python-deps.sh.)
#
# Two platforms: Linux with Debian's Apache and mod_wsgi
# (docs/deployment/apache.md), and macOS with gunicorn under launchd and
# Homebrew's nginx in front (docs/deployment/macos.md). Written for
# macOS's bash 3.2 as well: no mapfile or associative arrays, and empty
# arrays expanded as ${a[@]+"${a[@]}"}.
#
# Expects PYTHON and SCRIPT_DIR (the repo root) to be set, and to run as
# root.

VENV="${PLANETGEN_VENV:-/usr/local/planetgen/venv}"
MACOS_LOG_DIR=/usr/local/planetgen/log
GUNICORN_LABEL=org.planetgen.gunicorn
GUNICORN_PLIST="/Library/LaunchDaemons/$GUNICORN_LABEL.plist"

is_macos() {
    [[ "$(uname -s)" == Darwin ]]
}

# Says so when a Linux host has no swap on (PERF.181): with none, a burst of memory ends in the kernel killing MySQL,
# Redis or the web server. Only a warning; scripts/enable-swap.sh fixes it.
warn_if_no_swap() {
    is_macos && return 0
    command -v swapon >/dev/null 2>&1 || return 0
    if [[ -z "$(swapon --show=NAME --noheadings 2>/dev/null)" ]]; then
        echo "warning: this host has no swap. A memory spike then makes the kernel kill MySQL, Redis or the web" >&2
        echo "  server instead of slowing down. Fix: sudo $SCRIPT_DIR/scripts/enable-swap.sh  (see INSTALL.md, Swap)" >&2
    fi
    return 0
}

# The Python the venv is made from when PYTHON isn't given: on macOS
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

# After install-python-deps.sh: the site runs from the venv it made
# (OPS.45), so every later step uses that Python.
use_site_python() {
    PYTHON="$VENV/bin/python"
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
    sed -e "s|/var/lib/planetGen|$SCRIPT_DIR|g" -e "s|/usr/local/planetgen/venv|$VENV|g" "$@" \
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

# The static data the generator reads (GEN.548, planetgen.datasets): downloaded
# once, here, into a cache folder outside the checkout and made readable for
# Apache's user; nothing downloads while the site runs. A dataset whose host
# does not answer stays missing and DATASETS_MISSING is set so the end of the
# script says so (bulk generation refuses to start without it, GEN.590); it
# never stops an install or an update.
DATASETS_MISSING=0
ensure_datasets() {
    local dir="${PLANETGEN_DATA_DIR:-/var/cache/planetgen/datasets}"
    mkdir -p "$dir"
    if ! PLANETGEN_DATA_DIR="$dir" run_cli datasets fetch; then
        DATASETS_MISSING=1
        echo "warning: some data could not be downloaded (above); bulk generation refuses to start until the next" >&2
        echo "         update fetches it ('python3 -m planetgen.cli.datasets status' says which)." >&2
    fi
    chmod -R a+rX "$dir"
}

# The Redis server at config.json's redis.url (OPS.21), which the work
# queue, the Generate page's jobs and the rate limits use. On Linux with a
# local address it installs redis-server with apt and starts it when nothing
# answers; on macOS (Homebrew refuses to run as root) and for a remote address
# it only says what to run. A server that doesn't answer warns rather than
# stopping the install or update.
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
    echo "  Without it the Generate page starts no job and a run cannot use more than one worker (docs/deployment/README.md)." >&2
}

# Warns, with the exact line to add, when the running Redis server's
# settings don't suit planetGen (OPS.48): append-only file off, a memory
# policy that evicts, no memory limit, or listening beyond loopback. Never
# edits the Redis config (examples/redis/planetgen.conf is the example).
check_redis_settings() {
    local url
    url="$("$PYTHON" -c "from planetgen.util.settings import get_settings; print(get_settings().redis.url)")" || return 0
    redis_answers "$url" || return 0
    "$PYTHON" "$SCRIPT_DIR/scripts/check_redis_settings.py" "$url" || true
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
# APACHE_NEEDS_RESTART=1 when it installed or enabled mod_wsgi, so the caller can say
# a restart (only mod_wsgi does; the other modules load on a reload).
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
    # Only mod_wsgi needs a restart to load; a graceful reload takes the other modules (OPS.29).
    local was
    for was in ${enabled[@]+"${enabled[@]}"}; do
        [[ "$was" == wsgi ]] && APACHE_NEEDS_RESTART=1
    done
    check_apache_venv_rules
    check_mod_wsgi_python
}

# Reloads what serves the site, when it is running (OPS.8, OPS.29), and says what it did. It never starts anything:
# Apache (unit apache2 or httpd) is checked with its own configtest first and reloaded (a graceful reload loads new
# modules too), or restarted only when this update just installed or enabled mod_wsgi (APACHE_NEEDS_RESTART=1, which
# ensure_apache_modules sets for that alone); the update guard is held until the script ends, so no Generate job
# runs then. gunicorn (the planetgen-gunicorn unit on Linux, launchd's org.planetgen.gunicorn on macOS) gets a SIGHUP,
# which starts new workers and leaves a Generate job alone. Afterwards each answers /api/health (PLANETGEN_HEALTH_URL
# names another address; HEALTH_WAIT_SECONDS, default 20, is how long it waits).
# A reload aborts the requests a mod_wsgi daemon is serving after about 4 seconds (touching src/html/wsgi.py does
# not), which is why a reload is quick and a long request at that moment can fail.
# Returns 0 when everything reloaded answers, 1 when a step failed (the caller prints the command and ends non-zero),
# 2 when nothing that serves the site is running (the caller prints the command, the update is fine), 3 when every
# reload succeeded but the /api/health probe never got a good answer (the probe can be wrong for reasons that are not
# the new code's: another virtual host owns the address, no listener on the https port; the caller warns and opens the
# site, since the database migration has finished by then).
reload_web_if_running() {
    local acted=0 failed=0 unhealthy=0 unit="" verb=reload check="" health_url="" health_from=""
    # The one address every health check asks: config.json's api_base_url (or PLANETGEN_HEALTH_URL). No loopback default.
    resolve_api_health || return 1
    health_url="$API_HEALTH_URL"
    health_from="$API_HEALTH_FROM"
    if is_macos; then
        if launchctl print "system/$GUNICORN_LABEL" >/dev/null 2>&1; then
            if launchctl kill SIGHUP "system/$GUNICORN_LABEL" >/dev/null 2>&1; then
                echo "reloaded gunicorn."
                acted=1
                wait_for_health "$health_url" gunicorn "$health_from" || unhealthy=1
            else
                echo "warning: launchctl kill SIGHUP system/$GUNICORN_LABEL failed." >&2
                failed=1
            fi
        fi
    elif command -v systemctl >/dev/null 2>&1; then
        for unit in apache2 httpd; do
            systemctl is-active --quiet "$unit" 2>/dev/null && break
            unit=""
        done
        if [[ -n "$unit" ]]; then
            (( APACHE_NEEDS_RESTART )) && verb=restart
            acted=1
            echo "Testing Apache's configuration, then ${verb}ing it (a request still running is cut off):"
            if ! apache_configtest; then
                echo "warning: Apache's configuration test failed, so Apache was not ${verb}ed (the running one keeps serving)." >&2
                failed=1
            elif systemctl "$verb" "$unit" >/dev/null 2>&1; then
                echo "${verb}ed Apache."
                wait_for_health "$health_url" Apache "$health_from" || unhealthy=1
            else
                echo "warning: systemctl $verb $unit failed." >&2
                failed=1
            fi
        fi
        if systemctl is-active --quiet planetgen-gunicorn 2>/dev/null; then
            acted=1
            if systemctl reload planetgen-gunicorn >/dev/null 2>&1; then
                echo "reloaded gunicorn."
                wait_for_health "$health_url" gunicorn "$health_from" || unhealthy=1
            else
                echo "warning: systemctl reload planetgen-gunicorn failed." >&2
                failed=1
            fi
        fi
    fi
    (( failed )) && return 1
    (( unhealthy )) && return 3
    (( acted )) && return 0
    return 2
}

# The commands that reload what serves the site by hand, for the host this is (OPS.29).
print_reload_commands() {
    if is_macos; then
        echo "  sudo launchctl kill SIGHUP system/$GUNICORN_LABEL"
    else
        echo "  sudo apache2ctl configtest && sudo systemctl reload apache2     (Apache)"
        echo "  sudo systemctl reload planetgen-gunicorn                        (nginx or Caddy hosts)"
    fi
}

# Apache's own configuration test (apache2ctl, apachectl or httpd -t); true when there is none to run.
apache_configtest() {
    local tool
    for tool in apache2ctl apachectl; do
        if command -v "$tool" >/dev/null 2>&1; then
            "$tool" configtest >/dev/null 2>&1
            return
        fi
    done
    if command -v httpd >/dev/null 2>&1; then
        httpd -t >/dev/null 2>&1
        return
    fi
    return 0
}

# The one address the update asks whether the site answers (OPS.86): config.json's `api_base_url` (the API's /api mount,
# `PLANETGEN_API_BASE_URL` in the environment) with /health added, or PLANETGEN_HEALTH_URL when that is given. There is no
# loopback default and no second address to try: a host whose API is at 127.0.0.1 says so in config.json. Sets
# API_HEALTH_URL and API_HEALTH_FROM (the setting the address came from, named in every message about it) and returns 0;
# with no address configured it says so in plain words on stderr and returns 1, changing nothing.
resolve_api_health() {
    API_HEALTH_URL="" API_HEALTH_FROM=""
    if [[ -n "${PLANETGEN_HEALTH_URL:-}" ]]; then
        API_HEALTH_URL="$PLANETGEN_HEALTH_URL" API_HEALTH_FROM="PLANETGEN_HEALTH_URL"
        return 0
    fi
    local base=""
    base="$("${PYTHON:-python3}" -c 'from planetgen.util.settings import get_settings; print(get_settings().api_base_url)' 2>/dev/null)" || base=""
    base="${base%/}"
    if [[ -z "$base" ]]; then
        echo "error: no API address is configured: set api_base_url in config.json (for example \"https://your-site/api\", the address Apache serves the API at;" >&2
        echo "  PLANETGEN_API_BASE_URL in the environment does the same), or PLANETGEN_HEALTH_URL for the health check alone." >&2
        echo "  The update asks only that address whether the site answers and never tries 127.0.0.1 or localhost unless it is the one configured." >&2
        return 1
    fi
    API_HEALTH_URL="$base/health" API_HEALTH_FROM="api_base_url in config.json"
    return 0
}

# What curl's exit status says when it got no HTTP answer at all (the "HTTP 000" in a warning), in plain words.
curl_failure_words() {
    case "${1:-0}" in
        0) echo "no answer" ;;
        6) echo "the name could not be looked up" ;;
        7) echo "the connection was refused: nothing listens on that address and port" ;;
        28) echo "no answer in time" ;;
        35|51|53|54|58|59|60|77|80|82|83|90|91) echo "the TLS (https) handshake or certificate failed, curl exit $1" ;;
        47) echo "too many redirects" ;;
        52) echo "the server closed the connection without an answer" ;;
        56) echo "the connection was reset" ;;
        *) echo "curl exit $1" ;;
    esac
}

# Waits up to HEALTH_WAIT_SECONDS for $1 (an address) to answer below 400, naming $2 in the warning and $3 as the setting
# the address came from (OPS.86). It asks that one address and no other. A passing answer while the site is closed for the
# update says so (and why), on success and in the warning alike: /api/health stays open while the site is closed
# (OPS.75), so a reload is checked while the update's own maintenance mode is still on.
# Returns 0 when it answers, 1 when it never did (or without curl, 0).
wait_for_health() {
    command -v curl >/dev/null 2>&1 || return 0
    # The reload has just restarted the workers: they get HEALTH_SETTLE_SECONDS (default 10) to start before the first probe.
    echo "Giving the workers ${HEALTH_SETTLE_SECONDS:-10} seconds to start, then asking $1 (${3:-its setting}) for /api/health (up to ${HEALTH_WAIT_SECONDS:-20} seconds):"
    sleep "${HEALTH_SETTLE_SECONDS:-10}"
    local end=$(( SECONDS + ${HEALTH_WAIT_SECONDS:-20} )) reply code body rc
    while :; do
        rc=0
        reply="$(curl -s --max-time 5 -w '\n%{http_code}' "$1" 2>/dev/null)" || rc=$?
        code="${reply##*$'\n'}"
        body="${reply%$'\n'*}"
        if [[ "$code" =~ ^[23] ]]; then
            health_maintenance_note "$body"
            return 0
        fi
        if (( SECONDS >= end )); then
            local bare_note=""
            [[ "${code:-000}" == 000 ]] && bare_note=" (no HTTP answer at all: $(curl_failure_words "$rc"))"
            echo "warning: the health check of $2 failed: $1 (${3:-its setting}) answered HTTP ${code:-000}${bare_note} for ${HEALTH_WAIT_SECONDS:-20} seconds after the reload. Only that address was asked; check ${3:-the setting} if it is not where the site is served.$(health_maintenance_note "$body" brief)" >&2
            return 1
        fi
        (( SECONDS % 5 == 0 )) && echo "  no good answer yet (HTTP ${code:-000}); asking again, $(( end - SECONDS )) s left"
        sleep 1
    done
}

# Says whether a /api/health reply reports maintenance mode, and why; the reply is $1. With a second argument it only
# returns the sentence (for a warning), otherwise it prints it.
health_maintenance_note() {
    local body="$1" cause reason text
    [[ "$body" =~ \"maintenance\"[[:space:]]*:[[:space:]]*(true|\{[^}]*\"enabled\"[[:space:]]*:[[:space:]]*true) ]] || return 0
    cause="$(sed -n 's/.*"cause"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' <<<"$body" | head -n 1)"
    reason="$(sed -n 's/.*"reason"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' <<<"$body" | head -n 1)"
    text=" The site is in maintenance mode (cause: ${cause:-unknown}${reason:+, $reason})"
    if [[ -n "${2:-}" ]]; then
        echo -n "$text."
    else
        echo "${text# }; the update opens it again when it ends."
    fi
}

# The Apache files that can hold planetGen's vhost: the enabled sites and
# configuration of Debian/Ubuntu's apache2 and RHEL's httpd. The
# APACHE_CONF_DIRS and APACHE_MAIN_CONFS variables replace the defaults
# (the tests point them at a copy).
apache_config_files() {
    local dir file
    for dir in ${APACHE_CONF_DIRS:-/etc/apache2/sites-enabled /etc/apache2/conf-enabled /etc/httpd/conf.d}; do
        for file in "$dir"/*.conf; do
            [[ -r "$file" ]] && echo "$file"
        done
    done
    for file in ${APACHE_MAIN_CONFS:-/etc/apache2/apache2.conf /etc/httpd/conf/httpd.conf}; do
        [[ -r "$file" ]] && echo "$file"
    done
    return 0
}

# Runs a command as Apache's user (the one `detect_apache_group` finds), to
# see what the web server itself can reach. Returns 125 when that can't be
# done (no such user, no runuser or sudo).
as_apache_user() {
    local user="$1"
    shift
    id -u "$user" >/dev/null 2>&1 || return 125
    if command -v runuser >/dev/null 2>&1; then
        runuser -u "$user" -- "$@"
    elif command -v sudo >/dev/null 2>&1; then
        sudo -n -u "$user" "$@"
    else
        return 125
    fi
}

# Warns when Apache's user can't read the venv: every directory down to it
# has to be enterable (x), and the venv's own folders and pyvenv.cfg
# readable, by that user, or mod_wsgi's daemon starts with no Python
# libraries. Names the first path that blocks it and the command that opens
# it.
check_apache_user_reads_venv() {
    local user="" group="" line rc=0 blocked="" path="" args=()
    if [[ -r "$SCRIPT_DIR/examples/apache/apache-identity.sh" ]]; then
        # shellcheck source=../examples/apache/apache-identity.sh
        source "$SCRIPT_DIR/examples/apache/apache-identity.sh"
        line="$(detect_apache_group 2>/dev/null || true)"
        read -r user group <<< "$line"
    fi
    [[ -n "$user" ]] || return 0
    # "/", then each directory of the path, to the venv itself.
    args=(x /)
    local rest="${VENV#/}"
    while [[ -n "$rest" ]]; do
        path="$path/${rest%%/*}"
        args+=(x "$path")
        if [[ "$rest" == */* ]]; then rest="${rest#*/}"; else rest=""; fi
    done
    args+=(r "$VENV" r "$VENV/bin" x "$VENV/bin" r "$VENV/lib" x "$VENV/lib" r "$VENV/pyvenv.cfg")
    blocked="$(as_apache_user "$user" sh -c '
        while [ "$#" -gt 1 ]; do
            if [ "$1" = x ]; then test -x "$2" || { echo "$2"; exit 1; }
            else test -r "$2" || { echo "$2"; exit 1; }; fi
            shift 2
        done' sh "${args[@]}")" || rc=$?
    if (( rc == 125 )); then
        return 0
    elif (( rc != 0 )); then
        echo "warning: Apache's user $user can't read ${blocked:-the venv}, so mod_wsgi can't use $VENV." >&2
        echo "  Every directory down to the venv needs to be enterable by it (for example:" >&2
        echo "  sudo chmod o+rx ${blocked:-$VENV}), and the venv readable. Then restart Apache." >&2
    else
        echo "Apache's user $user can read $VENV."
    fi
}

# Checks planetGen's Apache vhost against mod_wsgi's rules for running an app
# in a virtual environment (OPS.45, OPS.46): the app runs in a daemon process
# group (WSGIDaemonProcess with `process-group=` on the WSGIScriptAlias, or
# a WSGIProcessGroup), that group's `python-home` is the venv's root (not
# bin/python or a parent), the app runs in the main interpreter
# (`application-group=%{GLOBAL}` or WSGIApplicationGroup %{GLOBAL}), and
# Apache's user can reach the venv. The vhost is the admin's own file: each
# problem is only a warning, with the line to add.
check_apache_venv_rules() {
    local conf="" file active alias_line group home daemon
    for file in $(apache_config_files); do
        if grep -vE '^[[:space:]]*#' "$file" | grep -qE '^[[:space:]]*WSGIScriptAlias[[:space:]].*src/html/wsgi\.py'; then
            conf="$file"
            break
        fi
    done
    [[ -n "$conf" ]] || return 0
    active="$(grep -vE '^[[:space:]]*#' "$conf")"
    alias_line="$(grep -E '^[[:space:]]*WSGIScriptAlias[[:space:]].*src/html/wsgi\.py' <<< "$active" | head -n1)"
    group="$(sed -n 's/.*[[:space:]]process-group=\([^[:space:]]*\).*/\1/p' <<< "$alias_line")"
    [[ -n "$group" ]] || group="$(sed -n 's/^[[:space:]]*WSGIProcessGroup[[:space:]]\{1,\}\([^[:space:]]*\).*/\1/p' <<< "$active" | head -n1)"
    if [[ -z "$group" || "$group" == %* ]]; then
        echo "warning: $conf runs planetGen in Apache's embedded mode, which can't use a virtual environment." >&2
        echo "  Declare a daemon process and use it (examples/apache/planetgen.conf.example):" >&2
        echo "    WSGIDaemonProcess planetgen-api processes=1 threads=5 request-timeout=60 python-home=$VENV" >&2
        echo "    WSGIScriptAlias / $SCRIPT_DIR/src/html/wsgi.py process-group=planetgen-api application-group=%{GLOBAL}" >&2
    else
        daemon="$(grep -E "^[[:space:]]*WSGIDaemonProcess[[:space:]]+$group([[:space:]]|\$)" <<< "$active" | head -n1 || true)"
        if [[ -z "$daemon" ]]; then
            echo "warning: $conf names the process group '$group' but has no WSGIDaemonProcess $group line." >&2
            echo "  Add: WSGIDaemonProcess $group processes=1 threads=5 request-timeout=60 python-home=$VENV" >&2
        else
            home="$(sed -n 's/.*[[:space:]]python-home=\([^[:space:]]*\).*/\1/p' <<< "$daemon")"
            [[ "$home" == / ]] || home="${home%/}"
            if [[ "$home" == "$VENV" ]]; then
                echo "Apache runs the site under $VENV (process group $group)."
            elif [[ -z "$home" ]]; then
                echo "warning: $conf does not run the site under planetGen's virtual environment." >&2
                echo "  Add python-home=$VENV to its WSGIDaemonProcess $group line (examples/apache/planetgen.conf.example)," >&2
                echo "  then restart Apache." >&2
            elif [[ "$home" == "$VENV"/bin || "$home" == "$VENV"/bin/* ]]; then
                echo "warning: $conf has python-home=$home, a path into the venv's bin/." >&2
                echo "  python-home is the venv's root, the folder that holds bin/ and lib/: python-home=$VENV" >&2
            else
                echo "warning: $conf has python-home=$home, not planetGen's virtual environment $VENV." >&2
                echo "  Set it to python-home=$VENV (the folder that holds bin/ and lib/), then restart Apache." >&2
            fi
        fi
        if ! grep -qF 'application-group=%{GLOBAL}' <<< "$alias_line" \
                && ! grep -qE '^[[:space:]]*WSGIApplicationGroup[[:space:]]+%\{GLOBAL\}' <<< "$active"; then
            echo "warning: $conf doesn't run planetGen in the main Python interpreter." >&2
            echo "  Add application-group=%{GLOBAL} to its WSGIScriptAlias line (or WSGIApplicationGroup %{GLOBAL})," >&2
            echo "  as numpy, scipy and numba's C code need it." >&2
        fi
    fi
    check_apache_user_reads_venv
}

# mod_wsgi embeds the Python it was built against, not whatever `python3`
# is. The libraries (which live in that Python's own
# site-packages) are set up for $PYTHON, so the two must be the same
# version or Apache won't see them (with the venv, $PYTHON is the venv's
# python, which is made from the system's). Warns rather than fails: the fix
# is a package choice for the admin.
check_mod_wsgi_python() {
    local so="" candidate built ours
    for candidate in /usr/lib/apache2/modules/mod_wsgi.so /usr/lib64/httpd/modules/mod_wsgi*.so; do
        [[ -e "$candidate" ]] && { so="$candidate"; break; }
    done
    [[ -n "$so" ]] && command -v ldd >/dev/null 2>&1 || return 0
    built="$(ldd "$so" 2>/dev/null | grep -o 'libpython[0-9]*\.[0-9]*' | head -n1 | sed 's/libpython//')"
    ours="$("$PYTHON" -c 'import sys; print("%d.%d" % sys.version_info[:2])')"
    if [[ -z "$built" ]]; then
        return 0
    elif [[ "$built" == "$ours" ]]; then
        echo "mod_wsgi runs Python $built, the same as $PYTHON."
    else
        echo "warning: mod_wsgi is built for Python $built, but the virtual environment $VENV" >&2
        echo "  is Python $ours ($PYTHON), so Apache can't use it. Either install the" >&2
        echo "  libapache2-mod-wsgi-py3 that matches Python $ours, or remake the venv from" >&2
        echo "  Python $built: sudo rm -rf $VENV && sudo PYTHON=/usr/bin/python$built ./update.sh" >&2
    fi
}

# The update guard (OPS.41, OPS.53; `planetgen.cli.update_guard`): lists the work
# that runs, stops unless the operator agrees to cancel it, closes the site
# (maintenance mode, cause "update") and takes the lock that keeps anything
# else from starting; `finish_update_guard` (update.sh's EXIT trap) opens the
# site again after a good update and leaves it closed, saying why, after a
# failed one. UPDATE_STEP names the step under way for that message.
UPDATE_GUARD_HELD=0
UPDATE_STEP="starting"

update_guard() {
    (cd / && PYTHONPATH="$SCRIPT_DIR/src${PYTHONPATH:+:$PYTHONPATH}" \
        "$VENV/bin/python" -m planetgen.cli.update_guard "$@")
}

# `update.sh --reset-jobs`: stops everything the update guard lists, kills the
# queue workers and empties the jobs table, the work lease and the Redis queues
# (`planetgen.cli.reset_jobs`), once, before the guard is taken. Like the guard
# it waits for a checkout and a virtual environment that have it.
reset_jobs_once() {
    if (( ! ${RESET_JOBS:-0} || ${RESET_JOBS_DONE:-0} )); then
        return 0
    fi
    if [[ ! -f "$SCRIPT_DIR/src/planetgen/cli/reset_jobs.py" || ! -x "$VENV/bin/python" ]]; then
        return 0
    fi
    if (cd / && PYTHONPATH="$SCRIPT_DIR/src${PYTHONPATH:+:$PYTHONPATH}" \
            "$VENV/bin/python" -m planetgen.cli.reset_jobs --yes); then
        RESET_JOBS_DONE=1
    else
        echo "error: the jobs and the queue could not be reset; nothing else was changed." >&2
        exit 1
    fi
}

# `update.sh --kill-all` (PERF.186): the full stop, `planetgen queue stop-all --no-hold` (`planetgen.queue.stopall`, the
# same routine as the Queue page's Stop everything): cancels every run and job, empties the queues, ends every planetgen
# worker and command process on the host and what they started. It runs before the guard and again after the pull
# if the checkout did not have it the first time (the first update that brings it).
kill_all_once() {
    if (( ! ${KILL_ALL:-0} || ${KILL_ALL_DONE:-0} )); then
        return 0
    fi
    if [[ ! -f "$SCRIPT_DIR/src/planetgen/cli/queue.py" || ! -x "$VENV/bin/python" ]]; then
        return 0
    fi
    if (cd / && PYTHONPATH="$SCRIPT_DIR/src${PYTHONPATH:+:$PYTHONPATH}" \
            "$VENV/bin/python" -m planetgen.cli.queue stop-all --no-hold); then
        KILL_ALL_DONE=1
    else
        echo "warning: some planetgen processes could not be ended (see above); the update goes on." >&2
        KILL_ALL_DONE=1
    fi
}

# Takes the guard if it is not held yet. Before the pull the checkout may not
# have the guard (the first update that brings it) and a new machine may have
# no virtual environment: then this does nothing, and is called again after
# the pull and after step 2. With the argument "required" it fails the update
# instead of doing nothing. Needs ANSWER_YES (an array, maybe empty).
ensure_update_guard() {
    kill_all_once
    if (( UPDATE_GUARD_HELD )); then
        return 0
    fi
    if [[ ! -f "$SCRIPT_DIR/src/planetgen/cli/update_guard.py" || ! -x "$VENV/bin/python" ]]; then
        if [[ "${1:-}" == "required" ]]; then
            echo "error: the update guard is missing from this checkout or the virtual environment is, so the update" >&2
            echo "       cannot make sure nothing else is running. Nothing was changed after the pull." >&2
            return 1
        fi
        return 0
    fi
    reset_jobs_once
    if update_guard begin --pid "$$" ${ANSWER_YES[@]+"${ANSWER_YES[@]}"}; then
        UPDATE_GUARD_HELD=1
    else
        echo "error: the update did not start." >&2
        exit 1
    fi
}

finish_update_guard() {
    local status=$?
    (( UPDATE_GUARD_HELD )) || return 0
    if (( status == 0 )); then
        update_guard end --pid "$$" || true
    else
        update_guard end --pid "$$" --failed --step "$UPDATE_STEP" || true
    fi
}

# Builds the site's caches (MAP.134, PERF.104, PERF.127): the Galaxy Map's opening view,
# the lists, every sector, system and phenomenon page (up to the warm-up cap)
# and the search, as Apache's user (the caches belong to it). It starts as soon as
# the code, the database, the permissions and Redis are in place (Boss, 10-30 EDT),
# in the background, and `wait_for_warm_opening` holds the update until its first
# phase (the opening views and the first page of each list) is built, before the
# reload and the /api/health check; the other phases go on after the update ends.
# `planetgen warm` builds the list of the warm-* commands (one per family) and
# hands it to the wrapper, which runs them in order on the Redis queue and prints
# one line per family; the update's own commands pass the update guard. At most
# 30 minutes; the output goes to $WARM_MAP_LOG. Never fails the update: the first
# visit builds what the warm-up did not.
warm_galaxy_map() {
    local user="" runner=()
    WARM_MAP_LOG="${WARM_MAP_LOG:-/tmp/planetgen-warm-map.log}"
    WARM_OPENING_FILE="${WARM_OPENING_FILE:-/tmp/planetgen-warm-opening.$$}"
    WARM_PID=""
    rm -f "$WARM_OPENING_FILE"
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
    WARM_PID="$(cd / && nohup ${runner[@]+"${runner[@]}"} ${limit[@]+"${limit[@]}"} "$PYTHON" -m planetgen.cli.generate warm \
        --opening-done "$WARM_OPENING_FILE" >"$WARM_MAP_LOG" 2>&1 </dev/null & echo $!)" || WARM_PID=""
    echo "Building the site's caches in the background (log: $WARM_MAP_LOG)."
}

# Waits (at most WARM_OPENING_WAIT_SECONDS, default 1200) until the warm-up's opening phase is built, or the warm-up
# ends, whichever comes first, so the reload and the health check come after the data the site opens with exists.
# Says what happened and never fails the update.
wait_for_warm_opening() {
    [[ -n "${WARM_PID:-}" ]] || return 0
    local end=$(( SECONDS + ${WARM_OPENING_WAIT_SECONDS:-1200} )) next=$(( SECONDS + 30 )) last
    echo "Waiting for the cache warm-up's opening phase (the opening views and the first page of each list) before reloading the site, at most $(( ${WARM_OPENING_WAIT_SECONDS:-1200} / 60 )) minutes:"
    while [[ ! -e "$WARM_OPENING_FILE" ]]; do
        if ! kill -0 "$WARM_PID" 2>/dev/null; then
            [[ -e "$WARM_OPENING_FILE" ]] || echo "warning: the cache warm-up ended before its opening phase was reported (see $WARM_MAP_LOG); the site is reloaded anyway." >&2
            return 0
        fi
        if (( SECONDS >= end )); then
            echo "warning: the cache warm-up's opening phase was not built in ${WARM_OPENING_WAIT_SECONDS:-1200} seconds; the site is reloaded anyway, the warm-up goes on." >&2
            return 0
        fi
        if (( SECONDS >= next )); then
            last="$(tail -n 1 "$WARM_MAP_LOG" 2>/dev/null | cut -c1-110)"
            echo "  still building, $(( end - SECONDS )) s left at most${last:+; last line of the warm-up log: $last}"
            next=$(( SECONDS + 30 ))
        fi
        sleep 1
    done
    echo "The opening phase is built; the rest of the warm-up goes on after the update."
}

# Tells the admin when an update changed an example config they copied
# into place. These are the admin's own files (the Apache site, the Redis
# server config): update.sh and install.sh never edit them, they only
# say that the example moved so the two can be compared.
# Args: the commit before the pull, the commit after it.
note_changed_example_configs() {
    local before="$1" after="$2"
    [[ "$before" == "$after" ]] && return 0
    if ! git diff --quiet "$before" "$after" -- examples/apache/planetgen.conf.example; then
        echo "note: examples/apache/planetgen.conf.example changed in this update."
        echo "  Your Apache site is yours and was not touched. To see what is new:"
        echo "    diff /etc/apache2/sites-available/planetgen.conf $SCRIPT_DIR/examples/apache/planetgen.conf.example"
        echo "  Copy over what you want, then: sudo apache2ctl configtest && sudo systemctl reload apache2"
    fi
    if ! git diff --quiet "$before" "$after" -- examples/redis/planetgen.conf; then
        echo "note: examples/redis/planetgen.conf changed in this update."
        echo "  Your Redis configuration was not touched. Compare it with $SCRIPT_DIR/examples/redis/planetgen.conf."
    fi
    return 0
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
    echo "Importing the web app as $user (a few seconds):"
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
    echo "Recording the version key (it runs the generator's fixed test cases once to fingerprint the code: 5 seconds to a minute, nothing prints until it ends):"
    run_cli version_history || echo "warning: could not record the version key (see above); the update carries on." >&2
}

# Works out, once and as part of the update, the mean system mass of every stored cluster age (`planetgen.cli.warm_shape`),
# which building the galaxy's density shape needs and no web request may spend the time on; the web app reads the results
# from the data cache folder. A failure only warns: the web app then works them out when it is asked.
warm_shape_means() {
    echo "Working out the cluster system masses the galaxy shape needs (a few minutes on a big galaxy, once per release):"
    run_cli warm_shape || echo "warning: could not work out the cluster system masses (see above); the update carries on." >&2
    chmod -R a+rX "${PLANETGEN_DATA_DIR:-/var/cache/planetgen/datasets}" 2>/dev/null || true
}

# The math check (`planetgen check-math`, TEST.68): known answers from
# real astronomy, identities and sampler distributions. A failure only
# warns -- the update carries on and the site keeps serving -- but sets
# MATH_CHECK_FAILED=1, so offer_population_pass skips the population pass
# and the closing message repeats the warning. Generation does not run
# the check itself (PERF.170): this is where it runs, and the admin page
# runs it on request.
MATH_CHECK_FAILED=0
check_math() {
    echo "Running the math check (a few seconds):"
    if run_cli generate check-math; then
        MATH_CHECK_FAILED=0
    else
        MATH_CHECK_FAILED=1
        echo "warning: the math check failed (above). Charted data may be wrong until it is fixed;" >&2
        echo "         run 'planetgen check-math -v' for every check." >&2
    fi
}

# Builds the data of the References pages (DOC.88): `scripts/build_reference_pages.py` writes docs/references/*.json from the
# reference list and the stand-in registry. The files are committed, so on a current checkout this changes nothing; it makes sure
# a server never shows numbers older than its own registries. A failure only warns: /references then says its pages are not built.
REFERENCE_PAGES_FAILED=0
build_reference_pages() {
    echo "Building the References pages:"
    if "$PYTHON" "$SCRIPT_DIR/scripts/build_reference_pages.py"; then
        REFERENCE_PAGES_FAILED=0
    else
        REFERENCE_PAGES_FAILED=1
        echo "warning: the References pages could not be built (above); /references may show old numbers." >&2
        echo "         run 'python3 scripts/build_reference_pages.py' to see why." >&2
    fi
}

# Optionally runs the population pass and the territories pass (`planetgen population`: species,
# civilizations and polities, then `planetgen territories`, docs/design/habitability/population-and-politics.md)
# over the stored galaxy. Off by default (Boss, 2026-10-01): it runs only
# when POPULATION=1 (or yes) is set, or when someone answers y to the
# prompt (y/N, 30 seconds, default N). With no terminal it is skipped.
# `planetgen population` and `planetgen territories` can always be run by hand later.
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
        run_cli generate territories
    else
        echo "Skipping the population pass (run 'planetgen population' any time, or set POPULATION=1)."
    fi
}
