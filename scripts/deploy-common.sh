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

# Every line this library prints is a log record (OPS.98): a time, a severity and a component (scripts/logrecord.sh).
# shellcheck source=scripts/logrecord.sh
declare -F log_record >/dev/null || source "$(dirname "${BASH_SOURCE[0]}")/logrecord.sh"

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
        log_warning "this host has no swap. A memory spike then makes the kernel kill MySQL, Redis or the web"
        log_warning "  server instead of slowing down. Fix: sudo $SCRIPT_DIR/scripts/enable-swap.sh  (see INSTALL.md, Swap)"
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
# From here on `update.sh` appends every line it prints to `update.log` in the log folder (OPS.101), after the lines it printed before the
# folder was known (`log_file_begin`). The folder is the code's to say (`config.json`'s `log_dir`), so the new checkout is asked; with no
# virtual environment yet, or no folder to write, the lines stay on the screen only.
start_update_log() {
    local path
    path="$(cd / && PYTHONPATH="$SCRIPT_DIR/src${PYTHONPATH:+:$PYTHONPATH}" "$VENV/bin/python" -m planetgen.util.logpaths update.log 2>/dev/null)" || return 0
    # SEC.44: the folder is config.json's `log_dir`, a file the web user can edit, and the lines are appended to it as root
    # (and the folder made): only where rootpaths.py allows it. A refused folder leaves the lines on the screen only.
    if [[ -n "$path" ]] && ! (cd / && "$VENV/bin/python" -I "$SCRIPT_DIR/examples/apache/deploy-paths.py" --check log_dir "$(dirname "$path")" "$SCRIPT_DIR" >/dev/null 2>&1); then
        log_warning "update.log was not written: log_dir in config.json is not a folder root may use; the lines stay on the screen only."
        return 0
    fi
    log_file_to "$path"
}

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
        log_info "launchd: $target present"
        return 0
    fi
    sed -e "s|/var/lib/planetGen|$SCRIPT_DIR|g" -e "s|/usr/local/planetgen/venv|$VENV|g" "$@" \
        "$example" > "$target"
    chown root:wheel "$target"
    chmod 644 "$target"
    launchctl bootstrap system "$target" 2>/dev/null || launchctl load -w "$target" || true
    log_info "launchd: installed and loaded $target ($label)"
}

# OPS.110: every systemd unit of planetGen is an example under examples/ with the default layout in it (the checkout at
# /var/lib/planetGen, its state folder at /var/lib/planetGen/state, the virtual environment at /usr/local/planetgen/venv). These
# helpers render a unit for this checkout, so the file systemd reads names the real paths. Nothing here writes outside the unit
# folder (SYSTEMD_UNIT_DIR overrides /etc/systemd/system, for tests).

# The state folder of the checkout $1, as deploy-paths.py gives it (its third line), with $2 as the Python. Read under -I so nothing
# the web user planted beside the script is imported. Prints nothing when it cannot be read.
deploy_state_folder() {
    (cd / && "$2" -I "$1/examples/apache/deploy-paths.py" "$1" 2>/dev/null | sed -n 3p) || true
}

# Succeeds when the checkout $1, the virtual environment $2 and the state folder $3 can be put into a unit file: not empty, no space
# or character that sed or systemd would read as more than text.
unit_paths_are_safe() {
    case "$1$2$3" in
        *[[:space:]\|\&\%\$\\]*|*'"'*) return 1 ;;
    esac
    [[ -n "$1" && -n "$2" && -n "$3" ]]
}

# Prints the example unit $1 as it should stand for the checkout $2, the virtual environment $3 and the state folder $4. The state
# folder is put in first (it lies under the default checkout path), then the checkout, then the environment.
rendered_unit() {
    sed -e "s|/var/lib/planetGen/state|%STATE%|g" -e "s|/var/lib/planetGen|$2|g" -e "s|/usr/local/planetgen/venv|$3|g" \
        -e "s|%STATE%|$4|g" "$1"
}

# Writes the example unit $1 rendered for the checkout $2, environment $3 and state folder $4 into the folder $5. A unit that is
# already there and equal is left alone; one that differs (an older copy of the example, or another layout) is replaced. Sets
# UNIT_ACTION to installed, replaced or unchanged; returns 1 (and writes nothing) when it could not write the file.
render_systemd_unit() {
    local source="$1" unit_dir="$5" name rendered
    name="$(basename "$source")"
    UNIT_ACTION=unchanged
    rendered="$(rendered_unit "$source" "$2" "$3" "$4")" || return 1
    if [[ -f "$unit_dir/$name" ]]; then
        [[ "$(cat "$unit_dir/$name")" == "$rendered" ]] && return 0
        UNIT_ACTION=replaced
    else
        UNIT_ACTION=installed
    fi
    if printf '%s\n' "$rendered" > "$unit_dir/$name.new" && chmod 644 "$unit_dir/$name.new" && mv -f "$unit_dir/$name.new" "$unit_dir/$name"; then
        return 0
    fi
    rm -f "$unit_dir/$name.new"
    return 1
}

# OPS.110, with SEC.43's Update button: installs every maintenance unit of planetGen (examples/maintenance) into /etc/systemd/system
# (SYSTEMD_UNIT_DIR overrides it, for tests), each rendered for this checkout, its virtual environment and its state folder, and
# enables the timers except the unattended update's, which a deploy only installs and says how to turn on. install.sh and update.sh call it with no arguments on every deploy; examples/maintenance/install-maintenance-
# timer.sh is a thin caller of this same function. The gunicorn unit (examples/systemd) is not part of it: it stays a manual copy.
#
#   the monthly orbit update       planetgen-orbits@.service/.timer (one timer instance per database)
#   the unattended update          planetgen-update.timer (and the service it and the Update button start)
#   the admin Update button        planetgen-update.path: the web user only creates <state folder>/update.request, this root-owned
#                                  unit starts planetgen-update.service, which removes the request and runs update.sh
#   the nightly map tile rebuild   planetgen-tiles@.service/.timer
#   the nightly smoke run          planetgen-smoke.service/.timer, only with --with-smoke-timer (it needs a database account that
#                                  may create and drop databases); once installed it is kept rendered like the others
#
# Options: --skip-update-timer, --skip-update-button, --skip-tiles-timer, --with-smoke-timer; --enable-all enables every requested
# timer, the update timer too (a deploy enables only the orbit, tile, Update button and smoke units it installs for the first time, so a timer an
# operator disabled stays disabled, and never the update timer);
# --credentials writes /etc/planetgen/maintenance.env (mode 600, MAINTENANCE_ENV_FILE overrides it) from the example when it is
# missing; the remaining arguments are the databases to update monthly (default $PLANETGEN_MYSQL_DATABASE, then planetgen).
# A unit already installed that differs from the rendered one (an older copy of the example, or another checkout) is replaced, and
# the replacement is logged. A checkout, environment or state path that is empty or has a space or a special character is refused
# (nothing is written). UNIT_STATE_DIR overrides the state folder deploy-paths.py gives. Linux with systemd only. Returns 1 when
# something was refused or could not be done (it has said what); a deploy only warns.
UPDATE_REQUEST_NAME=update.request
install_maintenance_units() {
    is_macos && return 0
    command -v systemctl >/dev/null 2>&1 || return 0
    local unit_dir="${SYSTEMD_UNIT_DIR:-/etc/systemd/system}" example="$SCRIPT_DIR/examples/maintenance"
    local state="${UNIT_STATE_DIR:-}" update_timer=1 update_button=1 tiles_timer=1 smoke_timer=0 enable_all=0 credentials=0
    local databases=() arg name db failed=0 changed=0 to_enable=() orbits_new=0 tiles_new=0
    for arg in "$@"; do
        case "$arg" in
            --skip-update-timer) update_timer=0 ;;
            --skip-update-button) update_button=0 ;;
            --skip-tiles-timer) tiles_timer=0 ;;
            --with-smoke-timer) smoke_timer=1 ;;
            --enable-all) enable_all=1 ;;
            --credentials) credentials=1 ;;
            *) databases+=("$arg") ;;
        esac
    done
    (( ${#databases[@]} )) || databases=("${PLANETGEN_MYSQL_DATABASE:-planetgen}")
    [[ -n "$state" ]] || state="$(deploy_state_folder "$SCRIPT_DIR" "$PYTHON")"
    if ! unit_paths_are_safe "$SCRIPT_DIR" "$VENV" "$state"; then
        log_warning "the systemd units were not installed: the checkout ($SCRIPT_DIR), the virtual environment ($VENV) or the state folder ($state) is empty or has a space or a special character, which a unit file cannot hold safely."
        return 1
    fi
    if [[ ! -d "$unit_dir" ]]; then
        log_info "no $unit_dir; the maintenance timers and the admin Update button need systemd and are not set up."
        return 0
    fi

    # Renders one unit; a unit that is not wanted but is already installed is still kept right (and left as the operator has it).
    # Sets UNIT_ACTION; a unit that is newly installed is queued for enabling.
    place_unit() {
        local unit="$1" wanted="$2" enable_as="${3:-}"
        if (( ! wanted )) && [[ ! -f "$unit_dir/$unit" ]]; then
            UNIT_ACTION=unchanged
            return 0
        fi
        if ! render_systemd_unit "$example/$unit" "$SCRIPT_DIR" "$VENV" "$state" "$unit_dir"; then
            log_warning "could not write $unit_dir/$unit."
            failed=1
            UNIT_ACTION=unchanged
            return 0
        fi
        if [[ "$UNIT_ACTION" != unchanged ]]; then
            log_info "systemd: $UNIT_ACTION $unit_dir/$unit$([[ $UNIT_ACTION == replaced ]] && echo " (it differed from the example rendered for $SCRIPT_DIR)")"
            changed=1
        fi
        if (( wanted )) && [[ -n "$enable_as" ]] && { (( enable_all )) || [[ "$UNIT_ACTION" == installed ]]; }; then
            to_enable+=("$enable_as")
        fi
        return 0
    }

    place_unit planetgen-orbits@.service 1
    place_unit planetgen-orbits@.timer 1
    [[ "$UNIT_ACTION" == installed ]] && orbits_new=1
    place_unit planetgen-update.service $(( update_timer || update_button ))
    place_unit planetgen-update.timer "$update_timer"
    (( update_timer )) && (( enable_all )) && to_enable+=(planetgen-update.timer)
    place_unit planetgen-update.path "$update_button" planetgen-update.path
    place_unit planetgen-tiles@.service "$tiles_timer"
    place_unit planetgen-tiles@.timer "$tiles_timer"
    [[ "$UNIT_ACTION" == installed ]] && tiles_new=1
    place_unit planetgen-smoke.service "$smoke_timer"
    place_unit planetgen-smoke.timer "$smoke_timer" planetgen-smoke.timer
    for db in "${databases[@]}"; do
        if (( enable_all || orbits_new )); then to_enable+=("planetgen-orbits@${db}.timer"); fi
        if (( tiles_timer )) && (( enable_all || tiles_new )); then to_enable+=("planetgen-tiles@${db}.timer"); fi
    done

    if (( credentials )); then
        local env_file="${MAINTENANCE_ENV_FILE:-/etc/planetgen/maintenance.env}"
        if [[ -f "$env_file" ]]; then
            log_info "$env_file already exists; left alone."
        elif mkdir -p "$(dirname "$env_file")" && install -m 600 "$example/maintenance.env.example" "$env_file"; then
            log_info "wrote $env_file (mode 600) from maintenance.env.example: edit it with the read-write MySQL account before the first scheduled run."
        else
            log_warning "could not write $env_file."
            failed=1
        fi
    fi
    if (( changed )); then
        systemctl daemon-reload || { log_warning "systemctl daemon-reload failed; run it by hand."; failed=1; }
    fi
    # The unattended update deploys whatever is on the tracked branch with no human in the loop, so a deploy installs the timer and never switches it
    # on (Boss, OPS.110); only the installer script, run by a person, enables it. Said at every deploy while it is off.
    if (( update_timer )) && ! (( enable_all )) && ! systemctl is-enabled --quiet planetgen-update.timer 2>/dev/null; then
        log_info "systemd: planetgen-update.timer is installed but off. When on, it runs update.sh at 03:00 on the 1st of each month, which pulls and deploys whatever is on the tracked branch with no review. To turn it on: sudo systemctl enable --now planetgen-update.timer"
    fi
    for name in ${to_enable[@]+"${to_enable[@]}"}; do
        if systemctl enable --now "$name" >/dev/null 2>&1; then
            log_info "systemd: enabled $name"
        else
            log_warning "could not enable $name; run: sudo systemctl enable --now $name"
            failed=1
        fi
    done
    return "$failed"
}

# macOS's counterpart of the Apache module check: gunicorn under launchd
# (examples/macos/org.planetgen.gunicorn.plist), serving on 127.0.0.1:8000
# for nginx. Installed and started the first time only. Its log
# directory belongs to _www.
ensure_gunicorn_daemon() {
    # gunicorn's own output; a folder that can't be made only warns (OPS.5).
    if ! { mkdir -p "$MACOS_LOG_DIR" && chown _www:_www "$MACOS_LOG_DIR"; }; then
        log_warning "couldn't set up gunicorn's log folder $MACOS_LOG_DIR; launchd can't start"
        log_warning "  gunicorn until it exists. Fix it with:"
        log_warning "    sudo mkdir -p $MACOS_LOG_DIR"
        log_warning "    sudo chown _www:_www $MACOS_LOG_DIR"
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
        log_info "NLTK 'words' corpus: present in $dir"
    else
        # A plain nltk.download() rather than `python -m nltk.downloader`,
        # which prints a spurious "found in sys.modules" RuntimeWarning
        # (nltk's __init__ already imported nltk.downloader).
        "$PYTHON" -c "import nltk, sys; sys.exit(0 if nltk.download('words', download_dir=sys.argv[1]) else 1)" "$dir"
        log_info "NLTK 'words' corpus: installed in $dir"
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
        log_warning "some data could not be downloaded (above); bulk generation refuses to start until the next"
        log_warning "         update fetches it ('python3 -m planetgen.cli.datasets status' says which)."
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
        log_warning "couldn't read redis.url from config.json; skipping the Redis check."
        return 0
    }
    host="$("$PYTHON" -c "import sys; from urllib.parse import urlsplit; print(urlsplit(sys.argv[1]).hostname or '')" "$url")"
    if redis_answers "$url"; then
        log_info "Redis: answering at $url"
        return 0
    fi
    if [[ "$host" == 127.0.0.1 || "$host" == localhost || "$host" == ::1 ]] && ! is_macos \
            && command -v apt-get >/dev/null 2>&1; then
        if ! command -v redis-server >/dev/null 2>&1; then
            log_info "Redis: installing redis-server."
            DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends redis-server || true
        fi
        if command -v systemctl >/dev/null 2>&1; then
            systemctl enable --now redis-server >/dev/null 2>&1 || true
        fi
        if redis_answers "$url"; then
            log_info "Redis: answering at $url"
            return 0
        fi
    fi
    log_warning "no Redis server answers at $url (config.json's redis.url)."
    if is_macos; then
        log_warning "  Install and start it as your own user: brew install redis && brew services start redis"
    else
        log_warning "  Install and start it: sudo apt install redis-server && sudo systemctl enable --now redis-server"
    fi
    log_warning "  Without it the Generate page starts no job and a run cannot use more than one worker (docs/deployment/README.md)."
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
        log_warning "a2enmod not found -- is apache2 installed?"
        log_warning "  Try: sudo apt install apache2"
        return 0
    fi
    if [[ ! -e /etc/apache2/mods-available/wsgi.load ]] && command -v apt-get >/dev/null 2>&1; then
        log_info "mod_wsgi is not installed: installing libapache2-mod-wsgi-py3."
        # The package enables the module itself, so a2query below finds
        # it already on; Apache still has to be restarted to load it.
        DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends libapache2-mod-wsgi-py3 \
            && APACHE_NEEDS_RESTART=1 || true
    fi
    local mod enabled=()
    for mod in headers deflate wsgi; do
        if a2query -q -m "$mod" 2>/dev/null; then
            log_info "Apache module $mod: enabled"
        elif a2enmod -q "$mod" >/dev/null 2>&1; then
            log_info "Apache module $mod: enabled now"
            enabled+=("$mod")
        else
            log_warning "could not enable Apache module $mod."
            if [[ "$mod" == wsgi ]]; then
                log_warning "  Install it first: sudo apt install libapache2-mod-wsgi-py3 && sudo a2enmod wsgi"
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
# names another address; HEALTH_WAIT_SECONDS, default 90, is how long it waits).
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
                log_info "reloaded gunicorn."
                acted=1
                wait_for_health "$health_url" gunicorn "$health_from" || unhealthy=1
            else
                log_warning "launchctl kill SIGHUP system/$GUNICORN_LABEL failed."
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
            log_info "Testing Apache's configuration, then ${verb}ing it (a request still running is cut off):"
            if ! apache_configtest; then
                log_warning "Apache's configuration test failed, so Apache was not ${verb}ed (the running one keeps serving)."
                failed=1
            elif systemctl "$verb" "$unit" >/dev/null 2>&1; then
                log_info "${verb}ed Apache."
                wait_for_health "$health_url" Apache "$health_from" || unhealthy=1
            else
                log_warning "systemctl $verb $unit failed."
                failed=1
            fi
        fi
        if systemctl is-active --quiet planetgen-gunicorn 2>/dev/null; then
            acted=1
            if systemctl reload planetgen-gunicorn >/dev/null 2>&1; then
                log_info "reloaded gunicorn."
                wait_for_health "$health_url" gunicorn "$health_from" || unhealthy=1
            else
                log_warning "systemctl reload planetgen-gunicorn failed."
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
        log_info "  sudo launchctl kill SIGHUP system/$GUNICORN_LABEL"
    else
        log_info "  sudo apache2ctl configtest && sudo systemctl reload apache2     (Apache)"
        log_info "  sudo systemctl reload planetgen-gunicorn                        (nginx or Caddy hosts)"
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
        log_error "no API address is configured: set api_base_url in config.json (for example \"https://your-site/api\", the address Apache serves the API at;"
        log_error "  PLANETGEN_API_BASE_URL in the environment does the same), or PLANETGEN_HEALTH_URL for the health check alone."
        log_error "  The update asks only that address whether the site answers and never tries 127.0.0.1 or localhost unless it is the one configured."
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

# Waits up to HEALTH_WAIT_SECONDS (default 90, OPS.94: workers still compiling and importing after new code can need over 20 s) for $1
# (an address) to answer below 400, naming $2 in the warning and $3 as the setting the address came from (OPS.86). It asks that one
# address and no other, and stops at the first good answer. Every retry line says what curl saw (a timeout, a refused connection, a
# failed name lookup and a TLS error read differently). A passing answer while the site is closed for the update says so (and why), on
# success and in the warning alike: /api/health stays open while the site is closed (OPS.75), so a reload is checked while the
# update's own maintenance mode is still on.
# Returns 0 when it answers, 1 when it never did (or without curl, 0).
wait_for_health() {
    command -v curl >/dev/null 2>&1 || return 0
    # The reload has just restarted the workers: they get HEALTH_SETTLE_SECONDS (default 10) to start before the first probe.
    log_info "Giving the workers ${HEALTH_SETTLE_SECONDS:-10} seconds to start, then asking $1 (${3:-its setting}) for /api/health (up to ${HEALTH_WAIT_SECONDS:-90} seconds):"
    sleep "${HEALTH_SETTLE_SECONDS:-10}"
    local end=$(( SECONDS + ${HEALTH_WAIT_SECONDS:-90} )) reply code body rc why last_said=-100 last_why=""
    while :; do
        rc=0
        reply="$(curl -s --max-time 5 -w '\n%{http_code}' "$1" 2>/dev/null)" || rc=$?
        code="${reply##*$'\n'}"
        body="${reply%$'\n'*}"
        if [[ "$code" =~ ^[23] ]]; then
            health_maintenance_note "$body"
            return 0
        fi
        why="HTTP ${code:-000}"
        [[ "${code:-000}" == 000 ]] && why+=": $(curl_failure_words "$rc")"
        if (( SECONDS >= end )); then
            local bare_note=""
            [[ "${code:-000}" == 000 ]] && bare_note=" (no HTTP answer at all: $(curl_failure_words "$rc"))"
            log_warning "the health check of $2 failed: $1 (${3:-its setting}) answered HTTP ${code:-000}${bare_note} for ${HEALTH_WAIT_SECONDS:-90} seconds after the reload.$(health_maintenance_note "$body" brief)"
            log_warning "  Only that address was asked; check ${3:-the setting} if it is not where the site is served."
            return 1
        fi
        # One line for the first miss, then one every 5 seconds or whenever the reason changes.
        if (( SECONDS - last_said >= 5 )) || [[ "$why" != "$last_why" ]]; then
            log_info "  no good answer yet ($why); asking again, $(( end - SECONDS )) s left"
            last_said=$SECONDS last_why="$why"
        fi
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
        log_warning "Apache's user $user can't read ${blocked:-the venv}, so mod_wsgi can't use $VENV."
        log_warning "  Every directory down to the venv needs to be enterable by it (for example:"
        log_warning "  sudo chmod o+rx ${blocked:-$VENV}), and the venv readable. Then restart Apache."
    else
        log_info "Apache's user $user can read $VENV."
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
        log_warning "$conf runs planetGen in Apache's embedded mode, which can't use a virtual environment."
        log_warning "  Declare a daemon process and use it (examples/apache/planetgen.conf.example):"
        log_warning "    WSGIDaemonProcess planetgen-api processes=1 threads=5 request-timeout=60 python-home=$VENV"
        log_warning "    WSGIScriptAlias / $SCRIPT_DIR/src/html/wsgi.py process-group=planetgen-api application-group=%{GLOBAL}"
    else
        daemon="$(grep -E "^[[:space:]]*WSGIDaemonProcess[[:space:]]+$group([[:space:]]|\$)" <<< "$active" | head -n1 || true)"
        if [[ -z "$daemon" ]]; then
            log_warning "$conf names the process group '$group' but has no WSGIDaemonProcess $group line."
            log_warning "  Add: WSGIDaemonProcess $group processes=1 threads=5 request-timeout=60 python-home=$VENV"
        else
            home="$(sed -n 's/.*[[:space:]]python-home=\([^[:space:]]*\).*/\1/p' <<< "$daemon")"
            [[ "$home" == / ]] || home="${home%/}"
            if [[ "$home" == "$VENV" ]]; then
                log_info "Apache runs the site under $VENV (process group $group)."
            elif [[ -z "$home" ]]; then
                log_warning "$conf does not run the site under planetGen's virtual environment."
                log_warning "  Add python-home=$VENV to its WSGIDaemonProcess $group line (examples/apache/planetgen.conf.example),"
                log_warning "  then restart Apache."
            elif [[ "$home" == "$VENV"/bin || "$home" == "$VENV"/bin/* ]]; then
                log_warning "$conf has python-home=$home, a path into the venv's bin/."
                log_warning "  python-home is the venv's root, the folder that holds bin/ and lib/: python-home=$VENV"
            else
                log_warning "$conf has python-home=$home, not planetGen's virtual environment $VENV."
                log_warning "  Set it to python-home=$VENV (the folder that holds bin/ and lib/), then restart Apache."
            fi
        fi
        if ! grep -qF 'application-group=%{GLOBAL}' <<< "$alias_line" \
                && ! grep -qE '^[[:space:]]*WSGIApplicationGroup[[:space:]]+%\{GLOBAL\}' <<< "$active"; then
            log_warning "$conf doesn't run planetGen in the main Python interpreter."
            log_warning "  Add application-group=%{GLOBAL} to its WSGIScriptAlias line (or WSGIApplicationGroup %{GLOBAL}),"
            log_warning "  as numpy, scipy and numba's C code need it."
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
        log_info "mod_wsgi runs Python $built, the same as $PYTHON."
    else
        log_warning "mod_wsgi is built for Python $built, but the virtual environment $VENV"
        log_warning "  is Python $ours ($PYTHON), so Apache can't use it. Either install the"
        log_warning "  libapache2-mod-wsgi-py3 that matches Python $ours, or remake the venv from"
        log_warning "  Python $built: sudo rm -rf $VENV && sudo PYTHON=/usr/bin/python$built ./update.sh"
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
        log_error "the jobs and the queue could not be reset; nothing else was changed."
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
        log_warning "some planetgen processes could not be ended (see above); the update goes on."
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
            log_error "the update guard is missing from this checkout or the virtual environment is, so the update"
            log_error "       cannot make sure nothing else is running. Nothing was changed after the pull."
            return 1
        fi
        return 0
    fi
    reset_jobs_once
    if update_guard begin --pid "$$" ${ANSWER_YES[@]+"${ANSWER_YES[@]}"}; then
        UPDATE_GUARD_HELD=1
    else
        log_error "the update did not start."
        exit 1
    fi
}

finish_update_guard() {
    local status=$?
    log_summary_forget
    (( UPDATE_GUARD_HELD )) || return 0
    if (( status == 0 )); then
        update_guard end --pid "$$" || true
    else
        update_guard end --pid "$$" --failed --step "$UPDATE_STEP" || true
    fi
}

# Builds the site's caches (MAP.134, PERF.104, PERF.127) in two stages (Boss, 11-01 EDT), as Apache's user (the caches
# belong to it), in the foreground with a progress bar (`planetgen warm --progress`; the update's own commands pass the
# update guard). Stage 1, `opening` (the galaxy's cluster field, the object index, the opening views and the first page
# of each list), is just what the Apache reload and the /api/health check need and runs BEFORE the reload, so the reload
# cannot interrupt it; it starts no other phase. Stage 2, the full warm-up (the opening phase is already stored), runs
# AFTER the reload and the health check, and update.sh stays up showing its bar until it ends. Args: the stage's
# name, then what `planetgen warm` is asked to build (nothing: every phase). Each stage has a time limit,
# WARM_STAGE_SECONDS (default 4 hours, OPS.95: a first warm-up after a new epoch takes far longer than the 30 minutes it
# once had; 0 means none), told in the stage's first line. The run writes its own record to warmup.log in the log folder (OPS.96: `planetgen warm` appends it
# itself, in local time with its zone, so nothing here pipes the output through tee and the bars draw on a terminal). A stage that fails or
# runs out of time only warns, and says that what it stored is kept and how to continue: the first visit builds what it
# did not.
WARM_STAGE_DEFAULT_SECONDS=14400

# The key that skips the cache warm-up while it runs (only with a terminal), and whether it was pressed: a skip ends the stage
# that runs and leaves out the stages after it.
WARM_SKIP_KEY="s"
WARM_SKIP_KEY_UPPER="S"    # both cases are spelled out: upper-casing by expansion is bash 4 only, and this script runs under macOS's 3.2
WARM_SKIPPED=0

# Sends SIGTERM to process $1 and all of its descendants (the warm-up runs under runuser, env and timeout).
kill_process_tree() {
    local child
    for child in $(pgrep -P "$1" 2>/dev/null); do
        kill_process_tree "$child"
    done
    kill -TERM "$1" 2>/dev/null || true
}

# A number of seconds in words: 14400 is "4 hours", 5400 is "90 minutes", 45 is "45 seconds".
duration_words() {
    local seconds="$1" count
    if (( seconds >= 3600 && seconds % 3600 == 0 )); then
        count=$(( seconds / 3600 ))
        echo "$count hour$(( count == 1 ? 0 : 1 ))" | sed 's/hour0/hour/; s/hour1/hours/'
    elif (( seconds >= 60 && seconds % 60 == 0 )); then
        count=$(( seconds / 60 ))
        echo "$count minute$(( count == 1 ? 0 : 1 ))" | sed 's/minute0/minute/; s/minute1/minutes/'
    else
        echo "$seconds second$(( seconds == 1 ? 0 : 1 ))" | sed 's/second0/second/; s/second1/seconds/'
    fi
}

run_warm_stage() {
    local name user="" runner=() limit=() status=0 seconds="${WARM_STAGE_SECONDS:-$WARM_STAGE_DEFAULT_SECONDS}" limit_words warm_started=$SECONDS
    name="$(log_lower "$(log_count Part "$1" "$2")"), $3"
    shift 3
    if (( WARM_SKIPPED )); then
        log_info "Skipping the cache warm-up, $name (skipped by a keypress)."
        return 0
    fi
    if ! [[ "$seconds" =~ ^[0-9]+$ ]]; then
        log_warning "WARM_STAGE_SECONDS=$seconds is not a number of seconds; using $WARM_STAGE_DEFAULT_SECONDS."
        seconds=$WARM_STAGE_DEFAULT_SECONDS
    fi
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
    if (( seconds > 0 )) && command -v timeout >/dev/null 2>&1; then
        limit=(timeout "$seconds" env "PLANETGEN_WARM_TIME_LIMIT_SECONDS=$seconds")
        limit_words="time limit $(duration_words "$seconds"), WARM_STAGE_SECONDS sets it"
    else
        limit_words="no time limit"
    fi
    log_info "Building the site's caches, $name ($limit_words; log: warmup.log in the log folder):"
    local can_skip=0
    [[ -t 0 && -t 1 ]] && can_skip=1
    (( can_skip )) && log_info "Press $WARM_SKIP_KEY to skip the cache warm-up (what it has stored is kept; the first visit builds the rest)."
    (cd / && ${runner[@]+"${runner[@]}"} env "PLANETGEN_WARM_STARTED_BY=update.sh, $name" ${limit[@]+"${limit[@]}"} \
        "$PYTHON" -m planetgen.cli.generate warm "$@" --progress </dev/null 2>&1) &
    local pid=$! key="" read_status=0
    if (( can_skip )); then
        while kill -0 "$pid" 2>/dev/null; do
            read_status=0
            read -rsn1 -t 1 key || read_status=$?
            if (( read_status == 0 )); then
                if [[ "$key" == "$WARM_SKIP_KEY" || "$key" == "$WARM_SKIP_KEY_UPPER" ]]; then
                    WARM_SKIPPED=1
                    kill_process_tree "$pid"
                    break
                fi
            elif (( read_status <= 128 )); then
                sleep 1         # end of input: nothing to read, do not spin
            fi
        done
    fi
    wait "$pid" || status=$?
    if (( WARM_SKIPPED )); then
        log_info "The cache warm-up ($name) was skipped by a keypress. What it stored is kept."
        log_info "  'planetgen warm --progress' (as the web server's user) or the next ./update.sh continues it, and the first visit builds what it did not."
        WARM_RESULT="cache warm-up ($name) skipped by a keypress"
        return 0
    fi
    if (( status != 0 )); then
        local why="failed with status $status (see warmup.log in the log folder)"
        (( status == 124 )) && why="ran out of its time limit of $(duration_words "$seconds") (WARM_STAGE_SECONDS; 0 means none)"
        log_warning "the cache warm-up ($name) did not finish: it $why."
        log_warning "  What it stored is kept. 'planetgen warm --progress' (as the web server's user) or the next ./update.sh continues where it stopped, and the first visit builds what it did not."
        WARM_RESULT="cache warm-up ($name) did not finish"
    else
        WARM_RESULT="cache warm-up ($name) done in $(log_duration $(( SECONDS - warm_started )))"
    fi
    return 0
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
        log_notice "examples/apache/planetgen.conf.example changed in this update."
        log_notice "  Your Apache site is yours and was not touched. To see what is new:"
        log_notice "    diff /etc/apache2/sites-available/planetgen.conf $SCRIPT_DIR/examples/apache/planetgen.conf.example"
        log_notice "  Copy over what you want, then: sudo apache2ctl configtest && sudo systemctl reload apache2"
    fi
    if ! git diff --quiet "$before" "$after" -- examples/redis/planetgen.conf; then
        log_notice "examples/redis/planetgen.conf changed in this update."
        log_notice "  Your Redis configuration was not touched. Compare it with $SCRIPT_DIR/examples/redis/planetgen.conf."
    fi
    return 0
}

# Imports the web app and the generator package with the interpreter the
# site runs under, as Apache's own user when that user exists, so a
# library that's installed but unreadable to www-data (or the checkout's
# own code failing to import) shows up here rather than as a 500.
set_secret_key() {
    local user="" runner=()
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
    (cd / && ${runner[@]+"${runner[@]}"} "$PYTHON" -m planetgen.util.secret_key) \
        || log_warning "the secret key could not be set (see above); the web app falls back to a random key per process."
}

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
    log_info "Importing the web app as $user:"
    if (cd / && ${runner[@]+"${runner[@]}"} "$PYTHON" - <<'EOF'
import planetgen.generation.system  # noqa: F401
from planetgen.web.app import create_app  # noqa: F401
EOF
    ); then
        log_info "The web app and the generator import cleanly as $user with $PYTHON."
    else
        log_error "the web app does not import as $user with $PYTHON (see above)."
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
    DB_SCHEMA_VERSION="$target"
    if (( pending > 0 )); then
        log_info "Database '$database' is at schema v$current; this version needs v$target ($pending migration step(s))."
        if [[ -t 0 ]]; then
            read -r -t 30 -p "Delete all galaxy data in '$database' instead of migrating it? [y/N] (default N in 30s): " answer \
                || { echo; answer=""; }
        else
            log_info "(No terminal to ask on: keeping the data and migrating it.)"
        fi
        if [[ "$answer" =~ ^[Yy]([Ee][Ss])?$ ]]; then
            log_info "Deleting the galaxy data in '$database' (admin logins are kept)."
            run_cli reset --yes
        else
            log_info "Keeping the data and migrating it."
        fi
    fi
    run_cli migrate
}

# Records the version key this update runs under, one row per galaxy in the
# control database's history (OPS.13; the last 10 are kept), after the
# schemas are current. A failure only warns: the update carries on.
record_version_key() {
    log_info "Recording the version key (it runs the generator's fixed test cases once to fingerprint the code):"
    run_cli version_history || log_warning "could not record the version key (see above); the update carries on."
}

# Works out, once and as part of the update, the mean system mass of every stored cluster age (`planetgen.cli.warm_shape`),
# which building the galaxy's density shape needs and no web request may spend the time on; the web app reads the results
# from the data cache folder. A failure only warns: the web app then works them out when it is asked.
warm_shape_means() {
    log_info "Working out the cluster system masses the galaxy shape needs:"
    run_cli warm_shape || log_warning "could not work out the cluster system masses (see above); the update carries on."
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
    log_info "Running the math check:"
    if run_cli generate check-math; then
        MATH_CHECK_FAILED=0
    else
        MATH_CHECK_FAILED=1
        log_warning "the math check failed (above). Charted data may be wrong until it is fixed;"
        log_warning "         run 'planetgen check-math -v' for every check."
    fi
}

# Builds the data of the References pages (DOC.88): `scripts/build_reference_pages.py` writes docs/references/*.json from the
# reference list and the stand-in registry. The files are committed, so on a current checkout this changes nothing; it makes sure
# a server never shows numbers older than its own registries. A failure only warns: /references then says its pages are not built.
REFERENCE_PAGES_FAILED=0
build_reference_pages() {
    log_info "Building the References pages:"
    if "$PYTHON" "$SCRIPT_DIR/scripts/build_reference_pages.py"; then
        REFERENCE_PAGES_FAILED=0
    else
        REFERENCE_PAGES_FAILED=1
        log_warning "the References pages could not be built (above); /references may show old numbers."
        log_warning "         run 'python3 scripts/build_reference_pages.py' to see why."
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
        log_info "Skipping the population pass: the math check failed."
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
        log_info "Running the population pass."
        run_cli generate population
        run_cli generate territories
    else
        log_info "Skipping the population pass (run 'planetgen population' any time, or set POPULATION=1)."
    fi
}
