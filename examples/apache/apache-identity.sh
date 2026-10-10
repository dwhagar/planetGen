#!/usr/bin/env bash
#
# examples/apache/apache-identity.sh
#
# Defines `detect_apache_group`, which prints the user and group the web
# app runs as: Apache2's workers on Linux ("www-data www-data" on
# Debian/Ubuntu), and "_www _www" on macOS, where gunicorn runs as
# macOS's built-in web server account (docs/deployment/macos.md; Homebrew's
# httpd and nginx use it too). Sourced by set-permissions.sh,
# create-cache-dir.sh and setup-debug-log.sh rather than run directly.

detect_apache_group() {
    if [[ "$(uname -s)" == Darwin ]]; then
        echo "_www _www"
        return 0
    fi

    # Debian/Ubuntu's apache2 package always defines APACHE_RUN_USER/
    # APACHE_RUN_GROUP in /etc/apache2/envvars (default www-data) -- this
    # is the authoritative source, since it's what apachectl itself reads
    # to decide who the worker processes run as.
    # APACHE_ENVVARS points the tests at a copy.
    local envvars="${APACHE_ENVVARS:-/etc/apache2/envvars}"
    if [[ -r "$envvars" ]]; then
        local user group
        # Run in a subshell so envvars' own `export`s and any unrelated
        # variables it sets don't leak into this script. `set +eu` there:
        # envvars reads $APACHE_CONFDIR unset, which under a caller's
        # `set -u` killed the subshell, and with it (silently) every
        # script that reads this function's output (TEST.62 found it).
        user="$(set +eu; source "$envvars" >/dev/null 2>&1; echo "${APACHE_RUN_USER:-}")" || user=""
        group="$(set +eu; source "$envvars" >/dev/null 2>&1; echo "${APACHE_RUN_GROUP:-}")" || group=""
        if [[ -n "$group" ]]; then
            echo "$user" "$group"
            return 0
        fi
    fi

    # Fallback: ask a running apache2 for its worker identity directly --
    # the master process runs as root, so pick the first non-root owner
    # among its children. Covers non-Debian layouts (no envvars file)
    # where apache2 is already running.
    local line
    line="$(ps -eo user:32,group:32,comm | awk '$3 == "apache2" && $1 != "root" {print $1, $2; exit}')"
    if [[ -n "$line" ]]; then
        echo "$line"
        return 0
    fi

    # Last resort: neither the config file nor a running process was
    # found -- most likely because apache2 hasn't been started/enabled
    # yet (e.g. running as part of install.sh, before the vhost exists).
    # www-data:www-data is the standard Debian/Ubuntu apache2 identity,
    # so assume it rather than hard-failing the whole install.
    echo "warning: could not detect Apache's user/group (no readable" >&2
    echo "  /etc/apache2/envvars and no running apache2 process found)." >&2
    echo "  Assuming the Debian/Ubuntu default: www-data:www-data" >&2
    echo "www-data www-data"
}
