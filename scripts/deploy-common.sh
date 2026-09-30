# scripts/deploy-common.sh
#
# Checks shared by install.sh and update.sh, sourced by both so the two
# can't drift apart. Each one looks first and only changes what is
# missing, so running it again on a working server does nothing.
# (Python libraries have their own script, install-python-deps.sh.)
#
# Expects PYTHON and SCRIPT_DIR (the repo root) to be set, and to run as
# root.

# The NLTK 'words' corpus, in a shared, world-readable directory (not a
# per-user home directory) so it works for every user that imports
# stellarObjects: a login shell running the generator, or Apache's
# www-data running the web app. stellarObjects/names.py checks
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
    else
        user=root
    fi
    if (cd / && "${runner[@]}" "$PYTHON" - "$SCRIPT_DIR" <<'EOF'
import os
import sys

root = sys.argv[1]
sys.path[:0] = [os.path.join(root, "src", "html"), os.path.join(root, "src")]
import stellarObjects  # noqa: E402,F401
from api.app import create_app  # noqa: E402,F401
import web  # noqa: E402,F401
EOF
    ); then
        echo "The web app and stellarObjects import cleanly as $user with $PYTHON."
    else
        echo "error: the web app does not import as $user with $PYTHON (see above)." >&2
        return 1
    fi
}
