#!/bin/bash
#
# setup-planetgen-env.sh -- paste this whole file into the "Setup script" field
# of a Claude Code cloud environment (Ubuntu 24.04, runs as root, before Claude
# starts). It installs everything planetGen's code and test suite need:
#
#   MySQL 8 server + client   (NOT MariaDB; root gets a throwaway password)
#   Apache 2 + mod_wsgi       (installed and configtest-checked, not serving the app)
#   Redis 7                   (preinstalled in the image; installed if missing)
#   Python dependencies       (requirements-server.lock, then the test/api/browser extras)
#   NLTK 'words' corpus       (/usr/local/share/nltk_data)
#   Playwright, pinned to the Chromium the image already carries (never downloads one)
#
# What a setup script can and cannot do (code.claude.com/docs/en/cloud-environments):
#   - It runs once, as root, and the result is a FILESYSTEM SNAPSHOT reused by later
#     sessions (rebuilt when the script changes or after about 7 days). It must exit 0
#     and finish in about 5 minutes or the snapshot is not kept.
#   - Running processes are NOT kept. So this script installs and configures, shuts
#     MySQL back down cleanly, and leaves a one-line command for each new session:
#
#         planetgen-env-up
#
#     which starts MySQL, Redis and Apache, and points the editable install at the
#     session's checkout. This script also adds a note to /root/.claude/CLAUDE.md
#     telling Claude to run it, and tries a SessionStart hook in
#     /root/.claude/settings.json so it runs by itself.
#
# Safe to run again (idempotent). No secrets: the database password below is a
# throwaway for a database that only ever exists inside the disposable VM. Nothing
# here touches your real server.
#
# Optional overrides (set as environment variables on the environment):
#   PLANETGEN_REPO          path to a checkout to read requirements from
#   PLANETGEN_REPO_URL      repo to clone if no checkout is found (default below)
#   PLANETGEN_TEST_DB_PASSWORD   MySQL root password (default planetgen-test)
#   PLAYWRIGHT_VERSION      pip playwright version (default 1.56.0 = Chromium 1194)

set -Eeuo pipefail

REPO_URL="${PLANETGEN_REPO_URL:-https://github.com/dwhagar/planetGen}"
DB_PASSWORD="${PLANETGEN_TEST_DB_PASSWORD:-planetgen-test}"
PLAYWRIGHT_VERSION="${PLAYWRIGHT_VERSION:-1.56.0}"
SEED_DIR=/opt/planetgen-seed
NLTK_DIR=/usr/local/share/nltk_data
LOG=/var/log/planetgen-setup.log
START=$(date +%s)

export DEBIAN_FRONTEND=noninteractive
export PIP_BREAK_SYSTEM_PACKAGES=1      # Ubuntu 24.04 marks its Python "externally managed"
export PIP_DISABLE_PIP_VERSION_CHECK=1
export NLTK_ALLOW_PROXIED_URLOPEN=1     # the agent proxy otherwise blocks NLTK's download

exec > >(tee -a "$LOG") 2>&1
say() { printf '\n==> [%3ss] %s\n' "$(( $(date +%s) - START ))" "$*"; }
trap 'echo "SETUP FAILED at line $LINENO (exit $?): $BASH_COMMAND" >&2; echo "Full log: $LOG" >&2' ERR

[ "$(id -u)" = 0 ] || { echo "must run as root" >&2; exit 1; }
. /etc/os-release
[ "${ID:-}" = ubuntu ] || echo "warning: written for Ubuntu 24.04, this is ${PRETTY_NAME:-unknown}" >&2
PY="$(command -v python3)"
"$PY" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 12) else 1)' \
    || { echo "planetGen needs Python 3.12 or newer, found $($PY -V)" >&2; exit 1; }

# --- 1. Find the requirements: a checkout if there is one, else a shallow clone ---
find_repo() {
    local d
    for d in "${PLANETGEN_REPO:-}" "$PWD" "${CLAUDE_PROJECT_DIR:-}" /home/user/planetGen \
             /workspace/planetGen /workspaces/planetGen /home/*/planetGen /root/planetGen "$SEED_DIR"; do
        [ -n "$d" ] && [ -f "$d/requirements-server.lock" ] && [ -f "$d/setup.py" ] && { echo "$d"; return 0; }
    done
    return 1
}
say "Locating planetGen requirements"
if ! REPO="$(find_repo)"; then
    echo "No checkout found (setup scripts may run before the clone); cloning $REPO_URL"
    rm -rf "$SEED_DIR"
    git clone --depth 1 "$REPO_URL" "$SEED_DIR"
    REPO="$SEED_DIR"
fi
echo "Using requirements from $REPO ($(git -C "$REPO" log -1 --format='%h %s' 2>/dev/null || echo 'not a git checkout'))"

# --- 2. System packages (apt), with services kept from starting inside the build ---
say "apt: MySQL 8, Apache 2 + mod_wsgi, Redis"
printf '#!/bin/sh\nexit 101\n' > /usr/sbin/policy-rc.d && chmod +x /usr/sbin/policy-rc.d
trap 'rm -f /usr/sbin/policy-rc.d' EXIT
apt-get update -qq
apt-get install -y -qq --no-install-recommends \
    mysql-server mysql-client \
    apache2 libapache2-mod-wsgi-py3 \
    redis-server redis-tools \
    ca-certificates curl git >/dev/null
[ -d /var/lib/mysql/mysql ] || mysqld --initialize-insecure --user=mysql

# --- 3. Python dependencies (parallel with the MySQL configuration below) ---
install_python() {
    say "pip: locked server requirements"
    # --ignore-installed: Debian's python3-* copies (blinker, cryptography...) have no RECORD
    # file, so pip cannot uninstall them.
    "$PY" -m pip install -q --ignore-installed -r "$REPO/requirements-server.lock"

    say "pip: test, api and browser extras (playwright pinned to the image's Chromium)"
    echo "playwright==$PLAYWRIGHT_VERSION" > /tmp/planetgen-constraints.txt
    ( cd "$REPO" && "$PY" -m pip install -q -c /tmp/planetgen-constraints.txt ".[test,api,browser]" )
    # The package itself is installed editable per session (planetgen-env-up); here only its
    # dependencies are wanted in the snapshot.
    "$PY" -m pip uninstall -y -q planetGen 2>/dev/null || true
    rm -rf "$REPO"/build "$REPO"/src/*.egg-info 2>/dev/null || true

    say "NLTK 'words' corpus"
    "$PY" -c "import nltk; ok = nltk.download('words', download_dir='$NLTK_DIR', quiet=True); raise SystemExit(0 if ok else 1)"
    [ -d "$NLTK_DIR/corpora/words" ] || [ -f "$NLTK_DIR/corpora/words.zip" ]

    if ls "${PLAYWRIGHT_BROWSERS_PATH:-/opt/pw-browsers}"/chromium-* >/dev/null 2>&1; then
        echo "Chromium present: $(ls -d "${PLAYWRIGHT_BROWSERS_PATH:-/opt/pw-browsers}"/chromium-* | xargs -n1 basename | tr '\n' ' ')"
    else
        echo "warning: no Chromium under ${PLAYWRIGHT_BROWSERS_PATH:-/opt/pw-browsers}; browser tests will skip." >&2
    fi
}
install_python > /tmp/planetgen-pip.log 2>&1 &
PIP_PID=$!

# --- 4. Configure MySQL once, then shut it down cleanly for the snapshot ---
say "MySQL: tuning and first-time root setup"
cat > /etc/mysql/mysql.conf.d/zz-planetgen-test.cnf <<'EOF'
# Throwaway test server: speed over durability. Written by setup-planetgen-env.sh.
[mysqld]
bind-address = 127.0.0.1
skip-log-bin
performance_schema = OFF
innodb_flush_log_at_trx_commit = 2
innodb_doublewrite = OFF
innodb_buffer_pool_size = 1G
max_connections = 500
EOF
mkdir -p /run/mysqld && chown mysql:mysql /run/mysqld
if ! mysqladmin ping --silent >/dev/null 2>&1; then
    nohup mysqld_safe --user=mysql >/var/log/planetgen-mysqld-setup.log 2>&1 &
    for _ in $(seq 60); do mysqladmin ping --silent >/dev/null 2>&1 && break; sleep 1; done
fi
mysqladmin ping --silent
# Root over TCP with a password (the tests connect to 127.0.0.1 as root and create/drop
# databases and one restricted user). Ubuntu's root starts on auth_socket, which ignores
# passwords, so the plugin is switched explicitly -- to MySQL 8's default
# caching_sha2_password, the same login path the real server uses. The 'planetgen' user
# (empty password, localhost only) is what conftest.py falls back to when no
# PLANETGEN_TEST_MYSQL_* variables are set at all.
SQL="ALTER USER 'root'@'localhost' IDENTIFIED WITH caching_sha2_password BY '$DB_PASSWORD';
CREATE USER IF NOT EXISTS 'planetgen'@'localhost' IDENTIFIED WITH caching_sha2_password BY '';
CREATE USER IF NOT EXISTS 'planetgen'@'127.0.0.1' IDENTIFIED WITH caching_sha2_password BY '';
GRANT ALL ON *.* TO 'planetgen'@'localhost' WITH GRANT OPTION;
GRANT ALL ON *.* TO 'planetgen'@'127.0.0.1' WITH GRANT OPTION;
FLUSH PRIVILEGES;"
mysql --no-defaults -uroot -e "$SQL" 2>/dev/null || mysql --no-defaults -uroot -p"$DB_PASSWORD" -e "$SQL"
printf '[client]\nuser=root\npassword=%s\n' "$DB_PASSWORD" > /root/.my.cnf && chmod 600 /root/.my.cnf
mysql --no-defaults -h127.0.0.1 -uroot -p"$DB_PASSWORD" -e "SELECT VERSION() AS mysql_version" 2>/dev/null
mysql --no-defaults -h127.0.0.1 -uplanetgen -e "SELECT CURRENT_USER() AS fallback_user" 2>/dev/null
mysqladmin --no-defaults -h127.0.0.1 -uroot -p"$DB_PASSWORD" shutdown 2>/dev/null || true
for _ in $(seq 30); do pgrep -x mysqld >/dev/null || break; sleep 1; done

# --- 5. Apache: modules planetGen's vhost example uses, and a config check ---
say "Apache: modules and configtest"
a2enmod -q headers deflate wsgi >/dev/null
apache2ctl configtest
a2dissite -q 000-default >/dev/null 2>&1 || true   # port 80 stays free for tests

# --- 6. The per-session command, profile variables and Claude's note ---
say "Installing planetgen-env-up and the environment variables"
cat > /usr/local/bin/planetgen-env-up <<'UP'
#!/bin/bash
# Start MySQL, Redis and Apache for this session and point the editable install of
# planetGen at the session's checkout. Safe to run any number of times.
#   planetgen-env-up            start everything, install planetGen editable
#   planetgen-env-up --quiet    same, one line of output (used by the SessionStart hook)
#   planetgen-env-up --check    only report what answers
set -uo pipefail
. /etc/profile.d/planetgen-test.sh
QUIET=0; CHECK=0
case "${1:-}" in --quiet) QUIET=1 ;; --check) CHECK=1 ;; esac
log() { [ "$QUIET" = 1 ] || echo "$*"; }
export PIP_BREAK_SYSTEM_PACKAGES=1

mysql_up()  { mysql --no-defaults -h127.0.0.1 -uroot -p"$PLANETGEN_TEST_MYSQL_PASSWORD" -e 'SELECT 1' >/dev/null 2>&1; }
redis_up()  { redis-cli -u "$PLANETGEN_TEST_REDIS_URL" ping 2>/dev/null | grep -q PONG; }
apache_up() { curl -s -o /dev/null --max-time 3 http://127.0.0.1/; }

if [ "$CHECK" = 0 ]; then
    if ! mysql_up; then
        mkdir -p /run/mysqld && chown mysql:mysql /run/mysqld
        nohup mysqld_safe --user=mysql >/var/log/planetgen-mysqld.log 2>&1 &
        for _ in $(seq 90); do mysql_up && break; sleep 1; done
    fi
    redis_up || { redis-server --daemonize yes --save "" --appendonly no >/dev/null; sleep 1; }
    apache_up || apachectl start >/dev/null 2>&1
    for d in "${PLANETGEN_REPO:-}" "${CLAUDE_PROJECT_DIR:-}" "$PWD" /home/user/planetGen /workspace/planetGen /workspaces/planetGen /root/planetGen; do
        if [ -n "$d" ] && [ -f "$d/setup.py" ] && grep -q "name='planetGen'" "$d/setup.py"; then
            SRC=$(python3 -c "import planetgen,sys; print(planetgen.__file__)" 2>/dev/null)
            if [ "${SRC#$d/}" = "$SRC" ]; then   # not already importing from this checkout
                pip install -q --no-deps --no-build-isolation -e "$d" >/dev/null 2>&1 \
                    || pip install -q --no-deps -e "$d" >/dev/null
            fi
            REPO_DIR="$d"; break
        fi
    done
fi

ok() { "$@" && echo up || echo DOWN; }
S="mysql: $(ok mysql_up) | redis: $(ok redis_up) | apache2: $(ok apache_up) | planetgen: ${REPO_DIR:-no checkout found}"
if [ "$QUIET" = 1 ]; then echo "planetgen-env-up -- $S"; else
    echo "$S"; echo "tests: $PLANETGEN_TEST_MYSQL_USER@$PLANETGEN_TEST_MYSQL_HOST:$PLANETGEN_TEST_MYSQL_PORT, $PLANETGEN_TEST_REDIS_URL"
fi
mysql_up && redis_up
UP
chmod +x /usr/local/bin/planetgen-env-up

cat > /etc/profile.d/planetgen-test.sh <<EOF
# Written by setup-planetgen-env.sh: what src/tests/conftest.py and docs/testing.md read.
export PLANETGEN_TEST_MYSQL_HOST=127.0.0.1
export PLANETGEN_TEST_MYSQL_PORT=3306
export PLANETGEN_TEST_MYSQL_USER=root
export PLANETGEN_TEST_MYSQL_PASSWORD='$DB_PASSWORD'
export PLANETGEN_TEST_REDIS_URL=redis://127.0.0.1:6379/0
export NLTK_ALLOW_PROXIED_URLOPEN=1
export PIP_BREAK_SYSTEM_PACKAGES=1
EOF
grep -q planetgen-test.sh /root/.bashrc 2>/dev/null \
    || echo '[ -f /etc/profile.d/planetgen-test.sh ] && . /etc/profile.d/planetgen-test.sh' >> /root/.bashrc

mkdir -p /root/.claude
touch /root/.claude/CLAUDE.md
if ! grep -q 'planetgen-env-up' /root/.claude/CLAUDE.md; then
    cat >> /root/.claude/CLAUDE.md <<'EOF'

## planetGen test environment (written by the environment's setup script)
MySQL 8, Redis, Apache 2, Flask and every test dependency are already installed, but a
new session starts with no services running. Before running tests, run `planetgen-env-up`
(starts MySQL, Redis and Apache and installs the checkout editable; prints one status
line; safe to repeat). The test variables (PLANETGEN_TEST_MYSQL_*, PLANETGEN_TEST_REDIS_URL)
are exported by /etc/profile.d/planetgen-test.sh. Do not run `playwright install`.
EOF
fi
# Best effort: also run it automatically at session start. Merges into an existing
# settings.json instead of replacing it.
"$PY" - <<'PYEOF'
import json, os
p = "/root/.claude/settings.json"
cfg = json.load(open(p)) if os.path.exists(p) and os.path.getsize(p) else {}
hook = {"type": "command", "command": "/usr/local/bin/planetgen-env-up --quiet", "timeout": 180}
entries = cfg.setdefault("hooks", {}).setdefault("SessionStart", [])
if not any(hook["command"] == h.get("command") for e in entries for h in e.get("hooks", [])):
    entries.append({"matcher": "startup|resume", "hooks": [hook]})
json.dump(cfg, open(p, "w"), indent=2)
PYEOF

# --- 7. Wait for pip, then prove the pieces import ---
say "Waiting for the Python install"
wait "$PIP_PID" || { echo "Python install failed; its log:" >&2; cat /tmp/planetgen-pip.log >&2; exit 1; }
cat /tmp/planetgen-pip.log
say "Verifying"
"$PY" - <<'PYEOF'
import flask, pymysql, redis, rq, numba, scipy, skimage, markdown, pytest, xdist, hypothesis, playwright, nltk
from importlib.metadata import version
print("flask", flask.__version__, "| numba", numba.__version__, "| markdown", version("markdown"),
      "| playwright", version("playwright"), "| pymysql", pymysql.__version__)
PYEOF
mysqld --version
apache2 -v | head -1
redis-server --version
echo
echo "planetGen environment ready in $(( $(date +%s) - START )) s. Start the services in each session with: planetgen-env-up"
exit 0
