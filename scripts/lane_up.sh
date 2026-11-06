#!/usr/bin/env bash
# Brings a build machine back after a restart: Redis, MySQL, the test environment variables and the datasets.
#
#     source scripts/lane_up.sh        starts what is down, exports the test variables into this shell, reports what is missing
#     bash scripts/lane_up.sh          the same, but the variables are only printed (a `source` keeps them)
#
# MySQL is started with --skip-log-bin (the binary logs fill the disk) and 500 connections. It takes one to three minutes after a
# restart for MySQL to recover; this waits for it. Nothing here installs packages: a missing one is reported with the command.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

export PLANETGEN_TEST_MYSQL_HOST=127.0.0.1 PLANETGEN_TEST_MYSQL_PORT=3306 PLANETGEN_TEST_MYSQL_USER=root PLANETGEN_TEST_MYSQL_PASSWORD=planetgen-test
export PLANETGEN_TEST_REDIS_URL=redis://127.0.0.1:6379/0
export PLANETGEN_MYSQL_HOST=127.0.0.1 PLANETGEN_MYSQL_PORT=3306 PLANETGEN_MYSQL_USER=root PLANETGEN_MYSQL_PASSWORD=planetgen-test
export PLANETGEN_DATA_DIR="${PLANETGEN_DATA_DIR:-/tmp/claude/data_base}"

say() { printf 'lane_up: %s\n' "$*"; }

if ! redis-cli ping >/dev/null 2>&1; then
  redis-server --daemonize yes >/dev/null 2>&1 && say "started Redis" || say "PROBLEM: Redis did not start (apt install redis-server)"
fi

mysql_alive() { mysqladmin -uroot -p"$PLANETGEN_TEST_MYSQL_PASSWORD" ping 2>/dev/null | grep -q alive; }
if ! mysql_alive; then
  if ! pgrep -x mysqld >/dev/null 2>&1; then
    (setsid nohup mysqld --user=mysql --max_connections=500 --skip-log-bin >/tmp/mysqld.out 2>&1 </dev/null &)
    say "started MySQL, waiting for it to recover"
  fi
  for _ in $(seq 1 90); do mysql_alive && break; sleep 3; done
  mysql_alive && say "MySQL is up" || say "PROBLEM: MySQL did not answer in 270 s (see /tmp/mysqld.out)"
fi

python3 - <<'PY'
import importlib.util
missing = [m for m in ("argon2", "playwright", "pytest_timeout", "xdist") if importlib.util.find_spec(m) is None]
print("lane_up: python packages ok" if not missing else f"lane_up: PROBLEM: missing {missing} (pip install -e '.[test,api]')")
PY
if [ -d "$PLANETGEN_DATA_DIR" ] && [ -n "$(ls -A "$PLANETGEN_DATA_DIR" 2>/dev/null)" ]; then
  say "datasets present in $PLANETGEN_DATA_DIR"
else
  say "PROBLEM: no datasets in $PLANETGEN_DATA_DIR (cd $ROOT/src && python -m planetgen.cli.datasets fetch)"
fi
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  say "variables for this shell:"; env | grep -E '^PLANETGEN_' | sed 's/^/export /'
fi
