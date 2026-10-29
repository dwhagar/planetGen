# tests/lane_services.py

"""
What a test run stands on, checked before it starts and watched while it runs (TEST.376).

A lane machine is a container that restarts, and when it does MySQL and Redis go with it; the MySQL binlogs fill the
disk; a fresh shell may not have the `PLANETGEN_TEST_*` variables. Each of those cost lanes hours of gate on
2026-10-24: a database test whose server is down is *skipped* (`conftest._mysql_server_available`), a run started
without `PLANETGEN_TEST_MYSQL_USER` logs in as `planetgen` and fails with `DatabaseLoginError`, and a server that dies in
the middle fails every test that runs until someone notices.

`scripts/test.py` therefore calls:

- `preflight(environ)` before a run: loads the `PLANETGEN_TEST_*` variables from `/etc/profile.d/planetgen-test.sh`
  when they are missing, starts MySQL and Redis when they do not answer (with the commands `planetgen-env-up` runs;
  not that script itself, which also starts Apache), purges the MySQL binlogs when the disk is nearly full, and checks the
  packages the tests import (`argon2-cffi`) and the one browser package that must match the installed Chromium
  (`playwright`, the version `examples/cloud-environment/setup-planetgen-env.sh` pins), and downloads the data the tests
  read as it is (`TEST_DATASETS`) when the machine has not (TEST.441). What it could not repair comes
  back as problems, and the run does not start.
- `Watchdog` during a run: a thread that asks MySQL and Redis every `INTERVAL` seconds. A server that does not answer is
  restarted at once and the watchdog tells the runner, which stops the phase's pytest, and starts it again without the
  tests that already passed (the progress log, `--skip-passed`), so the outage costs the tests that were in flight and
  not the run.

Every command and the service names can be replaced through `PLANETGEN_TEST_MYSQL_START` and
`PLANETGEN_TEST_REDIS_START` (a shell command each) on a machine set up another way.
"""

import importlib.metadata
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from urllib.parse import urlparse

PROFILE = "/etc/profile.d/planetgen-test.sh"
PLAYWRIGHT_VERSION = "1.56.0"
"""str: The Playwright the lane image's Chromium matches (`PLAYWRIGHT_VERSION` in setup-planetgen-env.sh)."""
REQUIRED = ("PLANETGEN_TEST_MYSQL_USER", "PLANETGEN_TEST_MYSQL_PASSWORD", "PLANETGEN_TEST_REDIS_URL")
"""tuple: Variables a run needs. A missing user is not an error to the tests, it is the wrong account (`planetgen`)."""
INTERVAL = 15
"""int: Seconds between two questions of the watchdog."""
HUNG_AFTER = 3
"""int: Answers in a row that must time out before a server that still holds its port counts as dead; a refused
connection counts at once, since nothing is listening."""
START_WAIT = 90
"""int: Seconds a restarted server has to answer."""
TEST_DATASETS = ("white_dwarf_models",)
"""tuple: The downloaded data the tests read as it is, not through a stand-in: the white-dwarf cooling sequences (a stand-in gives other
stars, so the same-seed and migration tests disagree with a machine that has them). `planetgen.cli.datasets fetch` is the setup step of
`update.sh`; a fresh test machine has not run it (TEST.441)."""
LOW_DISK_MB = 3000
"""int: Free space under which the MySQL binlogs are purged."""

_EXPORT = re.compile(r"""^\s*export\s+([A-Z_][A-Z0-9_]*)=(.*?)\s*$""")


def load_profile(environ, path=PROFILE):
    """Sets the variables `path` exports that `environ` lacks; returns their names. The file is the shell profile the
    cloud image writes: `export NAME=value` lines, a value optionally in single or double quotes."""
    loaded = []
    try:
        with open(path, encoding="utf-8") as profile:
            lines = profile.read().splitlines()
    except OSError:
        return loaded
    for line in lines:
        found = _EXPORT.match(line)
        if not found:
            continue
        name, value = found.groups()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in "'\"":
            value = value[1:-1]
        if not environ.get(name):
            environ[name] = value
            loaded.append(name)
    return loaded


def mysql_state(environ):
    """`"up"`, `"down"` (nothing listens or the login is refused outright) or `"hung"` (it holds the port and does
    not answer in time)."""
    import pymysql
    try:
        conn = pymysql.connect(
            host=environ.get("PLANETGEN_TEST_MYSQL_HOST", "127.0.0.1"),
            port=int(environ.get("PLANETGEN_TEST_MYSQL_PORT", 3306)),
            user=environ.get("PLANETGEN_TEST_MYSQL_USER", "planetgen"),
            password=environ.get("PLANETGEN_TEST_MYSQL_PASSWORD", ""), connect_timeout=10, read_timeout=10)
    except pymysql.err.OperationalError as exc:
        code = exc.args[0] if exc.args else 0
        if code == 1045:        # access denied: the server is up, the account is wrong (reported, not restarted)
            return "denied"
        if 1000 <= code < 2000:  # any other answer of the server itself (too many connections, ...): it is alive
            return "up"
        return "hung" if code == 2003 and "timed out" in str(exc) else "down"
    except (OSError, pymysql.MySQLError):
        return "down"
    try:
        with conn.cursor() as cursor:
            cursor.execute("SELECT 1")
        return "up"
    except (OSError, pymysql.MySQLError):
        return "hung"
    finally:
        conn.close()


def redis_state(environ):
    """`"up"` or `"down"` for the Redis server `PLANETGEN_TEST_REDIS_URL` names."""
    import redis
    try:
        return "up" if redis.Redis.from_url(environ["PLANETGEN_TEST_REDIS_URL"], socket_connect_timeout=5,
                                            socket_timeout=5).ping() else "down"
    except (redis.RedisError, OSError, KeyError):
        return "down"


def _run(command, **options):
    return subprocess.run(command, shell=isinstance(command, str), stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL, check=False, **options)


def _detached(command):
    subprocess.Popen(command, shell=isinstance(command, str), stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                     stderr=subprocess.DEVNULL, start_new_session=True)


def _wait_for(state, environ, seconds):
    deadline = time.time() + seconds
    while time.time() < deadline:
        if state(environ) == "up":
            return True
        time.sleep(1)
    return state(environ) == "up"


def _running(*names):
    """Whether a process with one of these names (the kernel's `comm`) exists."""
    try:
        for entry in os.listdir("/proc"):
            if entry.isdigit():
                try:
                    with open(f"/proc/{entry}/comm", encoding="utf-8") as comm:
                        if comm.read().strip() in names:
                            return True
                except OSError:
                    pass
    except OSError:
        pass
    return False


def start_mysql(environ):
    """Starts the MySQL server; True when it answers. `PLANETGEN_TEST_MYSQL_START`, else
    `mysqld_safe --user=mysql` (what `planetgen-env-up` runs). A server process that still exists (`mysqld_safe` brings a crashed
    `mysqld` back by itself, and a server in crash recovery answers late) is waited for, not started over."""
    custom = environ.get("PLANETGEN_TEST_MYSQL_START")
    if _running("mysqld_safe", "mysqld", "mariadbd"):
        pass
    elif custom:
        _detached(custom)
    elif shutil.which("mysqld_safe"):
        os.makedirs("/run/mysqld", exist_ok=True)
        _run(["chown", "mysql:mysql", "/run/mysqld"])
        _detached(["mysqld_safe", "--user=mysql"])
    else:
        return False
    return _wait_for(mysql_state, environ, START_WAIT)


def start_redis(environ):
    """Starts the Redis server; True when it answers (`PLANETGEN_TEST_REDIS_START`, else
    `redis-server --port <the URL's> --daemonize yes --save "" --appendonly no`)."""
    custom = environ.get("PLANETGEN_TEST_REDIS_START")
    if custom:
        _detached(custom)
    elif shutil.which("redis-server"):
        port = urlparse(environ.get("PLANETGEN_TEST_REDIS_URL", "")).port or 6379
        _run(["redis-server", "--port", str(port), "--daemonize", "yes", "--save", "", "--appendonly", "no"])
    else:
        return False
    return _wait_for(redis_state, environ, 30)


def purge_binlogs(environ):
    """`PURGE BINARY LOGS BEFORE NOW()` (the binlogs of a test server are of no use and fill the disk); True when it ran."""
    import pymysql
    try:
        conn = pymysql.connect(
            host=environ.get("PLANETGEN_TEST_MYSQL_HOST", "127.0.0.1"),
            port=int(environ.get("PLANETGEN_TEST_MYSQL_PORT", 3306)),
            user=environ.get("PLANETGEN_TEST_MYSQL_USER", "planetgen"),
            password=environ.get("PLANETGEN_TEST_MYSQL_PASSWORD", ""), connect_timeout=10)
    except pymysql.MySQLError:
        return False
    try:
        with conn.cursor() as cursor:
            cursor.execute("PURGE BINARY LOGS BEFORE NOW()")
        return True
    except pymysql.MySQLError:
        return False
    finally:
        conn.close()


def free_mb(path="/"):
    usage = shutil.disk_usage(path)
    return usage.free // (1024 * 1024)


def package_problems(browser):
    """What the tests import that is missing or wrong, as sentences with the command that fixes it."""
    problems = []
    try:
        importlib.metadata.version("argon2-cffi")
    except importlib.metadata.PackageNotFoundError:
        problems.append("argon2-cffi is not installed (pip install argon2-cffi); the password tests cannot run")
    if browser:
        try:
            have = importlib.metadata.version("playwright")
        except importlib.metadata.PackageNotFoundError:
            have = None
        if have != PLAYWRIGHT_VERSION:
            problems.append(f"playwright is {have or 'not installed'}, the image's Chromium needs {PLAYWRIGHT_VERSION} "
                            f"(pip install playwright=={PLAYWRIGHT_VERSION}; do not run `playwright install`)")
    return problems


def dataset_problems(say=print, run=None):
    """
    Fetches the `TEST_DATASETS` a machine lacks (`planetgen.cli.datasets status` exits 1 for a missing one, `fetch` downloads it).

    Returns:
        list: One sentence per dataset that is still missing after the fetch.
    """
    run = run or _run
    problems = []
    for key in TEST_DATASETS:
        command = [sys.executable, "-m", "planetgen.cli.datasets"]
        if run([*command, "status", key]).returncode == 0:
            continue
        say(f"preflight: the {key} data is not downloaded; fetching it")
        run([*command, "fetch", key])
        if run([*command, "status", key]).returncode != 0:
            problems.append(f"the {key} data could not be downloaded (python3 -m planetgen.cli.datasets fetch {key}); "
                            "the tests that read it would give other stars")
    return problems


def preflight(environ, browser=True, say=print):
    """
    Makes a run's footing sound.

    Returns:
        list: The problems it could not repair, each a sentence (empty when the run may start). What it repaired it
        reports through `say`.
    """
    problems = []
    loaded = load_profile(environ)
    if loaded:
        say(f"preflight: loaded {', '.join(loaded)} from {PROFILE}")
    missing = [name for name in REQUIRED if not environ.get(name)]
    if missing:
        problems.append(f"{', '.join(missing)} not set (the tests would log in as the wrong account or skip); "
                        f"export them, they are written to {PROFILE} on a lane image")
        return problems
    state = mysql_state(environ)
    if state in ("down", "hung"):
        say(f"preflight: MySQL is {state}; starting it")
        if not start_mysql(environ):
            problems.append("MySQL does not answer and could not be started (mysqld_safe --user=mysql, or "
                            "PLANETGEN_TEST_MYSQL_START)")
        else:
            say("preflight: MySQL is up again")
            state = "up"
    elif state == "denied":
        problems.append("MySQL refused the PLANETGEN_TEST_MYSQL_USER / PASSWORD login; every database test would fail")
    if state == "up" and free_mb() < LOW_DISK_MB:
        say(f"preflight: only {free_mb():,} MB free on the disk; purging the MySQL binary logs")
        purge_binlogs(environ)
        if free_mb() < LOW_DISK_MB:
            problems.append(f"only {free_mb():,} MB free on the disk after purging the binlogs; free space first")
    if redis_state(environ) != "up":
        say("preflight: Redis is down; starting it")
        if not start_redis(environ):
            problems.append("Redis does not answer and could not be started (redis-server --daemonize yes, or "
                            "PLANETGEN_TEST_REDIS_START)")
        else:
            say("preflight: Redis is up again")
    problems += package_problems(browser)
    problems += dataset_problems(say)
    return problems


class Watchdog:
    """
    Asks MySQL and Redis every `interval` seconds while a phase runs. A server that is down (or hung `HUNG_AFTER`
    answers in a row) is restarted; `tripped` is set with a sentence in `what`, and `on_trip` (the runner's way to stop
    the phase's pytest) is called once. `events` keeps every restart for the run's closing report.
    """

    def __init__(self, environ, on_trip, interval=INTERVAL, mysql=mysql_state, redis=redis_state,
                 restart_mysql=start_mysql, restart_redis=start_redis, say=print):
        self.environ, self.on_trip, self.interval = environ, on_trip, interval
        self._mysql, self._redis, self._restart_mysql, self._restart_redis = mysql, redis, restart_mysql, restart_redis
        self.say = say
        self.events = []
        self.what = ""
        self.failed = False
        self.tripped = threading.Event()
        self._stop = threading.Event()
        self._thread = None
        self._hung = 0

    def start(self):
        self.tripped.clear()
        self._stop.clear()
        self._thread = threading.Thread(target=self._loop, name="test-service-watchdog", daemon=True)
        self._thread.start()

    def stop(self):
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=self.interval + 30)

    def check(self):
        """One round of questions; returns the sentence about what it found down and restarted (or could not), else ""."""
        found = ""
        state = self._mysql(self.environ)
        self._hung = self._hung + 1 if state == "hung" else 0
        if state == "down" or self._hung >= HUNG_AFTER:
            self._hung = 0
            found = self._repair("MySQL", self._restart_mysql)
        if self._redis(self.environ) == "down":
            found = (found + "; " if found else "") + self._repair("Redis", self._restart_redis)
        return found

    def _repair(self, name, restart):
        at = time.strftime("%H:%M:%SZ", time.gmtime())
        ok = restart(self.environ)
        sentence = f"{name} stopped answering and was {'restarted' if ok else 'NOT restarted'} at {at}"
        self.events.append(sentence)
        if not ok:
            self.failed = True
        return sentence

    def _loop(self):
        while not self._stop.wait(self.interval):
            found = self.check()
            if found:
                self.what = found
                self.tripped.set()
                self.on_trip(found)
                return
