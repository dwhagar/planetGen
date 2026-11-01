# tests/conftest.py

"""
Shared pytest fixtures for database-backed tests.

Every test that touches the database needs a real MySQL server to run
against (TODO.md Phase 5's MySQL port dropped SQLite entirely, so there's
no more in-process/file-based fallback for a database-backed test the way
a per-test `tmp_path` SQLite file used to be). The `mysql_config` fixture
below is the MySQL equivalent: it creates a uniquely-named, throwaway
database before the test and drops it afterward, so tests never see each
other's data and never need their own manual cleanup -- the same
guarantee a fresh `tmp_path` file gave for free, just implemented against
a shared server instead of the filesystem.

Point these tests at a real MySQL server with `PLANETGEN_TEST_MYSQL_HOST`/
`_PORT`/`_USER`/`_PASSWORD` (falling back to the same `PLANETGEN_MYSQL_*`
variables every other entry point reads if the `_TEST_` ones aren't set,
so a single already-configured local MySQL server needs no extra
environment just to run the test suite against it). If neither is
reachable, every test depending on `mysql_config` is skipped rather than
failing the whole run -- this project's CI is expected to provision a
real MySQL service for this reason, but a contributor without one locally
can still run everything else.
"""

import functools
import os
import re
import shutil
import sys
import tempfile
import urllib.error
import uuid

import pymysql
import pytest

from planetgen import tuning
from planetgen.db import store as _db
from planetgen.db.store import MySQLConfig

from tests import areas, browser_waits, dbpool, failure_capture, inline_queue, private_redis, publicids, stray_workers

# Every program a test starts (the generate CLI, a worker, a stand-in script) imports `planetgen` from this
# checkout, not from whichever checkout the editable install points at (TEST.143, TEST.145): from a git
# worktree the install is the main checkout, so a child process tested old code against new tests.
os.environ["PYTHONPATH"] = os.pathsep.join(
    [os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))), *filter(None, [os.environ.get("PYTHONPATH")])])

# The always-on activity log (planetgen/admin/activity_log.py) would
# otherwise go to /var/log/planetgen/ (or be skipped there); tests write it
# to a throwaway folder instead. Tests that read it point
# PLANETGEN_LOG_DIR at their own tmp_path and call activitylog.reset().
os.environ.setdefault("PLANETGEN_LOG_DIR", tempfile.mkdtemp(prefix="planetgen-test-logs-"))

# Maintenance mode's state file (UX.112) stays out of the real /var/lib/planetGen, and the autouse fixture
# `_site_open` below removes it before every test, since a test that migrates an old database turns it on.
os.environ["PLANETGEN_STATE_DIR"] = tempfile.mkdtemp(prefix="planetgen-test-state-")
os.environ["PLANETGEN_MAINTENANCE_FILE"] = os.path.join(tempfile.mkdtemp(prefix="planetgen-test-maintenance-"),
                                                        "maintenance.json")

# A run's log is written to the jobs directory when it ends (PERF.99); the tests' go to a throwaway folder too.
os.environ.setdefault("PLANETGEN_JOBS_DIR", tempfile.mkdtemp(prefix="planetgen-test-jobs-"))

# Generation runs one sector at a time, in-process, unless a test asks for
# workers itself (tests/test_work_queue.py) or the run sets
# PLANETGEN_WORKERS: CI's worker legs (TEST.74) run the `generation` tests
# at 2 and 4 workers too. A test that patches planetgen for a run does
# it with tests/worker_patches.py, which reaches the workers as well.
os.environ.setdefault("PLANETGEN_WORKERS", "1")

# A run with more than one worker queues its tasks on Redis (PERF.24):
# the test server's, when PLANETGEN_TEST_REDIS_URL names one (CI does), in a
# database of this xdist worker's own (TEST.266, TEST.280): a process a test
# starts takes the machine's one live-run lock and its queues, which another
# worker's process would otherwise meet. Without a Redis server such a run
# falls back to one worker.
if os.environ.get("PLANETGEN_TEST_REDIS_URL"):
    os.environ["PLANETGEN_REDIS_URL"] = private_redis.url_for_worker(os.environ["PLANETGEN_TEST_REDIS_URL"])

# Request limits and login lockouts count on that Redis server in a real
# install (SEC.30). Tests count in memory, per process, so parallel workers
# and neighbouring tests can't see each other's counts; the tests of the
# Redis counting take `redis_server` and set their own storage and prefix.
os.environ.setdefault("PLANETGEN_RATELIMIT_STORAGE_URI", "memory://")

# The admin Generate page's jobs run on Redis too (PERF.24 step 3), so a
# test that starts one for real takes `redis_server`: it skips without one.

# Test runs never read or add to the control database's generation speed
# and size stats (PERF.3, PERF.10); tests/test_generation_stats.py turns
# them on where it checks them.
os.environ.setdefault("PLANETGEN_GENERATION_STATS", "0")
# A data change never starts the warm-up by itself in a test (PERF.141); tests/test_warm_start.py turns it on.
os.environ.setdefault("PLANETGEN_WARM_ON_CHANGE", "0")
# PERF.209: the generator battery digest takes about half a minute and the code does not change during a run, so the workers share one copy
# (kept under a key of every file it reads, so a changed generator works it out again).
os.environ.setdefault("PLANETGEN_BATTERY_DIGEST_CACHE", os.path.join(tempfile.gettempdir(), "planetgen-battery-digest"))
# PERF.218: a first sector fill builds hundreds of bright-star sampling tables (minutes), the same ones in every process of a run and the next
# run; they are kept under a key of every file, dataset and tuning constant they read, so a changed generator builds them again.
os.environ.setdefault("PLANETGEN_BAND_TABLE_CACHE", os.path.join(tempfile.gettempdir(), "planetgen-band-tables"))
# Table pages count exactly in the tests; tests/test_count_cache.py turns the stored counts (PERF.64) on.
os.environ.setdefault("PLANETGEN_COUNT_CACHE", "off")

# Every pooled connection runs at MySQL 8's default sql_mode (TEST.7) --
# ONLY_FULL_GROUP_BY and STRICT_TRANS_TABLES included, even on a server
# whose own mode is looser. Set PLANETGEN_MYSQL_SQL_MODE yourself to try
# another mode.
os.environ.setdefault(_db.SQL_MODE_ENV_VAR, _db.STRICT_SQL_MODE)

# The control schema (admin logins, sessions, the work queue's lease) is
# one fixed name by default (`planetgen_control`), shared by everything on
# a server. Code that falls back to it (migrateDb seeding, the work queue)
# would make parallel pytest-xdist workers share one control database, and
# a developer's own PLANETGEN_CONTROL_DATABASE could point tests at a real
# one. Each test process gets its own throwaway name instead, dropped when
# the session ends (`_drop_session_control_database` below). It doesn't
# start with the `planetgen` prefix, so `/api/databases` never lists it.
_XDIST_WORKER = os.environ.get("PYTEST_XDIST_WORKER", "main")
SESSION_CONTROL_DATABASE = f"pgtest_control_{_XDIST_WORKER}_{uuid.uuid4().hex[:8]}"
os.environ[_db.CONTROL_DB_ENV_VAR] = SESSION_CONTROL_DATABASE


def pytest_configure(config):
    """Failure capture (TEST.139): with PLANETGEN_FAILURE_LOG set (scripts/test.py does), every failure's full report is kept."""
    failure_capture.register(config)
    browser_waits.install()  # a longer default wait for every browser context (TEST.126)


dbpool.watch_checkouts()   # TEST.193: a database that cannot be dropped names the code that holds its connection


def _test_server_kwargs():
    """Connection kwargs (host/port/user/password -- no database) for the
    MySQL server tests run against."""
    return dict(
        host=os.environ.get("PLANETGEN_TEST_MYSQL_HOST", os.environ.get("PLANETGEN_MYSQL_HOST", "127.0.0.1")),
        port=int(os.environ.get("PLANETGEN_TEST_MYSQL_PORT", os.environ.get("PLANETGEN_MYSQL_PORT", 3306))),
        user=os.environ.get("PLANETGEN_TEST_MYSQL_USER", os.environ.get("PLANETGEN_MYSQL_USER", "planetgen")),
        password=os.environ.get("PLANETGEN_TEST_MYSQL_PASSWORD", os.environ.get("PLANETGEN_MYSQL_PASSWORD", "")),
    )


def pytest_sessionstart(session):
    """
    The math check gate (TEST.67): before any test runs, the whole
    `planetgen.physics.mathcheck` list runs once (well under a second), and a
    failure stops the run there -- with the math broken, every other
    failure would be noise. Under pytest-xdist only the controller runs it;
    the workers start after it passed. `test_math_check.py` repeats each
    check as its own test, for the report.
    """
    if hasattr(session.config, "workerinput"):
        return
    from planetgen.physics import mathcheck as mathCheck

    results = mathCheck.run_all()
    if mathCheck.failures(results):
        pytest.exit(
            "The math check failed, so the rest of the suite would only be noise. Fix these first "
            "(python -m planetgen.physics.mathcheck -v):\n" + mathCheck.format_report(results),
            returncode=1,
        )


def _generation_refs():
    """What a test file that runs generation has imported."""
    refs = [sys.modules.get("generate"), sys.modules.get("planetgen.queue.work")]
    support = sys.modules.get("tests.bughunt_support")
    if support is not None:
        refs.append(support.run_cli)
    return [ref for ref in refs if ref is not None]


PLAYWRIGHT_IMPORT = re.compile(r"""importorskip\(\s*["']playwright|^\s*(?:from|import)\s+playwright\b""", re.M)
"""re.Pattern: A test file's own import of Playwright."""

BROWSER_FIXTURES = frozenset({"browser", "gl_browser"})
"""frozenset: The Playwright browser fixtures (`test_web_a11y.py`, `test_web_browser_fixture_maps.py`), for a file
that borrows one without importing Playwright itself."""


@functools.lru_cache(maxsize=None)
def _imports_playwright(path):
    with open(path, encoding="utf-8") as handle:
        return bool(PLAYWRIGHT_IMPORT.search(handle.read()))


def _uses_playwright(item):
    """
    Whether `item` drives headless Chromium (TEST.304): its file imports Playwright or it takes a browser fixture.
    So the all-core phase of `scripts/test.py` never starts Chromium, whatever the file is called.
    """
    return bool(BROWSER_FIXTURES & set(getattr(item, "fixturenames", ()))) or _imports_playwright(str(item.path))


def pytest_collection_modifyitems(config, items):
    """Suite markers (registered in pytest.ini), so a run can pick a slice:
    `-m "not db"` without a database, `-m "not slow"` for a quick loop,
    `-m browser` for the Playwright checks alone, `-m mathcheck` for the
    math check gate.

    - `db`: the test uses a real MySQL database (`mysql_config`, directly or
      through another fixture).
    - `slow`: the brute-force and seeded-sweep files (`test_fuzz_*`,
      `test_bughunt_*`), whose cost is in their number of examples.
    - `browser`: every test that drives headless Chromium (`_uses_playwright`, TEST.304).
    - `mathcheck`: `test_math_check.py`, moved to the front of the run.
    - `generation`: files that run generation (`generate`, the work queue
      or `bughunt_support.run_cli`), which CI also runs at several worker
      counts (TEST.74).
    - `globaldocs`: a test that compares a generated global document with its source; skipped unless
      `PLANETGEN_GLOBAL_DOCS=1` (DOC.69).
    - `heavy`: the files that took over 100 s each in the 2026-10-14
      baseline (`tests/areas.py`); `scripts/test.py quick` leaves them out.
    - `area_<name>`: the part of planetGen a file covers (`tests/areas.py`),
      by file name, so `-m area_web` picks the web tests without moving a file.
    """
    generation_refs = _generation_refs()
    for item in items:
        fixturenames = getattr(item, "fixturenames", ())
        if "mysql_config" in fixturenames or "_mysql_server_available" in fixturenames:
            item.add_marker(pytest.mark.db)
        if item.get_closest_marker("globaldocs") and not os.environ.get("PLANETGEN_GLOBAL_DOCS"):
            item.add_marker(pytest.mark.skip(reason="a global document's freshness is the TODO thread's: "
                                                    "PLANETGEN_GLOBAL_DOCS=1 or scripts/test.py --global-docs (DOC.69)"))
        if item.get_closest_marker("budget") and not os.environ.get("PLANETGEN_BUDGET"):
            item.add_marker(pytest.mark.skip(reason="the big-galaxy query budget is opt-in: scripts/test.py budget (TEST.133)"))
        name = item.path.name
        item.add_marker(getattr(pytest.mark, "area_" + areas.area_of(name)))
        if areas.is_heavy(name):
            item.add_marker(pytest.mark.heavy)
        if name.startswith(("test_fuzz_", "test_bughunt_")):
            item.add_marker(pytest.mark.slow)
        if _uses_playwright(item):
            item.add_marker(pytest.mark.browser)
        module = getattr(item, "module", None)
        if module is not None and any(value is ref for value in vars(module).values() for ref in generation_refs):
            item.add_marker(pytest.mark.generation)
    # The math check's own tests run first (a stable sort keeps every
    # other test in its collected order).
    items.sort(key=lambda item: item.get_closest_marker("mathcheck") is None)


@pytest.fixture(scope="session")
def _mysql_server_available():
    """
    Session-scoped probe: attempts one connection to the configured test
    server and skips every test that (transitively, via `mysql_config`)
    depends on this fixture if it can't be reached. Pytest caches a
    session-scoped fixture's outcome -- including a skip -- for the whole
    run, so this only actually dials the server once regardless of how
    many tests need it.
    """
    try:
        conn = pymysql.connect(**_test_server_kwargs(), connect_timeout=3)
        conn.close()
    except pymysql.MySQLError as exc:
        pytest.skip(
            f"No MySQL test server reachable ({exc}) -- set PLANETGEN_TEST_MYSQL_HOST "
            f"(and _PORT/_USER/_PASSWORD as needed) to run database-backed tests."
        )
    yield
    # The per-process control database (SESSION_CONTROL_DATABASE above),
    # if anything created it.
    conn = pymysql.connect(**_test_server_kwargs(), connect_timeout=3)
    try:
        with conn.cursor() as cur:
            cur.execute(f"DROP DATABASE IF EXISTS `{SESSION_CONTROL_DATABASE}`")
        conn.commit()
    finally:
        conn.close()


@pytest.fixture(scope="session")
def _worker_database(_mysql_server_available):
    """The database this test process reuses (tests/dbpool.py, TEST.139); dropped when the process ends."""
    worker = dbpool.WorkerDatabase(_test_server_kwargs(), deep_check=bool(os.environ.get("PLANETGEN_TEST_DEEP_CHECK")))
    yield worker
    worker.close()


@pytest.fixture
def mysql_config(request, _mysql_server_available):
    """
    A MySQL database for one test: empty, with the schema in place.

    A test process keeps one such database and puts it back the way a new one is after each test
    (`tests/dbpool.py`), which costs a few milliseconds where creating, filling with the schema and dropping
    a database cost 0.66 s. A test that needs a database with no schema at all, or one that is not shared
    with the tests before it, is marked `fresh_database` and gets `fresh_mysql_config` below;
    `PLANETGEN_TEST_FRESH_DB=1` does that for every test.

    Yields:
        MySQLConfig: Points at the empty database -- pass straight to `get_connection`/`save_sector`/
            `save_system`/etc.
    """
    if os.environ.get("PLANETGEN_TEST_FRESH_DB") or request.node.get_closest_marker("fresh_database"):
        yield request.getfixturevalue("fresh_mysql_config")
        return
    worker = request.getfixturevalue("_worker_database")
    config = worker.lease()
    publicids.use(config)
    try:
        yield config
    finally:
        publicids.use(None)
        worker.release()


@pytest.fixture
def fresh_mysql_config(_mysql_server_available):
    """
    Creates a uniquely-named, empty MySQL database for the duration of one test -- no tables at all, until
    something connects and applies the schema -- and drops it afterward regardless of the test's outcome.

    Every test gets its own database name, so `planetgen.db.store`'s
    module-level pool cache (keyed by connection params, database
    included) would otherwise grow by one pool -- and its
    `mincached` real connections -- per test for the life of the process.
    Closing that pool here, once this test's database is being dropped for
    good, keeps a long test run from exhausting the server's
    `max_connections`.

    Yields:
        MySQLConfig: Points at the fresh, empty database.
    """
    kwargs = _test_server_kwargs()
    db_name = f"planetgen_test_{uuid.uuid4().hex[:16]}"

    admin_conn = pymysql.connect(**kwargs)
    try:
        with admin_conn.cursor() as cur:
            cur.execute(f"CREATE DATABASE `{db_name}`")
        admin_conn.commit()
    finally:
        admin_conn.close()

    config = MySQLConfig(database=db_name, **kwargs)
    publicids.use(config)
    try:
        yield config
    finally:
        publicids.use(None)
        _db.close_pool(config)
        dbpool.drop_database(kwargs, db_name)


# Registers and loads the hypothesis profiles for every `test_fuzz_*.py`
# file (see tests/fuzz_support.py) before any of them is collected.
from tests import fuzz_support  # noqa: E402,F401


FAST_PASSWORD_HASH_PARAMETERS = (8, 1, 1)
"""Argon2id's smallest legal cost: 8 KiB, one pass, one thread."""


@pytest.fixture(autouse=True)
def _fast_password_hashing(request, monkeypatch):
    """Every login, credential change and seeded admin user runs Argon2id at
    19 MiB and two passes in production (`adminAuth.hash_parameters`, SEC.32),
    tens of milliseconds per check -- thousands of fuzzed logins spend
    minutes on it without testing anything the cost changes. Tests use the
    smallest cost instead (same algorithm, same code path, same
    rehash-on-login logic), except those marked `real_password_hashing`
    (test_password_hardening.py), which check the production setting
    itself."""
    if request.node.get_closest_marker("real_password_hashing"):
        yield
        return
    from planetgen.admin import auth as adminAuth
    monkeypatch.setattr(adminAuth, "hash_parameters", lambda: FAST_PASSWORD_HASH_PARAMETERS)
    # The unknown-username dummy hash is cached; make the next one cheap
    # too (monkeypatch puts the production one back afterwards).
    monkeypatch.setattr(adminAuth, "_dummy_password_hash", None)
    yield


@pytest.fixture(autouse=True)
def _toy_galaxies_scatter_inside_their_outline(monkeypatch):
    """The star passes cover the grid's whole reach, 262 kpc by default (GEN.322), the bands beyond the outline
    included. A test galaxy is a few dozen parsecs across with an outline written by hand, and its fixed-size halo
    would fill the rest of the grid with stars, so a test draws inside the outline unless it asks for the reach
    (`tuning.SCATTER_COVERS_REACH`, set True by the tests of GEN.322)."""
    monkeypatch.setattr(tuning, "SCATTER_COVERS_REACH", False)
    yield


@pytest.fixture
def unit_density_plans(monkeypatch):
    """A test that plans a toy galaxy (a disk of a few dozen parsecs) from the command line gets it in relative units,
    density 1 at its reference radius (`tests/toy_shape.py`): the plan's own shape would hold the stellar mass its seed
    draws in that tiny volume (GEN.404)."""
    from planetgen.galaxy import shape_draw
    from tests import toy_shape
    monkeypatch.setattr(shape_draw, "plan_shape", toy_shape.plan_at_unit_density(shape_draw.plan_shape))
    yield


@pytest.fixture(autouse=True)
def _api_jobs_queue_of_its_own(monkeypatch):
    """The API jobs share one queue and its workers in production (PERF.39); a test, running beside others
    on one Redis server, gets a queue of its own, so its job is never picked up by another test's worker
    (which has that test's environment)."""
    import uuid
    from planetgen.queue import api_jobs
    monkeypatch.setattr(api_jobs, "QUEUE_NAME", f"planetgen-api-test-{uuid.uuid4().hex[:12]}")
    monkeypatch.setattr(api_jobs, "LONG_QUEUE_NAME", f"planetgen-api-test-long-{uuid.uuid4().hex[:12]}")
    yield


@pytest.fixture(autouse=True)
def _wanted_workers_of_its_own(monkeypatch):
    """The admin's wanted worker count is one key on the Redis server, which every test process reads when it counts
    the workers a queue should have; a test that sets it (the admin's page, `set_wanted`) beside one that counts them
    changed the other's answer (TEST.255: 3 where the test had made it 2). Each test gets a key of its own, emptied
    after."""
    import uuid
    from planetgen.queue import workers
    key = f"{workers.WANTED_KEY}:test-{uuid.uuid4().hex[:12]}"
    monkeypatch.setattr(workers, "WANTED_KEY", key)
    yield
    url = os.environ.get("PLANETGEN_REDIS_URL")
    if url:
        try:
            import redis
            redis.Redis.from_url(url).delete(key)
        except Exception:  # noqa: BLE001 -- no server: nothing was set
            pass


@pytest.fixture(autouse=True)
def _api_jobs_wait_patiently(monkeypatch):
    """A queued API edit (`run_queued`) starts a burst worker process and waits `SHORT_WAIT_SECONDS` for
    it; in a busy parallel run the worker can take longer than that to start, and the route answers
    202 instead of 200 (TEST.111, TEST.112). The tests that want the 202 patch `api_jobs.wait`
    themselves, so here the wait is long enough that a slow start is still a 200."""
    from planetgen.queue import api_jobs
    real_wait = api_jobs.wait
    monkeypatch.setattr(api_jobs, "wait", lambda job_id, seconds=120.0: real_wait(job_id, seconds))
    yield


def _counts_in_this_process(limiter):
    """Whether the limiter keeps its counts in this process's memory. Its `reset()` on a Redis storage deletes
    every limiter key on the server, including those of a test running at that moment in another process (a
    test of `test_ratelimit_redis.py` lost half its counts that way, TEST.146); a test that uses Redis keeps
    its keys under a prefix of its own and removes them itself (`key_prefix`)."""
    from limits.storage import MemoryStorage
    storage = getattr(limiter, "_storage", None)
    return storage is None or isinstance(storage, MemoryStorage)


@pytest.fixture(autouse=True)
def _reset_rate_limits():
    """Flask-Limiter's one shared instance keeps its counts for the whole
    process, so a test would otherwise start with whatever the tests
    before it (and so the order xdist deals them in) left: its counts and
    whether limits are on at all (TEST.83)."""
    def reset():
        try:
            from planetgen.api.limiter import limiter
            # `init_app` keeps the last app's RATELIMIT_ENABLED on the shared
            # instance, so a test of an app with limits off (RATELIMIT_ENABLED =
            # False) switched them off for every later test in its worker.
            limiter.enabled = True
            if getattr(limiter, "initialized", False) and _counts_in_this_process(limiter):
                limiter.reset()
        except Exception:  # noqa: BLE001 -- no API dependencies, or a storage that can't reset
            pass
    reset()
    yield
    reset()


@pytest.fixture(autouse=True)
def _reset_login_backoff():
    """The login lockouts' in-memory fallback (`api/loginguard.py`) is one
    process-wide store, so one test's failed logins (against an app with
    no real control database) would otherwise lock the same username or
    address in every later test. The database-backed counts live in each
    test's own throwaway database."""
    try:
        from planetgen.api.loginguard import memory_store
    except ImportError:  # the API's dependencies aren't installed
        yield
        return
    memory_store.clear()
    yield
    memory_store.clear()


_DATASETS_ARE_REAL = False
"""bool: Whether the running test is marked `real_datasets` (set by `_datasets_are_real`)."""


@pytest.fixture(scope="session", autouse=True)
def _downloaded_data_present():
    """
    Bulk generation refuses to start while a downloaded data file is missing (GEN.590, `run_galaxy.require_datasets`), and a
    test machine's cache folder is empty (a server's is filled by `update.sh`), so for the whole session the check reports
    nothing missing, including for the module- and session-scoped fixtures that generate a galaxy and for a test that calls
    `monkeypatch.undo()`. A test marked `real_datasets` gets the real check.
    """
    from unittest import mock

    from planetgen.datasets import catalogue

    real = catalogue.missing
    with mock.patch.object(catalogue, "missing", lambda: real() if _DATASETS_ARE_REAL else []):
        yield


@pytest.fixture(autouse=True)
def _datasets_are_real(request):
    global _DATASETS_ARE_REAL
    _DATASETS_ARE_REAL = request.node.get_closest_marker("real_datasets") is not None
    yield
    _DATASETS_ARE_REAL = False


@pytest.fixture(scope="session", autouse=True)
def _worker_redis_database():
    """The worker's own Redis database (see the top of this file) is empty when its tests start and when they end, and no
    queue worker process is left running on it (TEST.393)."""
    url = os.environ.get("PLANETGEN_REDIS_URL")
    if url:
        stray_workers.reap_database(url)    # TEST.393: a worker an earlier, stopped session left on this database
        private_redis.flush(url)
    yield
    if url:
        stray_workers.reap_database(url)    # and the ones this session started and did not stop
        private_redis.flush(url)


@pytest.fixture(autouse=True)
def _run_queue_inline(monkeypatch, request):
    """
    A stage on the new work queue (PERF.97, `planetgen/queue/units.py`) queues its units on Redis and has workers
    run them. A test gets a queue and a key prefix of its own, so parallel tests (pytest-xdist) neither share the
    machine's one active run nor take each other's units, and its `ensure_workers` runs the waiting units right
    here, in this process, so a patch the test made reaches them and no worker process is left running when the
    test's database goes. The workers read the MySQL server's credentials from the settings, so the environment
    carries the test server's. Nothing is patched, and nothing is touched on Redis, without a Redis server.
    A test that starts real workers or reads the shared queues' names is marked `real_queue` and gets none of this, and
    neither does a page test that visits the sample run (`test_web_a11y.sample_run`): the run is written once for the
    module, under the key prefix in force when the module's fixtures are made, and the page must find it under the same.
    """
    if request.node.get_closest_marker("real_queue") or "sample_run" in request.fixturenames:
        yield
        return
    with inline_queue.inline_queue(monkeypatch):
        yield


ROOMY_DISK_BYTES = 200 * 1024 ** 3
"""int: The free space (and twice the size) a test sees on the database's drive."""


@pytest.fixture(autouse=True)
def _roomy_database_disk(monkeypatch):
    """
    A generate run refuses to start when it would leave under 5 GB free on the database's drive (PERF.3), and `--strict`
    makes that an exit with status 1. A test machine whose drive is nearly full (the MySQL binary logs of a long
    session fill one) then failed every test that generates under `--strict` or reads the estimate's `refused`
    (TEST.268). The code that reaches `stats.database_disk` through the module sees a drive of at least 200 GB free;
    a test that wants another gives its own (`monkeypatch.setattr(generationStats, "database_disk", ...)`), and
    `test_database_disk.py` imports the function itself and measures the real thing.
    """
    import dataclasses

    from planetgen.generation import stats

    real = stats.database_disk

    def roomy(galaxy_conn, mysql_host, database=None):
        found = real(galaxy_conn, mysql_host, database)
        if isinstance(found, stats.DiskSpace):
            return dataclasses.replace(found, total_bytes=max(found.total_bytes, 2 * ROOMY_DISK_BYTES),
                                       free_bytes=max(found.free_bytes, ROOMY_DISK_BYTES))
        return found
    monkeypatch.setattr(stats, "database_disk", roomy)


@pytest.fixture
def redis_server():
    """The Redis URL jobs and workers use; skips the test without one."""
    url = os.environ.get("PLANETGEN_REDIS_URL")
    if not url:
        pytest.skip("no Redis server (PLANETGEN_TEST_REDIS_URL is not set)")
    from planetgen.queue import redisqueue
    try:
        redisqueue.connect(url).ping()
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"no Redis server at {url}: {exc}")
    return url


@pytest.fixture
def key_prefix(redis_server):
    """A key prefix of this test's own on the test server, emptied after."""
    import uuid

    import redis

    prefix = f"planetgen-test:{uuid.uuid4().hex}"
    yield prefix
    client = redis.Redis.from_url(redis_server)
    for key in client.scan_iter(match=f"{prefix}:*"):
        client.delete(key)


@pytest.fixture(autouse=True)
def _site_open():
    """Every test starts with maintenance mode off (UX.112) and no update running (OPS.41)."""
    try:
        os.unlink(os.environ["PLANETGEN_MAINTENANCE_FILE"])
    except FileNotFoundError:
        pass
    from planetgen.admin import updateguard
    try:
        os.unlink(updateguard.lock_path())
    except FileNotFoundError:
        pass
    shutil.rmtree(updateguard.running_dir(), ignore_errors=True)


@pytest.fixture
def large_spiral_fallback(monkeypatch):
    """The galaxy a unit with no galaxy bound belongs to is a 2.6 kpc spiral, as the hand-built shapes of the tests that
    place things at 8 kpc assume: the fallback seed's own galaxy is a small spiral since the Auto window narrowed (GEN.603)."""
    from planetgen.galaxy import scales
    monkeypatch.setattr(scales, "_drawn_scale_length_pc", lambda _seed: 2600.0)
