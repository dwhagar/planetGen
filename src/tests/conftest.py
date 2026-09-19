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

import os
import urllib.error
import urllib.request
import uuid

import pymysql
import pytest

from stellarObjects import _db
from stellarObjects._db import MySQLConfig


def _test_server_kwargs():
    """Connection kwargs (host/port/user/password -- no database) for the
    MySQL server tests run against."""
    return dict(
        host=os.environ.get("PLANETGEN_TEST_MYSQL_HOST", os.environ.get("PLANETGEN_MYSQL_HOST", "127.0.0.1")),
        port=int(os.environ.get("PLANETGEN_TEST_MYSQL_PORT", os.environ.get("PLANETGEN_MYSQL_PORT", 3306))),
        user=os.environ.get("PLANETGEN_TEST_MYSQL_USER", os.environ.get("PLANETGEN_MYSQL_USER", "planetgen")),
        password=os.environ.get("PLANETGEN_TEST_MYSQL_PASSWORD", os.environ.get("PLANETGEN_MYSQL_PASSWORD", "")),
    )


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


@pytest.fixture
def mysql_config(_mysql_server_available):
    """
    Creates a uniquely-named, empty MySQL database for the duration of
    one test, and drops it afterward regardless of the test's outcome.

    Every test gets its own database name, so `stellarObjects._db`'s
    module-level pool cache (keyed by connection params, database
    included) would otherwise grow by one `PooledDB` -- and its
    `mincached` real connections -- per test for the life of the process.
    Closing that pool here, once this test's database is being dropped for
    good, keeps a long test run from exhausting the server's
    `max_connections`.

    Yields:
        MySQLConfig: Points at the fresh, empty database -- pass straight
            to `get_connection`/`save_sector`/`save_system`/etc.
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
    try:
        yield config
    finally:
        _db.close_pool(config)
        admin_conn = pymysql.connect(**kwargs)
        try:
            with admin_conn.cursor() as cur:
                cur.execute(f"DROP DATABASE IF EXISTS `{db_name}`")
            admin_conn.commit()
        finally:
            admin_conn.close()


@pytest.fixture(scope="session")
def wikijs_config():
    """
    Connection details for a real, disposable Wiki.js instance to test
    `wikiClient.wikijs.WikiJsBackend` against end-to-end -- unlike
    `mysql_config` above, this project has no service of its own to
    provision one, so this is opt-in only: set `PLANETGEN_TEST_WIKIJS_BASE_URL`
    and `PLANETGEN_TEST_WIKIJS_TOKEN` (a Personal API Token from that
    instance's Admin -> API Access) to run `test_wikiclient_wikijs_integration.py`
    against it; every test depending on this fixture is skipped, not
    failed, when either is unset or the instance isn't reachable -- same
    "opt-in real service, skip without it" treatment `_mysql_server_available`
    gives MySQL.

    Yields:
        tuple[str, str]: `(base_url, api_token)`.
    """
    base_url = os.environ.get("PLANETGEN_TEST_WIKIJS_BASE_URL")
    api_token = os.environ.get("PLANETGEN_TEST_WIKIJS_TOKEN")
    if not base_url or not api_token:
        pytest.skip(
            "No Wiki.js test instance configured -- set PLANETGEN_TEST_WIKIJS_BASE_URL "
            "and PLANETGEN_TEST_WIKIJS_TOKEN to run wikiClient wikijs integration tests."
        )

    try:
        # A plain GET against the instance's own root, not /graphql --
        # only checking that *something* answers at this host/port at all,
        # kept separate from WikiJsBackend itself so a real auth/GraphQL
        # failure (a wrong token, wrong Wiki.js version) surfaces as a
        # genuine test failure later rather than being swallowed here as
        # "not reachable". An HTTPError (any status) still means the host
        # answered, so only a connection-level URLError counts as
        # unreachable.
        urllib.request.urlopen(base_url, timeout=5).close()
    except urllib.error.HTTPError:
        pass
    except urllib.error.URLError as exc:
        pytest.skip(f"Wiki.js test instance at {base_url} not reachable: {exc}")

    yield base_url, api_token


@pytest.fixture(scope="session")
def mediawiki_config():
    """
    Connection details for a real, disposable MediaWiki instance to test
    `wikiClient.mediawiki.MediaWikiBackend` against end-to-end -- the
    MediaWiki counterpart of `wikijs_config` above, same opt-in treatment.
    Set `PLANETGEN_TEST_MEDIAWIKI_BASE_URL`, `PLANETGEN_TEST_MEDIAWIKI_USERNAME`
    (a `Special:BotPasswords` username, `"User@BotName"` form), and
    `PLANETGEN_TEST_MEDIAWIKI_PASSWORD` to run
    `test_wikiclient_mediawiki_integration.py` against it; every test
    depending on this fixture is skipped, not failed, when any of the three
    is unset or the instance isn't reachable.

    Yields:
        tuple[str, str, str]: `(base_url, username, password)`.
    """
    base_url = os.environ.get("PLANETGEN_TEST_MEDIAWIKI_BASE_URL")
    username = os.environ.get("PLANETGEN_TEST_MEDIAWIKI_USERNAME")
    password = os.environ.get("PLANETGEN_TEST_MEDIAWIKI_PASSWORD")
    if not base_url or not username or not password:
        pytest.skip(
            "No MediaWiki test instance configured -- set PLANETGEN_TEST_MEDIAWIKI_BASE_URL, "
            "PLANETGEN_TEST_MEDIAWIKI_USERNAME, and PLANETGEN_TEST_MEDIAWIKI_PASSWORD to run "
            "wikiClient mediawiki integration tests."
        )

    try:
        # Same "something answers at all" pre-check wikijs_config does,
        # against the instance root rather than api.php -- a real
        # login/edit failure should surface from MediaWikiBackend itself,
        # not be swallowed here as "not reachable".
        urllib.request.urlopen(base_url, timeout=5).close()
    except urllib.error.HTTPError:
        pass
    except urllib.error.URLError as exc:
        pytest.skip(f"MediaWiki test instance at {base_url} not reachable: {exc}")

    yield base_url, username, password
