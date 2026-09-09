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
import uuid

import pymysql
import pytest

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

    try:
        yield MySQLConfig(database=db_name, **kwargs)
    finally:
        admin_conn = pymysql.connect(**kwargs)
        try:
            with admin_conn.cursor() as cur:
                cur.execute(f"DROP DATABASE IF EXISTS `{db_name}`")
            admin_conn.commit()
        finally:
            admin_conn.close()
