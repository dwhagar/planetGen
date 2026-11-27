# tests/test_abandoned_test_databases.py

"""
TEST.852: a test run leaves no database behind, and the databases of a run that was killed are dropped by the next one.

A gate kills the processes that pass its time limit, and the kernel kills the ones that take too much memory; neither reaches its teardown, and every pytest process made a control database
(`pgtest_control_*`) that only the tests that used the MySQL server dropped. One lane box held 1,216 of them (8.2 GB) and ran out of disk. Now every process drops its control database whether
or not its tests used the server, the names of the test databases say which machine and process made them (`dbpool.OWNER`), and the next run to start drops the ones of a process that is gone
(`dbpool.drop_abandoned`). Skipped without the MySQL test server.
"""

import subprocess
import sys

import pymysql
import pytest

from tests import dbpool
from tests.conftest import _test_server_kwargs


def _databases(like):
    conn = pymysql.connect(**_test_server_kwargs())
    try:
        with conn.cursor() as cur:
            cur.execute("SHOW DATABASES LIKE %s", (like,))
            return [row[0] for row in cur.fetchall()]
    finally:
        conn.close()


def _create(name):
    conn = pymysql.connect(**_test_server_kwargs())
    try:
        with conn.cursor() as cur:
            cur.execute(f"CREATE DATABASE `{name}`")
            cur.execute(f"CREATE TABLE `{name}`.t (a INT)")
        conn.commit()
    finally:
        conn.close()


def test_a_run_of_tests_that_use_no_database_leaves_none_behind(_mysql_server_available):
    """A real pytest process on one test that needs no server: its control database is made at its start and must be gone at its end."""
    child = subprocess.Popen([sys.executable, "-m", "pytest", "tests/test_version_key.py::test_the_design_example", "-n", "0", "-q", "-p", "no:randomly", "-p", "no:cacheprovider"],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out, _ = child.communicate(timeout=600)
    assert child.returncode == 0, out[-2000:]
    left = _databases(f"%h%p{child.pid}\\_%")
    assert left == [], f"the run (process {child.pid}) left {left}"


def test_the_databases_of_a_process_that_is_gone_are_dropped_and_no_other_is(_mysql_server_available):
    dead = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(600)"])
    dead.kill()
    dead.wait()
    gone = dbpool.owned_name("pgtest_control_gw0_").replace(dbpool.OWNER, f"{dbpool.OWNER[:6]}{dead.pid}")
    mine = dbpool.owned_name("pgtest_control_gw0_")
    elsewhere = dbpool.owned_name("pgtest_control_gw0_").replace(dbpool.OWNER, f"h0000p{dead.pid}")      # another machine's: not known to be gone
    names = (gone, mine, elsewhere)
    for name in names:
        _create(name)
    try:
        dropped = dbpool.drop_abandoned(_test_server_kwargs(), by_age=False, only=names)
        assert gone in dropped and mine not in dropped and elsewhere not in dropped
        remaining = [name for name in names if _databases(name)]
        assert remaining == [mine, elsewhere]
    finally:
        for name in names:
            dbpool.drop_database(_test_server_kwargs(), name)


def test_a_database_with_no_owner_in_its_name_is_dropped_when_its_tables_are_old_enough(_mysql_server_available, monkeypatch):
    name = "planetgen_a11y_ownerless_for_test_852"
    _create(name)
    try:
        assert name not in dbpool.drop_abandoned(_test_server_kwargs(), by_age=True, only=(name,))       # made just now: kept
        monkeypatch.setattr(dbpool, "ABANDONED_HOURS", -1.0)       # a table made a moment ago is older than minus one hour
        assert name in dbpool.drop_abandoned(_test_server_kwargs(), by_age=True, only=(name,))
        assert _databases(name) == []
    finally:
        dbpool.drop_database(_test_server_kwargs(), name)
