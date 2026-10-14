# tests/test_budget_support.py

"""
The big-galaxy budget tier (TEST.133, `test_budget_queries.py`) is opt-in and builds millions of systems, so
nothing in the ordinary run would notice it stop working. This builds a galaxy of a few thousand systems the same
way and runs every case of the tier on it, so a change to the schema or to a query function that breaks the tier
shows up here and not the next time somebody builds the big one.
"""

import uuid

import pymysql
import pytest

from planetgen.db import query
from planetgen.db.store import MySQLConfig

from tests import budget_support, test_budget_queries
from tests.conftest import _test_server_kwargs

pytestmark = pytest.mark.db

SYSTEMS, SECTORS = 3_000, 300


@pytest.fixture(scope="module")
def small_galaxy():
    server = _test_server_kwargs()
    name = f"planetgen_test_budget_{uuid.uuid4().hex[:12]}"
    budget_support.build(server, name, SYSTEMS, SECTORS, log=lambda *_: None, drop_others=False)
    yield MySQLConfig(database=name, **server)
    conn = pymysql.connect(**server)
    try:
        with conn.cursor() as cur:
            cur.execute(f"DROP DATABASE IF EXISTS `{name}`")
    finally:
        conn.close()


def test_the_budget_galaxy_holds_the_systems_it_was_asked_for(small_galaxy):
    conn = query.open_readonly(small_galaxy, statement_timeout_s=30)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE id > 1000000000").fetchone()["n"] == SYSTEMS - 1
        assert conn.execute("SELECT COUNT(*) AS n FROM stars").fetchone()["n"] == SYSTEMS
        assert conn.execute("SELECT COUNT(DISTINCT uid) AS n FROM star_systems").fetchone()["n"] == SYSTEMS
    finally:
        conn.close()


@pytest.mark.parametrize("name, call", [(name, call) for name, call, _, _ in test_budget_queries._cases()],
                         ids=[name for name, *_ in test_budget_queries._cases()])
def test_every_budget_case_runs_and_counts_the_rows_it_read(small_galaxy, name, call):
    conn = query.open_readonly(small_galaxy, statement_timeout_s=30)
    try:
        read, seconds = budget_support.rows_read(conn, call)
    finally:
        conn.close()
    assert read > 0
