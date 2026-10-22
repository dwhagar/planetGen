# tests/test_migration_request_support.py

"""
DB.39: migration 0126 makes `request_comments` and `request_notifications` and touches nothing else. A migration run again leaves
the tables and their rows as they are, a comment has a known author role, and both tables go with their plan.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

TABLES = ("request_comments", "request_notifications")


def _tables(conn):
    return {row["t"] for row in conn.execute(
        "SELECT TABLE_NAME AS t FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()").fetchall()}


def _plan(conn):
    return conn.execute("INSERT INTO change_plans (requester_name, requester_role, state) VALUES ('Boss', 'user', 'requested')").lastrowid


def _comment(conn, plan, role="user", body="Please"):
    return conn.execute("INSERT INTO request_comments (plan_id, author_name, author_role, body) VALUES (%s, 'Boss', %s, %s)",
                        (plan, role, body)).lastrowid


def test_a_v125_galaxy_gets_empty_tables(mysql_config):
    load_old_schema(mysql_config, 125)
    assert store.schema_status(mysql_config)[0] == 125
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert set(TABLES) <= _tables(conn)
        for table in TABLES:
            assert conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_tables_and_their_rows(mysql_config):
    load_old_schema(mysql_config, 125)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _comment(conn, _plan(conn), body="Make it an M planet")
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0125')")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT body FROM request_comments").fetchall() == [{"body": "Make it an M planet"}]
    finally:
        conn.close()


def test_a_comment_has_a_known_author_role_and_both_tables_go_with_their_plan(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        plan = _plan(conn)
        for role in ("user", "moderator", "admin"):
            _comment(conn, plan, role=role)
        with pytest.raises(Exception, match="chk_request_comments_author_role"):
            _comment(conn, plan, role="owner")
        conn.execute("INSERT INTO request_notifications (plan_id, recipient_id, kind) VALUES (%s, 7, 'comment')", (plan,))
        row = conn.execute("SELECT sent_at FROM request_notifications").fetchone()
        assert row == {"sent_at": None}
        conn.execute("DELETE FROM change_plans WHERE id = %s", (plan,))
        for table in TABLES:
            assert conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"] == 0
    finally:
        conn.close()
