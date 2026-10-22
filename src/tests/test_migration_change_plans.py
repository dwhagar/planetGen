# tests/test_migration_change_plans.py

"""
DB.36: migration 0122 makes `change_plans`, `change_plan_ops` and `edit_journal` and touches nothing else. A migration run again
leaves the tables and their rows as they are, a plan's operations go with the plan, and a plan that has journal rows cannot be
deleted.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

TABLES = ("change_plans", "change_plan_ops", "edit_journal")
UID = b"\x80" * 12


def _tables(conn):
    return {row["t"] for row in conn.execute(
        "SELECT TABLE_NAME AS t FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()").fetchall()}


def _plan(conn, state="draft", role="admin", **more):
    columns = {"requester_name": "Boss", "requester_role": role, "state": state, **more}
    return conn.execute(f"INSERT INTO change_plans ({', '.join(columns)}) VALUES ({', '.join(['%s'] * len(columns))})",
                        tuple(columns.values())).lastrowid


def test_a_v121_galaxy_gets_empty_tables_and_keeps_its_rows(mysql_config):
    load_old_schema(mysql_config, 121)
    assert store.schema_status(mysql_config)[0] == 121
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not set(TABLES) & _tables(conn)
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert set(TABLES) <= _tables(conn)
        for table in TABLES:
            assert conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_tables_and_their_rows(mysql_config):
    load_old_schema(mysql_config, 121)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _plan(conn, request_text="Make it an M planet")
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0121')")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT request_text FROM change_plans").fetchall() == [{"request_text": "Make it an M planet"}]
    finally:
        conn.close()


def test_a_plan_has_a_known_state_and_roles(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for state in ("draft", "requested", "proposed", "approved", "applied", "undone", "stale", "rejected", "withdrawn"):
            _plan(conn, state=state)
        row = conn.execute("SELECT state, approver_role, apply_job_id FROM change_plans ORDER BY id LIMIT 1").fetchone()
        assert row == {"state": "draft", "approver_role": None, "apply_job_id": None}
        with pytest.raises(Exception, match="chk_change_plans_state"):
            _plan(conn, state="pending")
        with pytest.raises(Exception, match="chk_change_plans_requester_role"):
            _plan(conn, role="owner")
        with pytest.raises(Exception, match="chk_change_plans_approver_role"):
            _plan(conn, approver_role="user")
        _plan(conn, role="user", state="requested", approver_role="moderator")
    finally:
        conn.close()


def test_the_operations_of_a_plan_are_in_order_and_go_with_the_plan(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        plan = _plan(conn)
        for order in (1, 2):
            conn.execute("INSERT INTO change_plan_ops (plan_id, op_order, kind, target_uid, before_json, after_json)"
                         " VALUES (%s, %s, 'set_class', %s, '{\"class\": \"A\"}', '{\"class\": \"M\"}')", (plan, order, UID))
        with pytest.raises(Exception, match="Duplicate entry"):
            conn.execute("INSERT INTO change_plan_ops (plan_id, op_order, kind) VALUES (%s, 2, 'set_class')", (plan,))
        with pytest.raises(Exception, match="chk_change_plan_ops_status"):
            conn.execute("INSERT INTO change_plan_ops (plan_id, op_order, kind, status) VALUES (%s, 3, 'x', 'done')", (plan,))
        row = conn.execute("SELECT status, JSON_EXTRACT(after_json, '$.class') AS c FROM change_plan_ops WHERE op_order = 1").fetchone()
        assert row["status"] == "pending" and row["c"] == '"M"'
        conn.execute("DELETE FROM change_plans WHERE id = %s", (plan,))
        assert conn.execute("SELECT COUNT(*) AS n FROM change_plan_ops").fetchone()["n"] == 0
    finally:
        conn.close()


def test_the_journal_is_keyed_by_object_and_sector_and_keeps_its_plan(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        plan = _plan(conn, state="applied")
        conn.execute("INSERT INTO edit_journal (plan_id, object_uid, sector_id, aspect, before_json, after_json, generator_epoch)"
                     " VALUES (%s, %s, 7, 'class', '\"A\"', '\"M\"', 44)", (plan, UID))
        found = conn.execute("SELECT aspect, generator_epoch FROM edit_journal WHERE object_uid = %s", (UID,)).fetchall()
        assert found == [{"aspect": "class", "generator_epoch": 44}]
        explain = conn.execute("EXPLAIN SELECT * FROM edit_journal WHERE object_uid = %s ORDER BY id", (UID,)).fetchone()
        assert explain["key"] == "idx_edit_journal_object"
        with pytest.raises(Exception, match="fk_edit_journal_plan"):
            conn.execute("DELETE FROM change_plans WHERE id = %s", (plan,))
        with pytest.raises(Exception, match="fk_edit_journal_plan"):
            conn.execute("INSERT INTO edit_journal (plan_id, object_uid, aspect, generator_epoch) VALUES (999999, %s, 'x', 1)", (UID,))
    finally:
        conn.close()
