# tests/test_migration_star_members.py

"""
GEN.630: migration 0133 lets `stars.role` be `member`, adds `stars.slot`, and gives a stored pair its slots (the
primary 0, the secondary 1); a single star keeps none. It can be run again.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import insert_filled_row as _insert, load_old_schema, table_columns as _columns

pytestmark = pytest.mark.fresh_database


def _old_database_with_stars(mysql_config):
    load_old_schema(mysql_config, 132)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        for system_id, configuration in ((1, "wide"), (2, None)):
            _insert(conn, "star_systems", id=system_id, binary_configuration=configuration,
                    is_binary=int(bool(configuration)))
        _insert(conn, "stars", star_system_id=1, role="primary", parent_node_id=0)
        _insert(conn, "stars", star_system_id=1, role="secondary", parent_node_id=0)
        _insert(conn, "stars", star_system_id=2, role="single")
        conn.execute("SET FOREIGN_KEY_CHECKS = 1")
        conn.commit()
        with pytest.raises(Exception):
            _insert(conn, "stars", star_system_id=1, role="member", parent_node_id=1)
    finally:
        conn.close()


def test_a_v132_database_takes_a_member_and_keeps_its_pair(mysql_config):
    _old_database_with_stars(mysql_config)
    assert store.schema_status(mysql_config)[0] == 132
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "slot" in _columns(conn, "stars")
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
        rows = {row["role"]: row["slot"] for row in conn.execute("SELECT role, slot FROM stars").fetchall()}
        assert rows == {"primary": 0, "secondary": 1, "single": None}
        _insert(conn, "stars", star_system_id=1, role="member", parent_node_id=1, slot=1)
        with pytest.raises(Exception, match="stars_chk_1"):
            _insert(conn, "stars", star_system_id=1, role="tertiary")
        conn.rollback()
    finally:
        conn.close()


def test_a_migration_run_again_changes_nothing(mysql_config):
    _old_database_with_stars(mysql_config)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0132')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert sorted(row["slot"] for row in conn.execute(
            "SELECT slot FROM stars WHERE role <> 'single'").fetchall()) == [0, 1]
    finally:
        conn.close()
