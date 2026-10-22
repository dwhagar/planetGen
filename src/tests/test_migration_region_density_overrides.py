# tests/test_migration_region_density_overrides.py

"""
DB.38: migration 0124 makes `region_density_overrides` and touches nothing else. A migration run again leaves the table and its
rows as they are, a row names a known shape and a multiplier above zero, and undoing is the one change a row gets.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _exists(conn):
    return conn.execute("SELECT COUNT(*) AS n FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
                        " AND TABLE_NAME = 'region_density_overrides'").fetchone()["n"]


def _override(conn, shape="sphere", multiplier=2, **more):
    columns = {"shape": shape, "params": '{"center_pc": [0, 0, 0], "radius_pc": 50}', "multiplier": multiplier,
               "set_by_name": "Boss", **more}
    return conn.execute(f"INSERT INTO region_density_overrides ({', '.join(columns)}) VALUES ({', '.join(['%s'] * len(columns))})",
                        tuple(columns.values())).lastrowid


def test_a_v123_galaxy_gets_an_empty_table_and_keeps_its_rows(mysql_config):
    load_old_schema(mysql_config, 123)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _exists(conn) == 0
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _exists(conn) == 1
        assert conn.execute("SELECT COUNT(*) AS n FROM region_density_overrides").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_table_and_its_rows(mysql_config):
    load_old_schema(mysql_config, 123)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _override(conn, multiplier=0.5, reason="thin the halo")
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0123')")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT multiplier, reason FROM region_density_overrides").fetchone()
        assert (float(row["multiplier"]), row["reason"]) == (0.5, "thin the halo")
    finally:
        conn.close()


def test_a_row_names_a_known_shape_and_a_multiplier_above_zero(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for shape in ("sphere", "box", "shell", "band", "wedge", "sectors"):
            _override(conn, shape=shape)
        _override(conn, multiplier=12.5)  # beyond 10 is the interface's warning, not the table's refusal
        with pytest.raises(Exception, match="chk_region_density_overrides_shape"):
            _override(conn, shape="cube")
        with pytest.raises(Exception, match="chk_region_density_overrides_multiplier"):
            _override(conn, multiplier=0)
        with pytest.raises(Exception, match="chk_region_density_overrides_multiplier"):
            _override(conn, multiplier=-1)
    finally:
        conn.close()


def test_the_settings_in_force_are_those_not_undone_and_undoing_stamps_one_row(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        first, second = _override(conn, multiplier=2), _override(conn, multiplier=0.5)
        conn.execute("UPDATE region_density_overrides SET undone_at = CURRENT_TIMESTAMP(3) WHERE id = %s", (first,))
        active = conn.execute("SELECT id FROM region_density_overrides WHERE undone_at IS NULL ORDER BY id").fetchall()
        assert active == [{"id": second}]
        row = conn.execute("SELECT set_at IS NOT NULL AS stamped, params->>'$.radius_pc' AS radius FROM region_density_overrides"
                           " WHERE id = %s", (second,)).fetchone()
        assert row == {"stamped": 1, "radius": "50"}
    finally:
        conn.close()
