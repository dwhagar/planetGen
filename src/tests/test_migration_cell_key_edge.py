# tests/test_migration_cell_key_edge.py

"""
DB.65: migration 0106 recomputes every generated `cell_key` with the expression that holds the cube number at 2^21 - 1,
so a point just under the reach no longer takes the opposite edge's key. A migration run again leaves the columns alone.
"""

import importlib
import math

import pytest

from planetgen.db import store
from planetgen.galaxy import cellkey
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

migration = importlib.import_module("planetgen.db.migrations.versions.0106_cell_key_edge")

JUST_UNDER = math.nextafter(cellkey.REACH_PC, 0.0)


def _insert_black_hole(conn, x):
    conn.execute("INSERT INTO black_holes (name, mass_solar, event_horizon_radius_km, spin, has_accretion_disk, temperature_k,"
                 " luminosity_w, age_gy, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc)"
                 " VALUES (?, 10, 30, 0, 0, 0, 0, 1, ?, 0, 0, ?)",
                 ("edge", x, x))
    conn.commit()


def _key(conn):
    return conn.execute("SELECT cell_key AS k FROM black_holes WHERE name = 'edge'").fetchone()["k"]


def test_a_v105_point_at_the_edge_takes_the_last_cube_after_the_migration(mysql_config):
    load_old_schema(mysql_config, 105)
    assert store.schema_status(mysql_config)[0] == 105
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _insert_black_hole(conn, JUST_UNDER)
        wrapped = _key(conn)
        assert wrapped != cellkey.cell_key(JUST_UNDER, 0.0, 0.0)  # the defect: the x axis lost its top bit
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _key(conn) == cellkey.cell_key(JUST_UNDER, 0.0, 0.0)
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION      # 106 at DB.65; later migrations move it on
        for table, _ in migration.CELL_TABLES:
            expression = conn.execute(
                "SELECT generation_expression AS e FROM information_schema.columns WHERE table_schema = DATABASE()"
                " AND table_name = ? AND column_name = 'cell_key'", (table,)).fetchone()["e"]
            assert "least(" in expression.lower(), table
    finally:
        conn.close()


def test_the_listed_tables_are_every_table_with_a_cell_key(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        held = {row["t"] for row in conn.execute(
            "SELECT table_name AS t FROM information_schema.columns WHERE table_schema = DATABASE()"
            " AND column_name = 'cell_key'").fetchall()}
    finally:
        conn.close()
    assert held == {table for table, _ in migration.CELL_TABLES}


def test_a_migration_run_again_leaves_the_columns_as_they_are(mysql_config):
    load_old_schema(mysql_config, 105)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _insert_black_hole(conn, JUST_UNDER)
        before = _key(conn)
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0105')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _key(conn) == before
    finally:
        conn.close()
