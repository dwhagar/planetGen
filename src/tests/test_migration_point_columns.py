# tests/test_migration_point_columns.py

"""
DB.55: migration 0102 adds the point columns the layouts of `db/point_columns.py` say an object lacks, backfills them
from what is stored, and can be run again.
"""

import importlib
import math

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

MIGRATION = importlib.import_module("planetgen.db.migrations.versions.0102_point_columns")


def _columns(conn, table):
    return {row["COLUMN_NAME"] for row in conn.execute(
        "SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?",
        (table,)).fetchall()}


def _insert(conn, table, **given):
    """A row of `table` with every NOT NULL column that has no default filled with a neutral value."""
    rows = conn.execute(
        "SELECT COLUMN_NAME, DATA_TYPE, IS_NULLABLE, COLUMN_DEFAULT, EXTRA FROM information_schema.COLUMNS"
        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", (table,)).fetchall()
    values = dict(given)
    for row in rows:
        name = row["COLUMN_NAME"]
        if name in values or row["IS_NULLABLE"] == "YES" or row["COLUMN_DEFAULT"] is not None or "auto_increment" in row["EXTRA"]:
            continue
        kind = row["DATA_TYPE"]
        values[name] = b"\x00" * 12 if kind in ("binary", "varbinary", "blob") else (
            "" if kind in ("char", "varchar", "text") else 0)
    names = ", ".join(f"`{name}`" for name in values)
    conn.execute(f"INSERT INTO {table} ({names}) VALUES ({', '.join('?' for _ in values)})", tuple(values.values()))


def test_a_v101_database_gets_every_new_column(mysql_config):
    load_old_schema(mysql_config, 101)
    assert store.schema_status(mysql_config)[0] == 101
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table, columns in MIGRATION.ADDITIONS.items():
            assert {name for name, _ in columns} <= _columns(conn, table), table
        assert {name for name, _ in MIGRATION.SECTOR_COLUMNS} <= _columns(conn, "sectors")
        assert store.detect_schema_version(conn) == store.SCHEMA_VERSION
    finally:
        conn.close()


def test_a_migration_run_again_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 101)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0101')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION


def test_a_black_hole_takes_the_tangent_at_its_old_speed(mysql_config):
    load_old_schema(mysql_config, 101)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _insert(conn, "black_holes", center_x_pc=3000.0, center_y_pc=4000.0, center_z_pc=0.0, galactic_radius_pc=5000.0,
                galactic_orbital_speed_kms=200.0, mass_class="stellar")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT velocity_x_kms, velocity_y_kms, velocity_z_kms FROM black_holes").fetchone()
    finally:
        conn.close()
    assert math.isclose(row["velocity_x_kms"], -160.0) and math.isclose(row["velocity_y_kms"], 120.0)
    assert row["velocity_z_kms"] == 0


def test_a_field_of_unknown_mass_takes_a_stated_zero_mass_and_vector(mysql_config):
    load_old_schema(mysql_config, 101)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _insert(conn, "asteroid_fields", center_x_pc=1.0, center_y_pc=0.0, center_z_pc=0.0, galactic_radius_pc=1.0,
                galactic_orbital_speed_kms=0.0, density="typical")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT mass_kg, velocity_x_kms, velocity_y_kms FROM asteroid_fields").fetchone()
    finally:
        conn.close()
    assert row["mass_kg"] == 0 and row["velocity_x_kms"] == 0 and row["velocity_y_kms"] == 0
