# tests/test_migration_survival_equipment.py

"""
GEN.655, DB.93: migration 0136 swaps `planets.equipment_tier` and `moons.equipment_tier` for the survival equipment levels.
A v135 galaxy gains the columns (NULL until `planetgen rescore-habitability`), loses the tier, is detected as current, and a
stopped run can be started again.
"""

import pytest

from planetgen.db import store
from planetgen.physics import equipment
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _columns(config, table):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        return {row["Field"] for row in conn.execute(f"SHOW COLUMNS FROM {table}").fetchall()}
    finally:
        conn.close()


def test_a_v135_galaxy_gains_the_equipment_levels_and_loses_the_tier(mysql_config):
    load_old_schema(mysql_config, 135)
    assert store.schema_status(mysql_config)[0] == 135
    assert "equipment_tier" in _columns(mysql_config, "planets")
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    for table in ("planets", "moons"):
        columns = _columns(mysql_config, table)
        assert set(equipment.FIELDS) <= columns and "equipment_tier" not in columns
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()


def test_a_stopped_run_can_be_started_again(mysql_config):
    load_old_schema(mysql_config, 135)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("ALTER TABLE planets ADD COLUMN eq_body TINYINT NULL")  # the first step done, the rest not
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert set(equipment.FIELDS) <= _columns(mysql_config, "planets")
    assert set(equipment.FIELDS) <= _columns(mysql_config, "moons")
