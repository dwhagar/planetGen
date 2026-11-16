# tests/test_migration_regional_temperature.py

"""
GEN.656: migration 0137 gives `planets` and `moons` the regional temperature columns. A v136 galaxy gains them (NULL
until `planetgen rescore-habitability`), is detected as current, and a stopped run can be started again.
"""

import pytest

from planetgen.db import store
from planetgen.physics import regional_temperature
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _columns(config, table):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        return {row["Field"] for row in conn.execute(f"SHOW COLUMNS FROM {table}").fetchall()}
    finally:
        conn.close()


def test_a_stopped_run_can_be_started_again(mysql_config):
    load_old_schema(mysql_config, 136)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("ALTER TABLE planets ADD COLUMN temp_polar_k DOUBLE NULL")  # the first step done, the rest not
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    for table in ("planets", "moons"):
        assert set(regional_temperature.FIELDS) <= _columns(mysql_config, table)
