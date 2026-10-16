# tests/test_migration_cell_resolution.py

"""
GEN.281: migration 0098 adds `cell_star_ceiling` and `fill_luminosity_cap_sol` to `galaxy_shape`, at the tunings'
defaults, and can be run again.
"""

import pytest

from planetgen import tuning
from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _columns(conn):
    return {row["COLUMN_NAME"] for row in conn.execute(
        "SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'galaxy_shape'").fetchall()}


def test_a_v97_galaxy_gets_the_two_settings_at_their_defaults(mysql_config):
    load_old_schema(mysql_config, 97)
    assert store.schema_status(mysql_config)[0] == 97
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert {"cell_star_ceiling", "fill_luminosity_cap_sol"} <= _columns(conn)
        stored = store.get_galaxy_shape(conn)
    finally:
        conn.close()
    assert stored is None or (stored.cell_star_ceiling, stored.fill_luminosity_cap_sol) == (
        tuning.CELL_STAR_CEILING, tuning.FILL_LUMINOSITY_CAP_SOL)


def test_a_migration_run_again_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 97)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0097')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert {"cell_star_ceiling", "fill_luminosity_cap_sol"} <= _columns(conn)
    finally:
        conn.close()
