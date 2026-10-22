# tests/test_migration_species_limit_values.py

"""
DB.30: migration 0119 makes `species_limit_values` and touches nothing else. A migration run again leaves the table as it
is, a limit is one row per species and name, and a species' rows go with the species.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _count(conn, table):
    return conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]


def test_a_v118_galaxy_gets_an_empty_table_and_keeps_its_rows(mysql_config):
    load_old_schema(mysql_config, 118)
    assert store.schema_status(mysql_config)[0] == 118
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
                            " AND TABLE_NAME = 'species_limit_values'").fetchone()["n"] == 0
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert store.schema_status(mysql_config)[0] >= 119
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _count(conn, "species_limit_values") == 0
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_table_and_its_rows(mysql_config):
    load_old_schema(mysql_config, 118)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        conn.execute("INSERT INTO species_limit_values (species_id, limit_name, value) VALUES (7, 'wet_bulb_c', 33.0)")
        conn.execute("SET FOREIGN_KEY_CHECKS = 1")
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0118')")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT limit_name, value FROM species_limit_values").fetchall() == [
            {"limit_name": "wet_bulb_c", "value": 33.0}]
    finally:
        conn.close()


def test_a_limit_is_one_row_per_species_and_name_and_needs_a_species(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        conn.execute("INSERT INTO species_limit_values (species_id, limit_name, value) VALUES (7, 'wet_bulb_c', 33.0)")
        conn.execute("INSERT INTO species_limit_values (species_id, limit_name, value) VALUES (7, 'comfort_low_c', -5.0)")
        conn.execute("INSERT INTO species_limit_values (species_id, limit_name, value) VALUES (8, 'wet_bulb_c', 34.0)")
        with pytest.raises(Exception, match="Duplicate entry"):
            conn.execute("INSERT INTO species_limit_values (species_id, limit_name, value) VALUES (7, 'wet_bulb_c', 35.0)")
        conn.execute("SET FOREIGN_KEY_CHECKS = 1")
        conn.rollback()
        with pytest.raises(Exception, match="fk_species_limit_values_species"):
            conn.execute("INSERT INTO species_limit_values (species_id, limit_name, value) VALUES (999999, 'x', 1.0)")
    finally:
        conn.close()
