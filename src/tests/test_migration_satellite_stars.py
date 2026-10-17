# tests/test_migration_satellite_stars.py

"""
GEN.324: migration 0105 lets `bright_stars.population` take `satellite`, the stars of a satellite galaxy that the
component pass places. No row is touched, and a migration run again leaves the check as it is.
"""

import importlib

import pytest

from planetgen.db import store
from planetgen.generation.bright_stars import POPULATIONS
from tests.db_schema_support import load_old_schema
from tests.test_migration_halo_populations import _insert_star

pytestmark = pytest.mark.fresh_database

migration = importlib.import_module("planetgen.db.migrations.versions.0105_satellite_stars")


def test_a_v104_galaxy_takes_satellite_stars_after_the_migration(mysql_config):
    load_old_schema(mysql_config, 104)
    assert store.schema_status(mysql_config)[0] == 104
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _insert_star(conn, "cluster")
        conn.commit()
        with pytest.raises(Exception, match="chk_bright_stars_population"):
            _insert_star(conn, "satellite")
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 1
        for population in migration.POPULATIONS:
            _insert_star(conn, population)
        with pytest.raises(Exception, match="chk_bright_stars_population"):
            _insert_star(conn, "disk")
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_check_as_it_is(mysql_config):
    load_old_schema(mysql_config, 104)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0104')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        clause = conn.execute(
            "SELECT check_clause AS c FROM information_schema.check_constraints WHERE constraint_schema = DATABASE()"
            " AND constraint_name = 'chk_bright_stars_population'").fetchone()["c"]
        assert clause.count("_utf8mb4") == len(POPULATIONS)    # the later migrations (satellite ages) add names
    finally:
        conn.close()
