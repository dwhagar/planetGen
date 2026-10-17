# tests/test_migration_satellite_age_populations.py

"""
GEN.526: migration 0107 lets `bright_stars.population` take `satellite_young` and `satellite_intermediate`, the stars a
satellite's history puts in the last 0.1 Gy and the 3 Gy before that, and widens the column for the longer name. No row
is touched, and a migration run again leaves the check as it is.
"""

import importlib

import pytest

from planetgen.db import store
from planetgen.generation.bright_stars import POPULATIONS
from tests.db_schema_support import load_old_schema
from tests.test_migration_halo_populations import _insert_star

pytestmark = pytest.mark.fresh_database

migration = importlib.import_module("planetgen.db.migrations.versions.0107_satellite_age_populations")


def test_a_v106_galaxy_takes_the_two_new_populations_after_the_migration(mysql_config):
    load_old_schema(mysql_config, 106)
    assert store.schema_status(mysql_config)[0] == 106
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _insert_star(conn, "satellite")
        conn.commit()
        with pytest.raises(Exception, match="chk_bright_stars_population|Data too long"):
            _insert_star(conn, "satellite_young")
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
    load_old_schema(mysql_config, 106)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0106')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        clause = conn.execute(
            "SELECT check_clause AS c FROM information_schema.check_constraints WHERE constraint_schema = DATABASE()"
            " AND constraint_name = 'chk_bright_stars_population'").fetchone()["c"]
        assert clause.count("_utf8mb4") == len(POPULATIONS) == len(migration.POPULATIONS)
    finally:
        conn.close()
