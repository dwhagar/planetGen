# tests/test_migration_star_count.py

"""
GEN.267: migration 0134 adds `system_configs.stars`, `max_stars` and `prevalence_multi_star`. A config saved and loaded keeps
its forced count, its cap and its multi-star prevalence; a config row made at v133 reads back with no forced count, the cap
of 12 and the usual chance, as it always generated.
"""

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.generation.config import SystemConfig
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def test_a_saved_config_keeps_its_count_its_cap_and_its_prevalence(mysql_config):
    config = SystemConfig()
    config.STARS, config.MAX_STARS = 6, 24
    config.PREVALENCE = {"multi_star": 75.0, "comets": 10.0}
    conn = store.get_connection(mysql_config)
    try:
        config_id = store.insert_system_config(conn, config)
        loaded = store.load_system_config(conn, config_id)
        assert (loaded.STARS, loaded.MAX_STARS) == (6, 24)
        assert loaded.PREVALENCE == {"multi_star": 75.0, "comets": 10.0}
        plain_id = store.insert_system_config(conn, SystemConfig())
        plain = store.load_system_config(conn, plain_id)
        assert (plain.STARS, plain.MAX_STARS, plain.PREVALENCE) == (None, tuning.MAX_STARS_DEFAULT, {})
    finally:
        conn.rollback()
        conn.close()


def test_a_v133_galaxy_keeps_its_configs_with_no_count_the_default_cap_and_the_usual_chance(mysql_config):
    load_old_schema(mysql_config, 133)
    assert store.schema_status(mysql_config)[0] == 133
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("INSERT INTO system_configs (binary_system, num_orbits) VALUES (1, 4)")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT stars, max_stars, prevalence_multi_star FROM system_configs").fetchone()
        assert (row["stars"], row["max_stars"], row["prevalence_multi_star"]) == (None, 12, 0)
        config = store.load_system_config(conn, conn.execute("SELECT id FROM system_configs").fetchone()["id"])
        assert (config.STARS, config.MAX_STARS, config.BINARY_SYSTEM, config.NUM_ORBITS) == (None, 12, True, 4)
        assert store.detect_schema_version(conn) == store.SCHEMA_VERSION
    finally:
        conn.close()
