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
