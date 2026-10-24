# tests/test_migration_life_scores.py

"""
GEN.649: a saved system keeps its four life scores; migration 0150 gives a galaxy made at v149 the scores worked out from the
two likelihoods and the evolutionary path its rows stored, and drops the three old habitability scores.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import life_scores
from planetgen.util import draw

pytestmark = pytest.mark.fresh_database

OLD_COLUMNS = ("phi_bio", "phi_cpx", "phi_tech")


def _system():
    for seed in range(60):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if sum(1 for p in system.planets if p.body_type == "t") >= 2:
            return system, cfg
    pytest.fail("no system with two rocky planets in 60 seeds")


def _downgrade(config):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            conn.execute(f"ALTER TABLE {table} DROP COLUMN life_microbial, DROP COLUMN life_flora, DROP COLUMN life_fauna,"
                         " DROP COLUMN life_sentient, ADD COLUMN phi_bio DOUBLE, ADD COLUMN phi_cpx DOUBLE,"
                         " ADD COLUMN phi_tech DOUBLE, ALGORITHM=COPY")
        conn.execute("UPDATE alembic_version SET version_num = '0149'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (149)")
        conn.commit()
        assert store._stored_version(conn, "schema_migrations") == 149
    finally:
        conn.close()


def test_a_saved_system_keeps_its_life_scores(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        rows = conn.execute("SELECT life_microbial, life_flora, life_fauna, life_sentient, phi4, body_type FROM planets"
                            " WHERE star_system_id = ?", (system_id,)).fetchall()
        assert any(row["phi4"] is not None for row in rows)
        for row in rows:
            scores = [row[name] for name in life_scores.FIELDS]
            assert all(value is None for value in scores) == (row["phi4"] is None)
        loaded = store.load_star_system(conn, system_id)
        for planet, original in zip((p for p in loaded.planets if p.body_type != "a"),
                                    (p for p in system.planets if p.body_type != "a")):
            for name in life_scores.FIELDS:
                assert getattr(planet, name) == getattr(original, name)
    finally:
        conn.rollback()
        conn.close()


def test_a_v149_galaxy_gets_its_scores_from_the_stored_likelihoods_and_path(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        ids = [row["id"] for row in conn.execute(
            "SELECT id FROM planets WHERE star_system_id = ? AND phi4 IS NOT NULL ORDER BY id", (system_id,)).fetchall()]
    finally:
        conn.close()
    assert len(ids) >= 2
    _downgrade(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        with conn:
            # a world with a civilization, one with microbes only and one without life, as a v149 row would hold them
            conn.execute("UPDATE planets SET phi_bio = 0.8, phi_cpx = 0.6, life_chemical = 'chlorophyll-a',"
                         " life_stage = 'technological_civilization', evolutionary_pace = 'Earth Norm',"
                         " life_milestone_age_gyr = 4.0, life_system_age_gyr = 4.6 WHERE id = ?", (ids[0],))
            conn.execute("UPDATE planets SET phi_bio = 0.4, phi_cpx = 0.0, life_chemical = 'chlorophyll-a',"
                         " life_stage = NULL, evolutionary_pace = NULL, life_milestone_age_gyr = NULL,"
                         " life_system_age_gyr = NULL WHERE id = ?", (ids[1],))
            conn.execute("UPDATE planets SET phi_tech = 0.7 WHERE id = ?", (ids[1],))
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        civilization, microbes = (conn.execute(
            "SELECT life_microbial, life_flora, life_fauna, life_sentient FROM planets WHERE id = ?", (i,)).fetchone()
            for i in ids[:2])
        assert tuple(civilization[name] for name in life_scores.FIELDS) == life_scores.scores(0.8, 0.6, "technological_civilization", 4.0, 4.6)
        assert civilization["life_sentient"] > 0
        assert tuple(microbes[name] for name in life_scores.FIELDS) == (4, 0, 0, 0)  # life with no timeline is microbial only
        for table in ("planets", "moons"):
            columns = {row["Field"] for row in conn.execute(f"SHOW COLUMNS FROM {table}").fetchall()}
            assert not columns & set(OLD_COLUMNS)
            assert set(life_scores.FIELDS) <= columns
    finally:
        conn.close()
