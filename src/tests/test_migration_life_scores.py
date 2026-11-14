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
