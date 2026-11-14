# tests/test_migration_bulk_composition.py

"""
DB.69, GEN.428: the bulk composition columns. A system saved and loaded keeps every body's composition; migration 0144 gives a galaxy made at
v143 the columns, empty, and a body without them reads back as having no composition.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import bulk_composition as bc
from planetgen.util import draw

pytestmark = pytest.mark.fresh_database


def _system():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(getattr(p, "moons", None) for p in system.planets) and any(p.body_type == "g" for p in system.planets):
            return system, cfg
    pytest.fail("no system with a giant and a moon in 60 seeds")


def _compositions(system):
    return [bc.bulk_fractions(body) for planet in system.planets if hasattr(planet, "moons") for body in [planet, *planet.moons]]


def test_a_saved_system_keeps_every_bodys_composition(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before, after = _compositions(system), _compositions(loaded)
        assert before and all(item is not None for item in before)
        assert after == before
        row = conn.execute("SELECT carbon_rock, composition_origin FROM moons WHERE star_system_id = ? LIMIT 1", (system_id,)).fetchone()
        assert row["carbon_rock"] in (0, 1) and row["composition_origin"] in bc.ORIGINS
    finally:
        conn.rollback()
        conn.close()
