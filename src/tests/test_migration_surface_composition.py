# tests/test_migration_surface_composition.py

"""
GEN.441: the surface composition columns. A system saved and loaded keeps every rocky body's ground; migration 0159 gives a galaxy made at
v151 the columns, empty, and a body without them reads back as having no ground.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import surface_composition as sc
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


def _grounds(system):
    return [{field: getattr(body, field) for field in sc.COLUMN_FIELDS}
            for planet in system.planets if hasattr(planet, "moons") for body in [planet, *planet.moons]]


def test_a_saved_system_keeps_every_bodys_ground(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before, after = _grounds(system), _grounds(loaded)
        assert any(item["ground_type"] for item in before) and any(item["ground_type"] is None for item in before)
        assert after == before
        for item in before:
            if item["ground_type"]:
                assert sum(item[field] for field in sc.FRACTION_FIELDS) == pytest.approx(1.0)
    finally:
        conn.rollback()
        conn.close()
