# tests/test_migration_gas_tables.py

"""
DB.71, GEN.424: the air as child table rows. A system saved and loaded keeps every body's gases, the balance gas included; migration 0155
gives a galaxy made at v154 the rows from its ten `p_<gas>_kpa` columns and drops the columns.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import atmosphere
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


def _airs(system):
    return [atmosphere.gas_partials_pa(body) for planet in system.planets if hasattr(planet, "moons") for body in [planet, *planet.moons]]


def test_a_saved_system_keeps_every_bodys_air(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before, after = _airs(system), _airs(loaded)
        assert any(before) and after == before
        giants = [p for p in system.planets if p.body_type == "g"]
        assert giants and all(set(atmosphere.gas_partials_pa(g)) >= {"H2", "He"} for g in giants)
        stored = sum(conn.execute(f"SELECT COUNT(*) AS n FROM {kind}_gases g JOIN {kind}s b ON b.id = g.{kind}_id WHERE b.star_system_id = ?",
                                  (system_id,)).fetchone()["n"] for kind in ("planet", "moon"))
        assert stored == sum(len(air) for air in before)
    finally:
        conn.rollback()
        conn.close()
