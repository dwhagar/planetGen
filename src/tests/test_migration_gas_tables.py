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
from tests import db_schema_support

pytestmark = pytest.mark.fresh_database

_GASES = ("o2", "co2", "co", "n2", "ar", "h2", "h2o", "ch4", "h2s", "so2")


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


def test_a_v154_galaxy_gets_its_air_in_the_child_tables(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
    finally:
        conn.close()
    db_schema_support.downgrade_to(mysql_config, 154)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            conn.execute(f"UPDATE {table} SET p_o2_kpa = 21.2, p_n2_kpa = 78.8, p_ch4_kpa = 0 WHERE star_system_id = ?", (system_id,))
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        columns = {row["COLUMN_NAME"] for row in conn.execute(
            "SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME IN ('planets', 'moons')").fetchall()}
        assert not any(f"p_{gas}_kpa" in columns for gas in _GASES)
        rows = conn.execute("SELECT g.gas, g.partial_pa FROM planet_gases g JOIN planets p ON p.id = g.planet_id WHERE p.star_system_id = ?",
                            (system_id,)).fetchall()
        assert rows and {row["gas"] for row in rows} == {"O2", "N2"}
        assert {row["gas"]: row["partial_pa"] for row in rows}["O2"] == pytest.approx(21200.0)
        loaded = store.load_star_system(conn, system_id)
        assert all(air in ({}, {"O2": pytest.approx(21200.0), "N2": pytest.approx(78800.0)}) for air in _airs(loaded))
    finally:
        conn.close()
