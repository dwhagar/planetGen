# tests/test_migration_giant_interiors.py

"""
DB.89, GEN.642: the interior of a giant. A system saved and loaded keeps each giant's summary and layers; editing a body deletes the old
layers; migration 0154 gives a galaxy made at v153 the columns and tables, empty, and a body without them reads back as having no interior.
"""

import pytest

from planetgen.db import edits, store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import giant_interior
from planetgen.util import draw
from tests import db_schema_support

pytestmark = pytest.mark.fresh_database


def _system():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(p.body_type == "g" and p.moons for p in system.planets):
            return system, cfg
    pytest.fail("no system with a giant that has a moon in 60 seeds")


def _bodies(system):
    return [body for planet in system.planets if hasattr(planet, "moons") for body in [planet, *planet.moons]]


def _interior_of(body):
    return (tuple(getattr(body, name) for name in giant_interior.BODY_FIELDS), getattr(body, giant_interior.LAYERS_ATTRIBUTE))


def test_a_giant_has_an_interior_whose_layers_add_up_and_a_rocky_world_has_none():
    system, _cfg = _system()
    giants = [body for body in _bodies(system) if body.body_type == "g"]
    assert giants
    for body in giants:
        kind, layers = body.interior_kind, body.interior_layers
        assert kind in giant_interior.INTERIOR_KINDS and layers
        assert sum(layer["mass_earth"] for layer in layers) == pytest.approx(body.mass / 5.972e24, rel=1e-3)
        assert giant_interior.interior_of(body).total_mass_earth == pytest.approx(sum(layer["mass_earth"] for layer in layers))
    for body in _bodies(system):
        if body.body_type != "g":
            assert body.interior_kind is None and body.interior_layers is None


def test_a_saved_system_keeps_every_giants_interior(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before, after = [_interior_of(body) for body in _bodies(system)], [_interior_of(body) for body in _bodies(loaded)]
        assert any(summary[0] is not None for summary, _layers in before)
        assert after == before
        rows = conn.execute("SELECT COUNT(*) AS n FROM planet_interior_layers pl JOIN planets p ON p.id = pl.planet_id"
                            " WHERE p.star_system_id = ?", (system_id,)).fetchone()["n"]
        assert rows == sum(len(layers) for planet in system.planets if hasattr(planet, "moons")
                           for layers in [planet.interior_layers] if layers)
    finally:
        conn.rollback()
        conn.close()


def test_a_planet_cannot_hold_two_layers_at_one_index_or_an_unknown_kind(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        planet_id = conn.execute("SELECT planet_id AS id FROM planet_interior_layers pl JOIN planets p ON p.id = pl.planet_id"
                                 " WHERE p.star_system_id = ? LIMIT 1", (system_id,)).fetchone()["id"]
        values = "(?, 0, ?, 'molecular_fluid', 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)"
        sql = ("INSERT INTO planet_interior_layers (planet_id, layer_index, kind, state, top_pressure_gpa, bottom_pressure_gpa, mass_earth,"
               " top_radius_fraction, bottom_radius_fraction, hydrogen, helium, iron, rock, water_liquid, ice_x, nh3_ice, ch4_ice)"
               f" VALUES {values}")
        with pytest.raises(Exception):
            conn.execute(sql, (planet_id, "molecular_envelope"))
        conn.rollback()
        with pytest.raises(Exception):
            conn.execute(sql, (planet_id, "granite"))
    finally:
        conn.rollback()
        conn.close()


def test_a_v153_galaxy_gets_the_columns_and_tables_empty(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
    finally:
        conn.close()
    db_schema_support.downgrade_to(mysql_config, 153)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = conn.execute("SELECT interior_kind, core_mass_earth FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()
        assert rows and all(row["interior_kind"] is None and row["core_mass_earth"] is None for row in rows)
        assert conn.execute("SELECT COUNT(*) AS n FROM planet_interior_layers").fetchone()["n"] == 0
        loaded = store.load_star_system(conn, system_id)
        assert all(_interior_of(body) == ((None,) * len(giant_interior.BODY_FIELDS), None) for body in _bodies(loaded))
    finally:
        conn.close()
