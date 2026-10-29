# tests/test_migration_rogue_habitability.py

"""
A saved rogue planet keeps its habitability score, and migration 0158 gives a v157 galaxy the score worked out from the surface
conditions its rows already hold (nothing is drawn).
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.physics import rogue_habitability
from planetgen.util import draw
from tests import db_schema_support

pytestmark = pytest.mark.fresh_database


def _rogue(mass_bin="terrestrial"):
    for seed in range(40):
        with draw.bound(seed):
            planet = RoguePlanet(SystemConfig(), mass_bin=mass_bin)
        if planet.planet_type == "t":
            return planet
    pytest.fail("no terrestrial rogue planet in 40 seeds")


def _downgrade(config):
    db_schema_support.downgrade_to(config, 157)


def test_a_saved_rogue_planet_keeps_its_score(mysql_config):
    planet = _rogue()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = store.insert_rogue_planet(conn, planet)
        row = conn.execute("SELECT * FROM rogue_planets WHERE id = ?", (row_id,)).fetchone()
        for name in rogue_habitability.ROGUE_FIELDS:
            assert row[name] == pytest.approx(getattr(planet, name)), name
    finally:
        conn.rollback()
        conn.close()


def test_a_gas_giant_keeps_no_score(mysql_config):
    for seed in range(60):
        with draw.bound(seed):
            planet = RoguePlanet(SystemConfig(), mass_bin="saturn")
        if planet.planet_type == "g":
            break
    else:
        pytest.skip("no gas giant in 60 seeds")
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = store.insert_rogue_planet(conn, planet)
        row = conn.execute("SELECT * FROM rogue_planets WHERE id = ?", (row_id,)).fetchone()
        assert all(row[name] is None for name in rogue_habitability.ROGUE_FIELDS)
    finally:
        conn.rollback()
        conn.close()


def test_a_v157_galaxy_gets_the_score_back_from_the_stored_surface(mysql_config):
    terrestrial, giant = _rogue(), None
    for seed in range(60):
        with draw.bound(seed):
            candidate = RoguePlanet(SystemConfig(), mass_bin="saturn")
        if candidate.planet_type == "g":
            giant = candidate
            break
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            ids = [store.insert_rogue_planet(conn, planet) for planet in (terrestrial, giant) if planet is not None]
    finally:
        conn.close()
    _downgrade(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT * FROM rogue_planets WHERE id = ?", (ids[0],)).fetchone()
        for name in rogue_habitability.ROGUE_FIELDS:
            assert row[name] == pytest.approx(getattr(terrestrial, name)), name
        if giant is not None:
            row = conn.execute("SELECT * FROM rogue_planets WHERE id = ?", (ids[1],)).fetchone()
            assert all(row[name] is None for name in rogue_habitability.ROGUE_FIELDS)
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION  # a second run changes nothing
