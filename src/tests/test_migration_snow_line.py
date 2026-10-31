# tests/test_migration_snow_line.py

"""
GEN.458: the stored snow line and the born-beyond flag. A system saved and loaded keeps its star's line and every body's flag;
migration 0132 gives a galaxy made at v131 the blackbody position for each star and the flag from each body's distance.
"""


import math

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants
from planetgen.util import draw
from tests import db_schema_support

pytestmark = pytest.mark.fresh_database


def _system_with_moons():
    for seed in range(60):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        flags = {p.born_beyond_snow_line for p in system.planets if p.body_type != 'a'}
        if flags == {True, False} and any(p.moons for p in system.planets):
            return system, cfg
    pytest.fail("no system with planets on both sides of the line and a moon in 60 seeds")


def test_a_saved_system_keeps_its_line_and_its_flags(mysql_config):
    system, cfg = _system_with_moons()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        star = conn.execute("SELECT snow_line_km FROM stars WHERE star_system_id = ?", (system_id,)).fetchone()
        assert star["snow_line_km"] == pytest.approx(system.star.snow_line_au * constants.AU_TO_KM)
        loaded = store.load_star_system(conn, system_id)
        assert loaded.star.snow_line_au == pytest.approx(system.star.snow_line_au)
        before = [(p.born_beyond_snow_line, [m.born_beyond_snow_line for m in p.moons]) for p in system.planets if p.body_type != 'a']
        after = [(p.born_beyond_snow_line, [m.born_beyond_snow_line for m in p.moons]) for p in loaded.planets if p.body_type != 'a']
        assert after == before
        moons = conn.execute("SELECT born_beyond_snow_line FROM moons WHERE star_system_id = ?", (system_id,)).fetchall()
        assert moons and all(row["born_beyond_snow_line"] in (0, 1) for row in moons)
    finally:
        conn.rollback()
        conn.close()


def test_a_v131_galaxy_gets_the_blackbody_line_and_the_flag_from_the_distance(mysql_config):
    system, cfg = _system_with_moons()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
    finally:
        conn.close()
    # make it a v131 galaxy again: the columns gone, the rows kept
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM galaxy_block_stats WHERE m > 65535")      # a block of an older galaxy never held more (MAP.258, v146)
        conn.commit()
    finally:
        conn.close()
    db_schema_support.downgrade_to(mysql_config, 131)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        star = conn.execute("SELECT snow_line_km, luminosity_w FROM stars WHERE star_system_id = ?", (system_id,)).fetchone()
        # the migration is frozen at the 170 K it was written with (the generator moved to 182 K in GEN.950): an older galaxy keeps its line
        centre_km = math.sqrt(star["luminosity_w"] / (16.0 * math.pi * constants.STEFAN_BOLTZMANN_CONSTANT)) / 170.0 ** 2 / 1000.0
        assert star["snow_line_km"] == pytest.approx(centre_km, rel=1e-9)
        planets = conn.execute("SELECT id, distance_km, born_beyond_snow_line FROM planets WHERE star_system_id = ?",
                               (system_id,)).fetchall()
        assert planets
        for planet in planets:
            assert bool(planet["born_beyond_snow_line"]) == (planet["distance_km"] >= centre_km)
        moons = conn.execute(
            "SELECT m.born_beyond_snow_line AS moon, p.born_beyond_snow_line AS planet FROM moons m JOIN planets p ON p.id = m.planet_id"
            " WHERE m.star_system_id = ?", (system_id,)).fetchall()
        assert moons and all(row["moon"] == row["planet"] for row in moons)
        null = conn.execute("SELECT is_nullable FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
                            " AND TABLE_NAME = 'stars' AND COLUMN_NAME = 'snow_line_km'").fetchone()
        assert next(iter(null.values())) == "NO"
    finally:
        conn.close()
