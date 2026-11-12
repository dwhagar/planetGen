# tests/test_migration_snow_line.py

"""
GEN.458: the stored snow line and the born-beyond flag. A system saved and loaded keeps its star's line and every body's flag;
migration 0132 gives a galaxy made at v131 the blackbody position for each star and the flag from each body's distance.
"""


import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants
from planetgen.util import draw

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
