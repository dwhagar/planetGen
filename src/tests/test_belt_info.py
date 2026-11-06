# tests/test_belt_info.py

"""
UX.241: the information page data of a stored belt (`describe/belt_info.py`): its stored members by designation, the size table of the rest of its mass and the density factor.
"""

import pytest

from planetgen.db import store
from planetgen.describe import belt_info
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw


def _config():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = True
    cfg.MOONS = False
    cfg.BINARY_SYSTEM = False
    return cfg


def _saved_belt(mysql_config):
    config = _config()
    for seed in range(1, 40):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=config)
        if not any(b.body_type == 'a' and b.members for b in system.planets):
            continue
        saved = store.save_system(system, config, config=mysql_config)
        with store.get_connection(mysql_config) as conn:
            loaded = store.load_star_system(conn, saved)
            belt = next(b for b in loaded.planets if b.body_type == 'a' and b.members)
            return conn.execute("SELECT id FROM asteroid_belts WHERE id = ?", (belt.db_id,)).fetchone()["id"], belt, mysql_config
    pytest.fail("no belt with stored members in 39 systems")


def test_a_stored_belt_has_a_page_of_its_members_its_bins_and_its_density(mysql_config):
    belt_id, belt, config = _saved_belt(mysql_config)
    with store.get_connection(config) as conn:
        page = belt_info.belt_info(conn, belt_id)
    assert [b["rank"] for b in page["bodies"]] == [m["rank"] for m in belt.members]
    assert all(b["designation"].count(".") == 1 for b in page["bodies"])
    assert page["bodies"][0]["mass"]["value"] == pytest.approx(belt.members[0]["mass_kg"])
    assert page["belt"]["star"] and page["belt"]["mass"]["unit"] == "kg"
    assert [f["above"]["value"] for f in page["density_factor"]["by_diameter"]] == [1.0, 10.0, 100.0]
    totals = page["totals"]
    assert abs(totals["difference"]["value"]) <= 1e-6 * totals["belt_mass"]["value"]       # the listed bodies and the bins add up to the belt


def test_a_belt_that_is_not_there_has_no_page(mysql_config):
    with store.get_connection(mysql_config) as conn:
        with pytest.raises(belt_info.BeltPageNotFound):
            belt_info.belt_info(conn, 999999999)
