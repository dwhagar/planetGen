# tests/test_rogue_topup_build.py

"""GEN.1120: the rogues a sector's top-up plans are built through the planet chain with no star, placed in the sector and marked as an apparent origin."""

import math

import pytest

from planetgen import tuning
from planetgen.galaxy import spherical
from planetgen.generation import rogue_topup_build
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, rogue_topup
from planetgen.util import draw

EDGE_PC = 4.0
ADDRESS = (2609, 7, 32)


@pytest.fixture(scope="module")
def host():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    with draw.bound("topup-host"):
        return StarSystem(config).star


def _plan(shortfall_earth=20.0, level="Normal", stars=4):
    return rogue_topup.plan(draw.Stream("topup-plan"), level, stars, shortfall_earth * constants.EARTH_MASS_TO_KG, [1.0e30], sector_edge_pc=EDGE_PC)


def test_each_planned_rogue_is_a_planet_of_the_planned_mass_with_no_star(host):
    plan = _plan()
    assert plan.rogues
    made = rogue_topup_build.build(plan.rogues, [host], SystemConfig(), "sector-A", ADDRESS, EDGE_PC)
    assert len(made) == len(plan.rogues)
    for item in made:
        assert item.planet.mass == pytest.approx(item.rogue.mass_earth * constants.EARTH_MASS_TO_KG, rel=1e-6)
        assert item.planet.star.luminosity == 0.0
        assert "apparent" in item.note and "pc" in item.note


def test_a_rogue_lands_inside_the_sector_it_is_placed_in(host):
    made = rogue_topup_build.build(_plan().rogues, [host], SystemConfig(), "sector-A", ADDRESS, EDGE_PC)
    center = spherical.sector_position_pc(*ADDRESS, EDGE_PC)
    for item in made:
        assert math.dist(item.landing.position_pc, center) <= math.sqrt(3.0) * EDGE_PC / 2.0 + 1e-9
        assert item.landing.address == ADDRESS
        assert item.planet.ejecta_velocity_kms == tuple(item.landing.velocity_km_s)


def test_the_same_sector_and_plan_give_the_same_rogues(host):
    plan = _plan()
    first = rogue_topup_build.build(plan.rogues, [host], SystemConfig(), "sector-A", ADDRESS, EDGE_PC)
    again = rogue_topup_build.build(plan.rogues, [host], SystemConfig(), "sector-A", ADDRESS, EDGE_PC)
    other = rogue_topup_build.build(plan.rogues, [host], SystemConfig(), "sector-B", ADDRESS, EDGE_PC)
    assert [(i.planet.radius, i.planet.planet_class) for i in first] == [(i.planet.radius, i.planet.planet_class) for i in again]
    assert [i.planet.radius for i in first] != [i.planet.radius for i in other]


def test_a_top_up_with_no_host_is_refused():
    with pytest.raises(rogue_topup.RogueTopUpError):
        rogue_topup_build.build(_plan().rogues, [], SystemConfig(), "sector-A", ADDRESS, EDGE_PC)


# --- stored ---------------------------------------------------------------------------------------------------------------------------------------

DB_ADDRESS = (2609, 7, 32)


@pytest.mark.fresh_database
def test_a_saved_top_up_has_no_donor_and_is_replaced_not_doubled(mysql_config, host):
    from planetgen.db import store
    from planetgen.galaxy.sector import SpaceSector

    sector = SpaceSector("Topup Sector", edge_ly=13.0)
    x, y, z = spherical.sector_position_pc(*DB_ADDRESS, 4.0)
    position = {"center_x_pc": x, "center_y_pc": y, "center_z_pc": z, "galactic_radius_pc": math.sqrt(x * x + y * y + z * z),
                "shell_index": DB_ADDRESS[0], "band_index": DB_ADDRESS[1], "wedge_index": DB_ADDRESS[2]}
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position=position)
    made = rogue_topup_build.build(_plan().rogues, [host], SystemConfig(), "sector-DB", DB_ADDRESS, 4.0)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            assert store.save_topup_rogues(conn, sector_id, made) == len(made)
            assert store.save_topup_rogues(conn, sector_id, made) == len(made)       # saved again: replaced
        rows = conn.execute("SELECT uid, donor_uid, event_index, rogue_note, since_ejection_gy FROM rogue_planets WHERE sector_id = ?", (sector_id,)).fetchall()
        assert len(rows) == len(made) == store.count_rogues_in_sector(conn, sector_id)
        assert all(r["donor_uid"] is None and r["uid"] is not None and "apparent" in r["rogue_note"] for r in rows)
    finally:
        conn.close()
