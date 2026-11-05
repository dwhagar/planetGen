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


@pytest.mark.fresh_database
def test_a_rogue_born_elsewhere_is_a_resident_but_one_born_here_is_not(mysql_config, host):
    from planetgen.db import store
    from planetgen.galaxy import object_uid
    from planetgen.galaxy.sector import SpaceSector

    sector = SpaceSector("Resident Sector", edge_ly=13.0)
    x, y, z = spherical.sector_position_pc(*DB_ADDRESS, 4.0)
    position = {"center_x_pc": x, "center_y_pc": y, "center_z_pc": z, "galactic_radius_pc": math.sqrt(x * x + y * y + z * z),
                "shell_index": DB_ADDRESS[0], "band_index": DB_ADDRESS[1], "wedge_index": DB_ADDRESS[2]}
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position=position)
    made = rogue_topup_build.build(_plan().rogues, [host], SystemConfig(), "sector-RES", DB_ADDRESS, 4.0)
    conn = store.get_connection(mysql_config)
    try:
        assert store.arrived_residents(conn, DB_ADDRESS) == []
        with conn:
            store.save_topup_rogues(conn, sector_id, made)
        assert store.arrived_residents(conn, DB_ADDRESS) == []                       # born here: this sector's own fill replaces it
        other = object_uid.pack(DB_ADDRESS[0], DB_ADDRESS[1], DB_ADDRESS[2] + 1, object_uid.SERIAL_GENERATED, 1)
        with conn:
            conn.execute("UPDATE rogue_planets SET uid = ? WHERE sector_id = ? LIMIT 1", (object_uid.to_bytes(other), sector_id))
        found = store.arrived_residents(conn, DB_ADDRESS)
        assert len(found) == 1 and found[0].kind == "rogue" and found[0].mass_kg > 0.0
    finally:
        conn.close()


@pytest.mark.fresh_database
def test_the_stage_tops_a_charted_sector_up_once(mysql_config, redis_server):
    from types import SimpleNamespace

    from planetgen.db import chart_ledger, store
    from planetgen.galaxy.sector import SpaceSector
    from planetgen.generation import topup_stage

    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    with draw.bound("stage-system"):
        system = StarSystem(config)
    sector = SpaceSector("Stage Sector", edge_ly=13.0)
    sector.add_system(system, position=(0.0, 0.0, 0.0))
    x, y, z = spherical.sector_position_pc(*DB_ADDRESS, 4.0)
    position = {"center_x_pc": x, "center_y_pc": y, "center_z_pc": z, "galactic_radius_pc": math.sqrt(x * x + y * y + z * z),
                "shell_index": DB_ADDRESS[0], "band_index": DB_ADDRESS[1], "wedge_index": DB_ADDRESS[2]}
    store.save_sector(sector, config=mysql_config, galaxy_position=position)
    # a profile that expects far more mass than the sector holds: the shortfall is large, so the level is what stops the top-up
    profile = SimpleNamespace(expected_systems=50.0, mix={"old": 1.0}, shape=_shape(), position_pc=(x, y, z))
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            chart_ledger.record(conn, DB_ADDRESS, None, 1)
        with conn:
            result = topup_stage.topup_sector(conn, DB_ADDRESS, profile, "few", SystemConfig(), 4.0)
        stars = conn.execute("SELECT COUNT(*) AS n FROM stars").fetchone()["n"]
        assert result.skipped is None and result.written == int(2.5 * stars) and result.stopped_by_level      # the Few level, 2.5 a star
        assert result.fields >= 1 and conn.execute("SELECT COUNT(*) AS n FROM asteroid_fields WHERE mass_kg > 0").fetchone()["n"] == result.fields      # the left-over mass is fields
        assert chart_ledger.row(conn, DB_ADDRESS)["topped_up_at"] is not None
        with conn:
            again = topup_stage.topup_sector(conn, DB_ADDRESS, profile, "few", SystemConfig(), 4.0)
        assert again.skipped == "already topped up" and again.written == 0
        assert store.count_rogues_in_sector(conn, store.kept_sector_id_at(conn, *DB_ADDRESS)) == result.written
    finally:
        conn.close()


def _shape():
    from planetgen.galaxy import shape_draw
    return shape_draw.default_shape(b"topup-stage-seed")
