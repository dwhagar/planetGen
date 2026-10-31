# tests/test_designation_history.py

"""
GEN.1138: the `designations` table and its writers. A rename (a system, a star, a planet or moon, a phenomenon) and the collision bump of
`_rank_chain` write the retired name in the transaction of the change; an ID changes only through `designation_history.renumber`, which
reserves the old one for ever and refuses a new one that is or ever was any object's.
"""

import pytest

from planetgen.db import designation_history, store
from planetgen.galaxy import designations, object_uid
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.system import StarSystem
from planetgen.names import object_id as objectId
from planetgen.util import draw

from tests.test_object_id import _POSITION

pytestmark = pytest.mark.db


def _planetary_system():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(p.moons for p in system.planets):
            return system, cfg
    pytest.fail("no system with a moon in 60 seeds")


PLACE = {"center_x_pc": 1000.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 1000.0,
         "shell_index": 250, "band_index": 0, "wedge_index": 5}


@pytest.fixture
def stored(mysql_config):
    """`(conn, system_id)` of one saved single-star system with planets and moons, in a placed sector (so every object has an ID)."""
    system, cfg = _planetary_system()
    sector = SpaceSector("Designation Sector", edge_ly=11.5)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config, galaxy_position=PLACE)
    conn = store.get_connection(mysql_config)
    try:
        yield conn, conn.execute("SELECT id FROM star_systems").fetchone()["id"]
    finally:
        conn.rollback()
        conn.close()


def _names(conn, table, row_id):
    return [(row.designation, row.reason) for row in designation_history.history(conn, table, row_id)]


def test_a_system_rename_keeps_the_old_name_and_its_bodies_follow(stored):
    conn, system_id = stored
    old = conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"]
    planets = {row["id"]: row["name"] for row in conn.execute("SELECT id, name FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()}
    assert any(name.startswith(old) for name in planets.values())
    with conn:
        assert store.rename_star_system(conn, system_id, "Zzyzx Prime")
    assert _names(conn, "star_systems", system_id) == [(old, "rename")]
    followed = [row_id for row_id, name in planets.items() if name.startswith(old)]
    assert followed
    for row_id in followed:
        retired = designation_history.history(conn, "planets", row_id)
        assert [(r.designation, r.reason) for r in retired] == [(planets[row_id], "parent change")]
    stars = conn.execute("SELECT id FROM stars WHERE star_system_id = ?", (system_id,)).fetchall()
    assert all(designation_history.history(conn, "stars", row["id"])[0].reason in ("rename", "parent change") for row in stars)


def test_changes_are_kept_in_order_and_the_first_one_has_no_start(stored):
    conn, system_id = stored
    first = conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"]
    with conn:
        store.rename_star_system(conn, system_id, "Zzyzx Prime")
        store.rename_star_system(conn, system_id, "Zzyzx Major", reason="collision bump")
    rows = designation_history.history(conn, "star_systems", system_id)
    assert [(r.designation, r.reason) for r in rows] == [("Zzyzx Prime", "collision bump"), (first, "rename")]
    assert rows[1].valid_from_epoch is None and rows[0].valid_from_epoch == rows[1].retired_epoch
    assert designations.former_lines(rows)[0][1] == "Zzyzx Prime"


def test_a_rename_to_the_same_name_writes_nothing(stored):
    conn, system_id = stored
    name = conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"]
    with conn:
        store.rename_star_system(conn, system_id, name)
    assert designation_history.history(conn, "star_systems", system_id) == []


def test_a_planet_rename_keeps_its_name_and_leaves_its_moons(stored):
    conn, system_id = stored
    planet = conn.execute("SELECT id, name FROM planets WHERE star_system_id = ? ORDER BY id LIMIT 1", (system_id,)).fetchone()
    with conn:
        assert store.rename_body(conn, "planets", planet["id"], "Hand-named")
    assert _names(conn, "planets", planet["id"]) == [(planet["name"], "rename")]
    assert conn.execute("SELECT COUNT(*) AS n FROM designations WHERE object_kind = 'moons'").fetchone()["n"] == 0


def test_a_reason_off_the_list_is_refused(stored):
    conn, system_id = stored
    with pytest.raises(designations.DesignationError):
        designation_history.retire_name(conn, "star_systems", system_id, "A", "B", "whim")


def test_the_collision_bump_keeps_the_name_it_moved_from(mysql_config):
    for wedge in (2, 1, 0):
        sector = SpaceSector(f"Sector {wedge}", edge_ly=11.5)
        sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet", position=(1.0, 2.5, 1.0))
        store.save_sector(sector, config=mysql_config, galaxy_position={**_POSITION, "shell_index": 1, "band_index": 0, "wedge_index": wedge})
    conn = store.get_connection(mysql_config)
    try:
        current = {row["id"]: row["name"] for row in conn.execute("SELECT id, name FROM rogue_planets").fetchall()}
        moved = conn.execute("SELECT object_row_id, name, reason FROM designations WHERE object_kind = 'rogue_planets'"
                             " ORDER BY id").fetchall()
        assert len(current) == 3 and moved and all(row["reason"] == "collision bump" for row in moved)
        for row in moved:
            assert row["name"] != current[row["object_row_id"]]
            assert objectId.parse_id(row["name"])["kind"] == "rogue-planet"
    finally:
        conn.close()


def _some_moon(conn, system_id):
    return conn.execute("SELECT id, uid FROM moons WHERE star_system_id = ? ORDER BY id LIMIT 1", (system_id,)).fetchone()


def _fresh_uid(conn, near):
    """An ID no object has held, in the same system as `near` (an integer)."""
    for body in range(object_uid.max_body(), 0, -1):
        candidate = object_uid.with_body(near, body)
        if not designation_history.is_live_uid(conn, candidate) and not designation_history.is_retired_uid(conn, candidate):
            return candidate
    pytest.fail("no free ID")


def test_a_renumbered_object_reserves_its_old_id_for_ever(stored):
    conn, system_id = stored
    moon = _some_moon(conn, system_id)
    old = object_uid.from_bytes(moon["uid"])
    new = _fresh_uid(conn, old)
    with conn:
        designation_history.renumber(conn, "moons", moon["id"], new)
    assert object_uid.from_bytes(conn.execute("SELECT uid FROM moons WHERE id = ?", (moon["id"],)).fetchone()["uid"]) == new
    assert [(r.designation, r.reason) for r in designation_history.history(conn, "moons", moon["id"])] == [(old, "admin renumber")]
    assert designation_history.is_retired_uid(conn, old) and not designation_history.is_live_uid(conn, old)
    with pytest.raises(designations.DesignationError, match="reserved for ever"):
        designation_history.renumber(conn, "moons", moon["id"], old)
    other = conn.execute("SELECT id FROM moons WHERE star_system_id = ? AND id <> ? LIMIT 1", (system_id, moon["id"])).fetchone()
    if other is not None:
        with pytest.raises(designations.DesignationError, match="reserved for ever"):
            designation_history.renumber(conn, "moons", other["id"], old)


def test_a_new_id_equal_to_a_live_one_is_refused(stored):
    conn, system_id = stored
    moon = _some_moon(conn, system_id)
    planet_uid = object_uid.from_bytes(conn.execute("SELECT uid FROM planets WHERE star_system_id = ? LIMIT 1", (system_id,)).fetchone()["uid"])
    with pytest.raises(designations.DesignationError, match="already an object's ID"):
        designation_history.renumber(conn, "moons", moon["id"], planet_uid)
    with pytest.raises(designations.DesignationError, match="already an object's ID"):
        designation_history.renumber(conn, "moons", moon["id"], object_uid.from_bytes(moon["uid"]))
    assert designation_history.history(conn, "moons", moon["id"]) == []


def test_an_id_others_embed_does_not_renumber_alone(stored):
    conn, system_id = stored
    row = conn.execute("SELECT id, uid FROM star_systems WHERE id = ?", (system_id,)).fetchone()
    with pytest.raises(designations.DesignationError, match="embedded"):
        designation_history.renumber(conn, "star_systems", row["id"], _fresh_uid(conn, object_uid.from_bytes(row["uid"])))


def test_a_reverted_id_is_a_new_one(stored):
    conn, system_id = stored
    moon = _some_moon(conn, system_id)
    first = object_uid.from_bytes(moon["uid"])
    second = _fresh_uid(conn, first)
    with conn:
        designation_history.renumber(conn, "moons", moon["id"], second)
    with pytest.raises(designations.DesignationError):
        designation_history.renumber(conn, "moons", moon["id"], first)
    third = _fresh_uid(conn, first)
    with conn:
        designation_history.renumber(conn, "moons", moon["id"], third, reason="regeneration")
    rows = designation_history.history(conn, "moons", moon["id"])
    assert [(r.designation, r.reason) for r in rows] == [(second, "regeneration"), (first, "admin renumber")]
    assert rows[0].valid_from_epoch == rows[1].retired_epoch


def test_the_table_refuses_a_row_that_is_neither_an_id_nor_a_name(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with pytest.raises(Exception):
            conn.execute("INSERT INTO designations (designation_kind, object_kind, object_row_id, retired_epoch, reason)"
                         " VALUES ('name', 'moons', 1, 1, 'rename')")
    finally:
        conn.rollback()
        conn.close()
