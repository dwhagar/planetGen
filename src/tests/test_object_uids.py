# tests/test_object_uids.py

"""
Unique IDs for every object (GEN.69, `galaxy/uid.py`, `store.assign_uids`):
the pure functions, then a saved sector's rows checked against them.
"""

import pytest

from planetgen.db import store
from planetgen.galaxy import uid as galaxy_uid
from planetgen.galaxy.geometry import provisional_sector_designation
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.system import StarSystem
from planetgen.names import object_id

SEED = bytes(range(16))


def test_a_sector_id_is_its_designation_and_needs_no_row():
    uid = galaxy_uid.sector_uid(12, 3, 7)
    assert galaxy_uid.format_sector_uid(uid) == provisional_sector_designation(12, 3, 7)
    assert galaxy_uid.sector_address(uid) == (12, 3, 7)


@pytest.mark.parametrize("bits", [galaxy_uid.GALAXY_BITS, galaxy_uid.LOCAL_BITS])
def test_a_derived_id_is_fixed_by_seed_kind_parent_and_index(bits):
    uid = galaxy_uid.derived_uid(SEED, "planet", "ABC", 2, bits)
    assert uid == galaxy_uid.derived_uid(SEED, "planet", "ABC", 2, bits)
    assert 0 <= uid < (1 << bits)
    for other in (galaxy_uid.derived_uid(bytes(16), "planet", "ABC", 2, bits),
                  galaxy_uid.derived_uid(SEED, "moon", "ABC", 2, bits),
                  galaxy_uid.derived_uid(SEED, "planet", "ABD", 2, bits),
                  galaxy_uid.derived_uid(SEED, "planet", "ABC", 3, bits)):
        assert other != uid


def test_a_galaxy_wide_id_can_never_equal_a_position_id():
    assert all(galaxy_uid.derived_uid(SEED, "system", "S", n) >> 95 == 1 for n in range(50))
    assert (1 << object_id.ID_BITS) <= (1 << 95), "a position ID has its top 20 bits clear"


def test_no_seed_draws_from_the_zero_seed():
    assert galaxy_uid.derived_uid(None, "system", "S", 0) == galaxy_uid.derived_uid(galaxy_uid.NO_SEED, "system", "S", 0)


def test_ids_round_trip_through_bytes_and_hex():
    uid = galaxy_uid.derived_uid(SEED, "system", "S", 1)
    assert galaxy_uid.uid_from_bytes(galaxy_uid.uid_bytes(uid)) == uid
    assert len(galaxy_uid.format_uid(uid)) == 24
    local = galaxy_uid.derived_uid(SEED, "star", "S", 0, galaxy_uid.LOCAL_BITS)
    assert len(galaxy_uid.format_uid(local, galaxy_uid.LOCAL_BITS)) == 16


def _sector(name, count, rogue):
    sector = SpaceSector(name, edge_ly=11.5)
    for index in range(count):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        system.name = f"Uidstar {name} {index}"
        sector.add_system(system, position=(float(index), 0.0, 0.0), system_config=cfg)
    sector.add_phenomenon(RoguePlanet(SystemConfig(), name=rogue), "rogue-planet", position=(0.0, 0.5, 0.0))
    return sector


PLACE = {"center_x_pc": 1000.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 1000.0,
         "ring_index": 250, "layer_index": 0, "ring_slot_index": 5}


def _rows(config, sql):
    conn = store.get_connection(config)
    try:
        return conn.execute(sql).fetchall()
    finally:
        conn.close()


def test_a_saved_sector_gives_every_object_its_id(mysql_config):
    store.save_sector(_sector("Uid One", 3, "Uid Wanderer"), config=mysql_config, galaxy_position=PLACE)

    [sector] = _rows(mysql_config, "SELECT uid FROM sectors")
    assert sector["uid"] == galaxy_uid.sector_uid(250, 0, 5)

    parent = galaxy_uid.format_sector_uid(sector["uid"])
    systems = _rows(mysql_config, "SELECT id, uid FROM star_systems ORDER BY id")
    assert [galaxy_uid.uid_from_bytes(row["uid"]) for row in systems] == [
        galaxy_uid.derived_uid(None, "system", parent, rank) for rank in range(3)]

    texts = {row["id"]: galaxy_uid.format_uid(galaxy_uid.uid_from_bytes(row["uid"])) for row in systems}
    for table, kind in (("stars", "star"), ("planets", "planet"), ("asteroid_belts", "belt")):
        for row in _rows(mysql_config, f"SELECT star_system_id, uid FROM {table}"):
            assert row["uid"] is not None, table
    planets = _rows(mysql_config, "SELECT id, star_system_id, uid FROM planets ORDER BY id")
    assert planets
    seen = {}
    for row in planets:
        rank = seen[row["star_system_id"]] = seen.get(row["star_system_id"], -1) + 1
        assert row["uid"] == galaxy_uid.derived_uid(None, "planet", texts[row["star_system_id"]], rank,
                                                    galaxy_uid.LOCAL_BITS)
    assert all(row["uid"] is not None for row in _rows(mysql_config, "SELECT uid FROM moons"))

    # A phenomenon named by hand takes a derived ID.
    [rogue] = _rows(mysql_config, "SELECT name, uid FROM rogue_planets")
    assert galaxy_uid.uid_from_bytes(rogue["uid"]) == galaxy_uid.derived_uid(None, "rogue_planets", parent, 0)


def test_an_interstellar_object_keeps_the_position_id_it_is_named_by(mysql_config):
    sector = _sector("Uid Pos", 1, "Placeholder")
    sector.phenomena[0].phenomenon.name_given = False
    store.save_sector(sector, config=mysql_config, galaxy_position=PLACE)
    [rogue] = _rows(mysql_config, "SELECT name, uid FROM rogue_planets")
    assert object_id.is_object_id(rogue["name"])
    assert galaxy_uid.uid_from_bytes(rogue["uid"]) == int(rogue["name"], 16)


def test_a_bright_sweep_system_is_named_by_the_registry_and_keeps_its_position_id(mysql_config):
    """GEN.72: once its sector is generated the system gets a word-salad name
    and keeps the position ID it was known by as its unique ID."""
    uids = []
    for slot in (5, 6):
        sector = _sector(f"Uid Bright {slot}", 1, "Placeholder")
        sector.entries[0].bright_star_id = 900 + slot
        # Both sectors put their system at the same galactic point.
        place = dict(PLACE, ring_slot_index=slot)
        expected = object_id.pack("bright-star", store._placement_center(
            store._galaxy_placement_from_sector_offset(place, sector.entries[0].position)))
        store.save_sector(sector, config=mysql_config, galaxy_position=place)
        uids.append((expected, sector.entries[0].star_system.name))
    rows = _rows(mysql_config, "SELECT name, uid FROM star_systems ORDER BY id")
    assert [row["name"] for row in rows] == [name for _expected, name in uids]
    assert not any(object_id.is_object_id(row["name"]) for row in rows)
    first = galaxy_uid.uid_from_bytes(rows[0]["uid"])
    assert first == uids[0][0]
    # The same point twice: the second is bumped past the stored ID.
    second = galaxy_uid.uid_from_bytes(rows[1]["uid"])
    assert uids[0][0] == uids[1][0] and second == object_id.bump(first)


def test_saving_the_same_sector_again_gives_the_same_ids(mysql_config):
    from tests.db_schema_support import scratch_database

    ids = []
    for config in (mysql_config, None):
        if config is None:
            with scratch_database(mysql_config) as other:
                store.save_sector(_sector("Uid Two", 2, "Uid Wanderer"), config=other, galaxy_position=PLACE)
                ids.append(_rows(other, "SELECT uid FROM sectors") and
                           [galaxy_uid.uid_from_bytes(r["uid"]) for r in _rows(other, "SELECT uid FROM star_systems ORDER BY id")])
        else:
            store.save_sector(_sector("Uid Two", 2, "Uid Wanderer"), config=config, galaxy_position=PLACE)
            ids.append([galaxy_uid.uid_from_bytes(r["uid"]) for r in _rows(config, "SELECT uid FROM star_systems ORDER BY id")])
    assert ids[0] == ids[1] and len(ids[0]) == 2


def test_a_system_saved_on_its_own_gets_ids(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    system.name = "Uidsolo One"
    system_id = store.save_system(system, cfg, config=mysql_config)
    [row] = _rows(mysql_config, f"SELECT uid FROM star_systems WHERE id = {system_id}")
    assert row["uid"] is not None
    assert all(r["uid"] is not None for r in _rows(mysql_config, f"SELECT uid FROM stars WHERE star_system_id = {system_id}"))
