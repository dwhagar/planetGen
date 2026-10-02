# tests/test_object_id.py

"""
GEN.64: the 64-bit position ID (`stellarObjects/objectId.py`) every
interstellar object and bright-sweep system is named by, and how a sector
save claims those IDs (`_db._claim_object_ids`, `insert_sector`).

The database tests take the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import pytest

from stellarObjects import _db, objectId
from stellarObjects.config import SystemConfig
from stellarObjects.navigation import course_between
from stellarObjects.roguePlanetData import RoguePlanet
from stellarObjects.spaceSector import SpaceSector

_POSITION = {"center_x_pc": 8000.0, "center_y_pc": 20.0, "center_z_pc": 5.0, "galactic_radius_pc": 8000.03}


def test_fields_pack_in_their_bit_ranges():
    name = objectId.format_id(objectId.pack("rogue-planet", (8000.0, 20.0, 5.0)))
    assert len(name) == 16 and name == name.upper()
    fields = objectId.parse_id(name)
    assert fields["kind"] == "rogue-planet" and fields["unit"] == "pc" and fields["distance"] == 8000
    course = course_between((0.0, 0.0, 0.0), (8000.0, 20.0, 5.0))
    assert fields["bearing_deg"] == pytest.approx(course.bearing_deg, abs=360 / 2 ** 20)
    assert fields["mark_deg"] == pytest.approx(course.mark_deg, abs=360 / 2 ** 20)
    assert int(name, 16) >> 60 == objectId.KIND_CODES["rogue-planet"]


@pytest.mark.parametrize("distance_pc, unit, value", [
    (0.0, "mpc", 0), (12.5, "mpc", 12500), (131.071, "mpc", 131071), (500.0, "cpc", 50000),
    (1310.71, "cpc", 131071), (8000.4, "pc", 8000), (250_000.0, "kpc", 250), (3e9, "Mpc", 3000),
    (1.3e14, "Gpc", 130000),
])
def test_distance_takes_the_smallest_unit_it_fits(distance_pc, unit, value):
    fields = objectId.parse_id(objectId.format_id(objectId.pack("nebula", (distance_pc, 0.0, 0.0))))
    assert (fields["unit"], fields["distance"]) == (unit, value)


def test_kind_sets_the_type_bits_only():
    position = (100.0, -40.0, 3.0)
    ids = {kind: objectId.pack(kind, position) for kind in objectId.KIND_CODES}
    low = {object_id & ((1 << 60) - 1) for object_id in ids.values()}
    assert len(low) == 1 and len(set(ids.values())) == len(objectId.KIND_CODES)


def test_unknown_kind_and_bad_names_are_refused():
    with pytest.raises(ValueError):
        objectId.pack("star-system", (1.0, 0.0, 0.0))
    for name in ("Vestara", "140AF0FEC9CFFFB", "0000000000000000", "G40AF0FEC9CFFFB8", "1E0AF0FEC9CFFFB8"):
        assert not objectId.is_object_id(name)
    assert objectId.is_object_id("140af0fec9cffFB8")


def test_bump_moves_one_mark_step_and_wraps_inside_the_mark_field():
    object_id = objectId.pack("rogue-planet", (8000.0, 20.0, 5.0))
    assert objectId.bump(object_id) == object_id + 1 or object_id & 0xFFFFF == 0xFFFFF
    top = (object_id | 0xFFFFF)
    assert objectId.bump(top) == top & ~0xFFFFF


def test_same_point_ids_are_bumped_in_generation_order(mysql_config):
    first, second = RoguePlanet(SystemConfig()), RoguePlanet(SystemConfig())
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db._claim_object_ids(conn, [(first, "rogue-planet", (10.0, 0.0, 0.0)),
                                         (second, "rogue-planet", (10.0, 0.0, 0.0))])
    finally:
        conn.close()
    plain = objectId.pack("rogue-planet", (10.0, 0.0, 0.0))
    assert first.name == objectId.format_id(plain)
    assert second.name == objectId.format_id(objectId.bump(plain))


def test_a_stored_id_is_skipped(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            stored = RoguePlanet(SystemConfig())
            _db.insert_rogue_planet(conn, stored, placement={"center_x_pc": 10.0, "center_y_pc": 0.0,
                                                             "center_z_pc": 0.0, "galactic_radius_pc": 10.0})
            later = RoguePlanet(SystemConfig())
            _db.insert_rogue_planet(conn, later, placement={"center_x_pc": 10.0, "center_y_pc": 0.0,
                                                            "center_z_pc": 0.0, "galactic_radius_pc": 10.0})
            registry = conn.execute("SELECT COUNT(*) AS n FROM system_name_registry").fetchone()["n"]
    finally:
        conn.close()
    plain = objectId.pack("rogue-planet", (10.0, 0.0, 0.0))
    assert (stored.name, later.name) == (objectId.format_id(plain), objectId.format_id(objectId.bump(plain)))
    assert registry == 0


def test_unplaced_objects_and_given_names_keep_their_names(mysql_config):
    unplaced, given = RoguePlanet(SystemConfig()), RoguePlanet(SystemConfig(), name="Drifter")
    generated = unplaced.name
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_rogue_planet(conn, unplaced)
            _db.insert_rogue_planet(conn, given, placement={"center_x_pc": 10.0, "center_y_pc": 0.0,
                                                            "center_z_pc": 0.0, "galactic_radius_pc": 10.0})
    finally:
        conn.close()
    assert (unplaced.name, given.name) == (generated, "Drifter")


def test_a_placed_sector_names_its_phenomena_by_id(mysql_config):
    sector = SpaceSector("Halfway Sector", edge_ly=11.5)
    for index in range(3):
        sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet",
                              position=(1.0, float(index) + 2.5, 1.0))
    _db.save_sector(sector, config=mysql_config, galaxy_position=_POSITION)
    names = [entry.phenomenon.name for entry in sector.phenomena]
    expected = [objectId.format_id(objectId.pack(
        "rogue-planet", _db._placement_center(_db._galaxy_placement_from_sector_offset(_POSITION, entry.position))))
        for entry in sector.phenomena]
    assert names == expected
    conn = _db.get_connection(mysql_config)
    try:
        stored = sorted(row["name"] for row in conn.execute("SELECT name FROM rogue_planets").fetchall())
        registry = conn.execute("SELECT COUNT(*) AS n FROM system_name_registry").fetchone()["n"]
    finally:
        conn.close()
    assert stored == sorted(names) and registry == 0


def test_a_remnant_core_follows_its_remnant_id(mysql_config):
    from stellarObjects.supernovaRemnantData import SupernovaRemnant
    for _ in range(200):
        remnant = SupernovaRemnant(SystemConfig())
        if remnant.compact_remnant is not None:
            break
    else:
        pytest.skip("no remnant with a core drawn")
    assert remnant.compact_remnant.name == f"{remnant.name} Core"
    sector = SpaceSector("Halfway Sector", edge_ly=11.5)
    sector.add_phenomenon(remnant, "supernova-remnant", position=(1.0, 2.5, 1.0))
    _db.save_sector(sector, config=mysql_config, galaxy_position=_POSITION)
    assert objectId.is_object_id(remnant.name)
    assert remnant.compact_remnant.name == f"{remnant.name} Core"
