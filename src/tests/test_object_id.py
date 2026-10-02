# tests/test_object_id.py

"""
GEN.64: the 76-bit position ID (`stellarObjects/objectId.py`) every
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
    assert len(name) == 19 and name == name.upper()
    fields = objectId.parse_id(name)
    assert fields["kind"] == "rogue-planet" and fields["unit"] == "pc" and fields["distance"] == 8000
    course = course_between((0.0, 0.0, 0.0), (8000.0, 20.0, 5.0))
    assert fields["bearing_deg"] == pytest.approx(course.bearing_deg, abs=360 / 2 ** 22)
    assert fields["mark_deg"] == pytest.approx(course.mark_deg, abs=360 / 2 ** 22)
    assert int(name, 16) >> 70 == objectId.KIND_CODES["rogue-planet"]
    assert fields["collision"] == 0 and int(name, 16) & 0xF == 0


@pytest.mark.parametrize("distance_pc, unit, value", [
    (0.0, "mpc", 0), (12.5, "mpc", 12500), (524.287, "mpc", 524287), (600.0, "cpc", 60000),
    (5242.87, "cpc", 524287), (8000.4, "pc", 8000), (600_000.0, "kpc", 600), (3e9, "Mpc", 3000),
    (5e14, "Gpc", 500000),
])
def test_distance_takes_the_smallest_unit_it_fits(distance_pc, unit, value):
    fields = objectId.parse_id(objectId.format_id(objectId.pack("nebula", (distance_pc, 0.0, 0.0))))
    assert (fields["unit"], fields["distance"]) == (unit, value)


def test_kind_sets_the_type_bits_only():
    position = (100.0, -40.0, 3.0)
    ids = {kind: objectId.pack(kind, position) for kind in objectId.KIND_CODES}
    low = {object_id & ((1 << 70) - 1) for object_id in ids.values()}
    assert len(low) == 1 and len(set(ids.values())) == len(objectId.KIND_CODES)


def test_unknown_kind_and_bad_names_are_refused():
    with pytest.raises(ValueError):
        objectId.pack("star-system", (1.0, 0.0, 0.0))
    good = objectId.format_id(objectId.pack("rogue-planet", (8000.0, 20.0, 5.0)))
    bad_unit = objectId.format_id((objectId.KIND_CODES["nebula"] << 70) | (7 << 67))
    for name in ("Vestara", good[:-1], "0" * 17, "G" + good[1:], bad_unit, "F" + good[1:]):
        assert not objectId.is_object_id(name)
    assert objectId.is_object_id(good.lower())


def test_bump_counts_sixteen_collisions_then_moves_one_mark_step():
    object_id = objectId.pack("rogue-planet", (8000.0, 20.0, 5.0))
    seen = [object_id]
    for _ in range(16):
        seen.append(objectId.bump(seen[-1]))
    assert [objectId.collision(i) for i in seen] == list(range(16)) + [0]
    assert seen[15] == object_id + 15
    fields, next_step = (objectId.parse_id(objectId.format_id(i)) for i in (object_id, seen[16]))
    assert next_step["mark_deg"] == pytest.approx(fields["mark_deg"] + 360 / 2 ** 22)
    assert ({k: v for k, v in next_step.items() if k != "mark_deg"}
            == {k: v for k, v in fields.items() if k != "mark_deg"})
    mark_and_collision = (((1 << 22) - 1) << 4) | 0xF
    assert objectId.bump(object_id | mark_and_collision) == object_id & ~mark_and_collision


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


def test_a_remnant_core_gets_its_own_core_id(mysql_config):
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
    assert objectId.parse_id(remnant.name)["kind"] == "supernova-remnant"
    core = objectId.parse_id(remnant.compact_remnant.name)
    assert core["kind"] in ("black-hole-core", "neutron-star-core")
