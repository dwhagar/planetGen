# tests/test_phenomena_placement.py

"""
Tests for schema v18's galaxy-frame placement of nebulae/asteroid fields --
`stellarObjects._db.compute_phenomenon_placement`/`insert_nebula`/
`insert_asteroid_field`/`save_phenomenon`, and the read side
(`queryDb.phenomena_near_sector`/`galaxy_placed_phenomena`) `html/
lib/starmap.py`'s Sector Map and `html/lib/galaxymap.py`'s Galaxy Map
build on. See `schema.sql`'s "v18" header note for the full design: a
nebula/asteroid field gets a real galaxy-frame sphere (`center_x/y/z_pc`,
`radius_ly`) centered near a given sector, not a sector-relative offset,
since it's frequently far larger than any one sector's own cube.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
they're skipped, not failed, when no MySQL test server is configured/
reachable.
"""

import math

import pytest

import queryDb
from stellarObjects import _db
from stellarObjects.asteroidFieldData import AsteroidField
from stellarObjects.config import SystemConfig
from stellarObjects.nebulaData import Nebula
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.utils import ly_to_pc, mpc_to_pc, pc_to_ly


def _place_sector(mysql_config, name, center_pc, edge_ly=11.5):
    """Persists an empty, galaxy-placed sector at an exact `center_pc` --
    a plain synthetic placement (no real Fibonacci-sphere/shell address),
    which is all `compute_phenomenon_placement`/`phenomena_near_sector`
    need (they only ever read `sectors.center_x/y/z_pc`/`edge_mpc`, never
    `shell_index`/`shell_slot_index`)."""
    sector = SpaceSector(name, edge_ly=edge_ly)
    cx, cy, cz = center_pc
    galaxy_position = {
        "center_x_pc": cx, "center_y_pc": cy, "center_z_pc": cz,
        "galactic_radius_pc": math.sqrt(cx * cx + cy * cy + cz * cz),
        "vertices_pc": {"inner": [], "outer": []},
    }
    return _db.save_sector(sector, config=mysql_config, galaxy_position=galaxy_position)


def test_compute_phenomenon_placement_center_falls_within_sector_cube_half_extent(mysql_config):
    sector_id = _place_sector(mysql_config, "Placement Test Sector", (100.0, -50.0, 25.0), edge_ly=11.5)
    edge_pc = ly_to_pc(11.5)

    conn = _db.get_connection(mysql_config)
    try:
        placement = _db.compute_phenomenon_placement(conn, sector_id)
    finally:
        conn.close()

    for axis_value, sector_center in zip(
        (placement["center_x_pc"], placement["center_y_pc"], placement["center_z_pc"]),
        (100.0, -50.0, 25.0),
    ):
        assert abs(axis_value - sector_center) <= edge_pc / 2 + 1e-9

    expected_radius = math.sqrt(
        placement["center_x_pc"] ** 2 + placement["center_y_pc"] ** 2 + placement["center_z_pc"] ** 2
    )
    assert placement["galactic_radius_pc"] == pytest.approx(expected_radius)


def test_compute_phenomenon_placement_rejects_an_unplaced_sector(mysql_config):
    sector = SpaceSector("Unplaced Sector", edge_ly=11.5)
    sector_id = _db.save_sector(sector, config=mysql_config)  # no galaxy_position

    conn = _db.get_connection(mysql_config)
    try:
        with pytest.raises(ValueError, match="no galaxy placement"):
            _db.compute_phenomenon_placement(conn, sector_id)
    finally:
        conn.close()


def test_compute_phenomenon_placement_rejects_a_nonexistent_sector(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with pytest.raises(ValueError):
            _db.compute_phenomenon_placement(conn, 999999)
    finally:
        conn.close()


def test_save_phenomenon_persists_placement_and_nearest_sector_for_nebula_and_asteroid_field(mysql_config):
    sector_id = _place_sector(mysql_config, "Nebula Home Sector", (10.0, 20.0, 30.0))

    cfg = SystemConfig()
    nebula = Nebula(cfg)
    field = AsteroidField(cfg)

    nebula_id = _db.save_phenomenon(nebula, cfg, "nebula", config=mysql_config, sector_id=sector_id)
    field_id = _db.save_phenomenon(field, cfg, "asteroid-field", config=mysql_config, sector_id=sector_id)

    conn = _db.get_connection(mysql_config)
    try:
        nebula_row = conn.execute("SELECT * FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
        field_row = conn.execute("SELECT * FROM asteroid_fields WHERE id = ?", (field_id,)).fetchone()
    finally:
        conn.close()

    for row in (nebula_row, field_row):
        assert row["sector_id"] == sector_id
        assert row["center_x_pc"] is not None
        assert row["center_y_pc"] is not None
        assert row["center_z_pc"] is not None
        assert row["galactic_radius_pc"] == pytest.approx(
            math.sqrt(row["center_x_pc"] ** 2 + row["center_y_pc"] ** 2 + row["center_z_pc"] ** 2)
        )


def test_save_phenomenon_without_sector_id_leaves_it_unplaced(mysql_config):
    cfg = SystemConfig()
    nebula = Nebula(cfg)
    nebula_id = _db.save_phenomenon(nebula, cfg, "nebula", config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        row = conn.execute("SELECT * FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
    finally:
        conn.close()

    assert row["sector_id"] is None
    assert row["center_x_pc"] is None
    assert row["center_y_pc"] is None
    assert row["center_z_pc"] is None
    assert row["galactic_radius_pc"] is None


def test_phenomena_near_sector_finds_a_nebula_placed_at_that_sector(mysql_config):
    sector_id = _place_sector(mysql_config, "Nearby Sector", (200.0, 0.0, 0.0))

    cfg = SystemConfig()
    nebula = Nebula(cfg)
    nebula.radius_ly = 5.0  # small and deterministic, for a clean "definitely inside" case
    nebula_id = _db.save_phenomenon(nebula, cfg, "nebula", config=mysql_config, sector_id=sector_id)

    conn = _db.get_connection(mysql_config)
    try:
        matches = queryDb.phenomena_near_sector(conn, sector_id)
    finally:
        conn.close()

    assert len(matches) == 1
    match = matches[0]
    assert match["id"] == nebula_id
    assert match["type"] == "nebula"
    assert match["name"] == nebula.name
    assert match["descriptor"] == nebula.nebula_type
    assert match["radius_ly"] == pytest.approx(5.0)
    # offset_*_ly must actually be the phenomenon's real galaxy-frame
    # offset from the sector's own center, not a placeholder.
    expected_distance = math.sqrt(
        match["offset_x_ly"] ** 2 + match["offset_y_ly"] ** 2 + match["offset_z_ly"] ** 2
    )
    assert match["distance_ly"] == pytest.approx(expected_distance)


def test_phenomena_near_sector_excludes_a_small_far_nebula_but_includes_a_huge_one(mysql_config):
    # Two sectors far enough apart (500 pc ~= 1630 ly) that a small nebula
    # placed at sector A's own position can't plausibly reach sector B's
    # cube -- but a nebula whose own radius exceeds that separation must
    # still be found from B too (the whole point of a real sphere-vs-cube
    # overlap test rather than trusting the "nearest sector" link alone).
    sector_a = _place_sector(mysql_config, "Sector A", (0.0, 0.0, 0.0))
    sector_b = _place_sector(mysql_config, "Sector B", (500.0, 0.0, 0.0))

    cfg = SystemConfig()
    small_nebula = Nebula(cfg)
    small_nebula.radius_ly = 5.0
    small_id = _db.save_phenomenon(small_nebula, cfg, "nebula", config=mysql_config, sector_id=sector_a)

    huge_nebula = Nebula(cfg, name="Huge Nebula")
    huge_nebula.radius_ly = 3000.0  # bigger than the whole 500 pc (~1630 ly) separation
    huge_id = _db.save_phenomenon(huge_nebula, cfg, "nebula", config=mysql_config, sector_id=sector_a)

    conn = _db.get_connection(mysql_config)
    try:
        near_a = {m["id"] for m in queryDb.phenomena_near_sector(conn, sector_a)}
        near_b = {m["id"] for m in queryDb.phenomena_near_sector(conn, sector_b)}
    finally:
        conn.close()

    assert near_a == {small_id, huge_id}
    assert near_b == {huge_id}  # the small one is real but too far away


def test_phenomena_near_sector_is_empty_for_an_unplaced_sector(mysql_config):
    sector = SpaceSector("Unplaced Sector", edge_ly=11.5)
    sector_id = _db.save_sector(sector, config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        assert queryDb.phenomena_near_sector(conn, sector_id) == []
    finally:
        conn.close()


def test_galaxy_placed_phenomena_returns_only_placed_ones(mysql_config):
    sector_id = _place_sector(mysql_config, "Galaxy Dot Sector", (42.0, 42.0, 42.0))

    cfg = SystemConfig()
    placed_nebula = Nebula(cfg)
    placed_id = _db.save_phenomenon(placed_nebula, cfg, "nebula", config=mysql_config, sector_id=sector_id)
    unplaced_field = AsteroidField(cfg)
    _db.save_phenomenon(unplaced_field, cfg, "asteroid-field", config=mysql_config)  # no sector_id

    conn = _db.get_connection(mysql_config)
    try:
        placed = queryDb.galaxy_placed_phenomena(conn)
    finally:
        conn.close()

    ids = {(p["id"], p["type"]) for p in placed}
    assert (placed_id, "nebula") in ids
    assert len(placed) == 1  # the unplaced asteroid field must not appear

    entry = next(p for p in placed if p["id"] == placed_id and p["type"] == "nebula")
    assert entry["galactic_radius_pc"] == pytest.approx(
        math.sqrt(entry["x"] ** 2 + entry["y"] ** 2 + entry["z"] ** 2)
    )


def test_sector_detail_includes_nearby_phenomena(mysql_config):
    sector_id = _place_sector(mysql_config, "Detail Sector", (7.0, 8.0, 9.0))
    cfg = SystemConfig()
    field = AsteroidField(cfg)
    field_id = _db.save_phenomenon(field, cfg, "asteroid-field", config=mysql_config, sector_id=sector_id)

    conn = _db.get_connection(mysql_config)
    try:
        detail = queryDb.sector_detail(conn, sector_id)
    finally:
        conn.close()

    assert any(p["id"] == field_id and p["type"] == "asteroid_field" for p in detail["phenomena"])


def test_pc_ly_round_trip_used_by_placement_math():
    # Sanity-check the conversion helpers this module's own math leans on --
    # not new code, but a cheap guard against a future accidental change to
    # either direction silently breaking the placement/overlap math above.
    for value in (0.1, 1.0, 11.5, 500.0):
        assert pc_to_ly(ly_to_pc(value)) == pytest.approx(value)
        assert mpc_to_pc(value * 1000) == pytest.approx(value)
