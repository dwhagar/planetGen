# tests/test_query_db_sector_neighbors.py

"""
Coverage for `queryDb.sector_neighbors` (and `sector_detail`'s own
`neighbors` key) -- the Sector Map's (`html/lib/starmap.py`) neighboring-
sector indicators, one clickable marker per immediately surrounding
sector address, each tagged with whether a real `sectors` row already
exists there.

Sectors are inserted directly via `stellarObjects._db.insert_sector` with
an explicit `galaxy_position` (real `galaxyGeometry` addresses and
positions, not hand-picked placeholders -- see `test_galaxy_gen.py`'s
`test_sectors_table_rejects_duplicate_address` for the same pattern) rather than through the much
heavier `generate.py galaxy` path, since these tests only need real
`sectors` rows at known addresses, not fully generated systems.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable,
same as every other database-backed test in this suite.
"""

import pytest

import queryDb
from stellarObjects import _db
from stellarObjects.galaxyGeometry import galactic_radius_pc, neighbor_addresses, sector_position_pc
from stellarObjects.spaceSector import SpaceSector

EDGE_PC = 3.526  # ~DEFAULT_SECTOR_EDGE_LY (11.5 ly) converted to parsecs.
ADDRESS = (5, 2, 17)


def _insert_placed_sector(conn, name, address):
    position = sector_position_pc(*address, EDGE_PC)
    galaxy_position = {
        "center_x_pc": position[0], "center_y_pc": position[1], "center_z_pc": position[2],
        "galactic_radius_pc": galactic_radius_pc(position),
        "ring_index": address[0], "layer_index": address[1], "ring_slot_index": address[2],
    }
    return _db.insert_sector(conn, SpaceSector(name=name), galaxy_position=galaxy_position)


def _address(neighbor):
    return neighbor["ring_index"], neighbor["layer_index"], neighbor["ring_slot_index"]


def test_sector_neighbors_reports_an_existing_neighbor_as_such(mysql_config):
    lateral = (ADDRESS[0], ADDRESS[1], ADDRESS[2] + 1)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            sector_id = _insert_placed_sector(conn, "This Sector", ADDRESS)
            neighbor_id = _insert_placed_sector(conn, "Neighbor Sector", lateral)

        sector = conn.execute("SELECT * FROM sectors WHERE id = ?", (sector_id,)).fetchone()
        neighbors = queryDb.sector_neighbors(conn, sector)
    finally:
        conn.close()

    matches = [n for n in neighbors if _address(n) == lateral]
    assert len(matches) == 1
    assert matches[0]["exists"] is True
    assert matches[0]["sector_id"] == neighbor_id
    assert matches[0]["sector_name"] == "Neighbor Sector"


def test_sector_neighbors_reports_every_not_yet_generated_neighbor_with_its_address(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            sector_id = _insert_placed_sector(conn, "This Sector", ADDRESS)
        sector = conn.execute("SELECT * FROM sectors WHERE id = ?", (sector_id,)).fetchone()
        neighbors = queryDb.sector_neighbors(conn, sector)
    finally:
        conn.close()

    # Nothing else was ever inserted -- every reported neighbor must be
    # "not yet generated".
    assert neighbors
    for neighbor in neighbors:
        assert neighbor["exists"] is False
        assert neighbor["sector_id"] is None
        assert neighbor["sector_name"] is None
        assert neighbor["designation"]

    assert [_address(n) for n in neighbors] == neighbor_addresses(*ADDRESS)
    reported = {_address(n) for n in neighbors}
    assert (ADDRESS[0], ADDRESS[1] + 1, ADDRESS[2]) in reported
    assert (ADDRESS[0], ADDRESS[1] - 1, ADDRESS[2]) in reported
    assert any(ring == ADDRESS[0] + 1 for ring, _layer, _slot in reported)
    assert any(ring == ADDRESS[0] - 1 for ring, _layer, _slot in reported)


def test_sector_neighbors_direction_points_at_the_neighbor(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            sector_id = _insert_placed_sector(conn, "This Sector", ADDRESS)
        sector = conn.execute("SELECT * FROM sectors WHERE id = ?", (sector_id,)).fetchone()
        neighbors = queryDb.sector_neighbors(conn, sector)
    finally:
        conn.close()

    above = next(n for n in neighbors if _address(n) == (ADDRESS[0], ADDRESS[1] + 1, ADDRESS[2]))
    assert above["direction_pc"] == pytest.approx([0.0, 0.0, EDGE_PC], abs=1e-3)


def test_sector_neighbors_is_empty_for_a_sector_with_no_galaxy_placement(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            sector_id = _db.insert_sector(conn, SpaceSector(name="Unplaced Sector"))
        sector = conn.execute("SELECT * FROM sectors WHERE id = ?", (sector_id,)).fetchone()
        neighbors = queryDb.sector_neighbors(conn, sector)
    finally:
        conn.close()

    assert neighbors == []


def test_sector_detail_includes_neighbors_key(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            sector_id = _insert_placed_sector(conn, "This Sector", ADDRESS)
        detail = queryDb.sector_detail(conn, sector_id)
    finally:
        conn.close()

    assert len(detail["neighbors"]) > 0
    for neighbor in detail["neighbors"]:
        assert "direction_pc" in neighbor and len(neighbor["direction_pc"]) == 3
