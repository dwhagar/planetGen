# tests/test_query_db_sector_neighbors.py

"""
Coverage for `queryDb.sector_neighbors` (and `sector_detail`'s own
`neighbors` key) -- the Sector Map's (`html/lib/starmap.py`) neighboring-
sector indicators, one clickable marker per immediately surrounding
sector address, each tagged with whether a real `sectors` row already
exists there.

Sectors are inserted directly via `stellarObjects._db.insert_sector` with
an explicit `galaxy_position` (real `sectorGeometry`/`galaxyGeometry`
addresses and positions, not hand-picked placeholders -- see
`test_galaxy_gen.py`'s identical `test_sectors_table_rejects_duplicate_
shell_slot_address` for the same pattern) rather than through the much
heavier `generate.py galaxy` path, since these tests only need real
`sectors` rows at known addresses, not fully generated systems.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable,
same as every other database-backed test in this suite.
"""

import pytest

import queryDb
from stellarObjects import _db
from stellarObjects.galaxyGeometry import galactic_radius_pc, sector_position_pc
from stellarObjects.sectorGeometry import lateral_neighbor_slots, radial_neighbor_slot
from stellarObjects.spaceSector import SpaceSector

EDGE_PC = 3.526  # ~DEFAULT_SECTOR_EDGE_LY (11.5 ly) converted to parsecs.


def _insert_placed_sector(conn, name, shell_index, shell_slot_index):
    position = sector_position_pc(shell_index, shell_slot_index, EDGE_PC)
    galaxy_position = {
        "center_x_pc": position[0], "center_y_pc": position[1], "center_z_pc": position[2],
        "galactic_radius_pc": galactic_radius_pc(position),
        "shell_index": shell_index, "shell_slot_index": shell_slot_index,
        "vertices_pc": {"inner": [], "outer": []},
    }
    return _db.insert_sector(conn, SpaceSector(name=name), galaxy_position=galaxy_position)


def test_sector_neighbors_reports_an_existing_lateral_neighbor_as_such(mysql_config):
    shell_index, shell_slot_index = 5, 100
    lateral_slot = lateral_neighbor_slots(shell_index, shell_slot_index, EDGE_PC)[0]

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            sector_id = _insert_placed_sector(conn, "This Sector", shell_index, shell_slot_index)
            neighbor_id = _insert_placed_sector(conn, "Neighbor Sector", shell_index, lateral_slot)

        sector = conn.execute("SELECT * FROM sectors WHERE id = ?", (sector_id,)).fetchone()
        neighbors = queryDb.sector_neighbors(conn, sector)
    finally:
        conn.close()

    matches = [n for n in neighbors if n["shell_index"] == shell_index and n["shell_slot_index"] == lateral_slot]
    assert len(matches) == 1
    assert matches[0]["exists"] is True
    assert matches[0]["sector_id"] == neighbor_id
    assert matches[0]["sector_name"] == "Neighbor Sector"


def test_sector_neighbors_reports_a_not_yet_generated_neighbor_with_its_address(mysql_config):
    shell_index, shell_slot_index = 5, 100
    lateral_slots = lateral_neighbor_slots(shell_index, shell_slot_index, EDGE_PC)
    assert lateral_slots

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            sector_id = _insert_placed_sector(conn, "This Sector", shell_index, shell_slot_index)
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

    reported = {(n["shell_index"], n["shell_slot_index"]) for n in neighbors}
    for slot in lateral_slots:
        assert (shell_index, slot) in reported
    outward = radial_neighbor_slot(shell_index, shell_slot_index, EDGE_PC, 1)
    assert (shell_index + 1, outward) in reported


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
    shell_index, shell_slot_index = 5, 100
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            sector_id = _insert_placed_sector(conn, "This Sector", shell_index, shell_slot_index)
        detail = queryDb.sector_detail(conn, sector_id)
    finally:
        conn.close()

    assert len(detail["neighbors"]) > 0
    for neighbor in detail["neighbors"]:
        assert "direction_pc" in neighbor and len(neighbor["direction_pc"]) == 3
