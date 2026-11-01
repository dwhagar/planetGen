"""
GEN.1127: the chart ledger (`db/chart_ledger.py`), the sector row that cannot vanish under its residents, and the delete that takes only
what was born in a sector. A saved sector charts its cell; deleting it uncharts the cell; an object that drifted into a sector from another
(`sector_id` is where it is filed now, not where it was born) is not deleted with the sector it is filed in, and that sector is not deleted
under it; a sector emptied to be made again keeps its row when such a resident points at it, and the fill reuses the row.
"""

import math

import pytest

from planetgen.db import chart_ledger, edits, store
from planetgen.galaxy import object_uid, spherical
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.sector import SpaceSector
from planetgen import tuning

pytestmark = pytest.mark.db


def _addresses():
    low, high = spherical.band_range(40)
    first = (40, (low + high) // 2, 0)
    return first, sorted(spherical.neighbor_addresses(*first))[0]


def _sector(address, rogue):
    sector = SpaceSector(f"Ledger {address}", edge_ly=20.0)
    for _ in range(rogue):
        sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    return sector


def _save(config, address, rogue=1):
    centre = spherical.sector_position_pc(*address, float(tuning.DEFAULT_SECTOR_EDGE_PC))
    return store.save_sector(_sector(address, rogue), config=config, galaxy_position={
        "center_x_pc": float(centre[0]), "center_y_pc": float(centre[1]), "center_z_pc": float(centre[2]),
        "galactic_radius_pc": math.dist(centre, (0.0, 0.0, 0.0)),
        "shell_index": address[0], "band_index": address[1], "wedge_index": address[2]})


def _sector_id(conn, address):
    row = conn.execute("SELECT id FROM sectors WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", address).fetchone()
    return row["id"] if row else None


def _rogues(conn):
    return conn.execute("SELECT id, sector_id, uid FROM rogue_planets ORDER BY id").fetchall()


def test_saving_a_sector_charts_its_cell_and_deleting_it_uncharts_the_cell(mysql_config):
    first, second = _addresses()
    conn = store.get_connection(mysql_config)
    try:
        assert chart_ledger.state(conn, first, 1) == chart_ledger.UNCHARTED
        assert store.get_sector_id_at(conn, *first) is None
    finally:
        conn.close()
    _save(mysql_config, first)
    conn = store.get_connection(mysql_config)
    try:
        assert chart_ledger.state(conn, first, 1) == chart_ledger.CHARTED
        assert chart_ledger.state(conn, second, 1) == chart_ledger.UNCHARTED
        assert first in store.filled_sector_addresses(conn)
        assert store.get_sector_id_at(conn, *first) == _sector_id(conn, first)
        record = chart_ledger.row(conn, first)
        assert record["charted_at"] is not None and record["topped_up_at"] is None
        assert chart_ledger.mark_topped_up(conn, first) is True and chart_ledger.mark_topped_up(conn, first) is False
        edits.delete_sector_with_contents(conn, _sector_id(conn, first))
        conn.commit()
        assert chart_ledger.state(conn, first, 1) == chart_ledger.UNCHARTED
        assert first not in store.filled_sector_addresses(conn) and store.get_sector_id_at(conn, *first) is None
        assert _rogues(conn) == []
    finally:
        conn.close()


def test_a_resident_from_another_sector_is_not_deleted_with_the_sector_it_is_filed_in(mysql_config):
    first, second = _addresses()
    _save(mysql_config, first)
    _save(mysql_config, second)
    conn = store.get_connection(mysql_config)
    try:
        here, there = _sector_id(conn, first), _sector_id(conn, second)
        born_there = [r for r in _rogues(conn) if r["sector_id"] == there][0]
        assert bytes(born_there["uid"])[:7] == object_uid.sector_bytes(*second)
        conn.execute("UPDATE rogue_planets SET sector_id = ? WHERE id = ?", (here, born_there["id"]))    # it drifted in
        conn.commit()
        assert edits.sector_residents(conn, here) == {"rogue_planets": 1}
        with pytest.raises(edits.SectorHasResidents):
            edits.delete_sector_with_contents(conn, here)
        conn.rollback()
        assert len(_rogues(conn)) == 2
        # The sector it drifted in from is deleted with what was born there, wherever it is filed now: the resident goes, the host stays.
        assert edits.delete_sector_with_contents(conn, there) is not None
        conn.commit()
        assert [r["sector_id"] for r in _rogues(conn)] == [here]
        assert _sector_id(conn, first) == here
    finally:
        conn.close()


def test_a_sector_row_cannot_be_deleted_under_a_resident(mysql_config):
    first, second = _addresses()
    _save(mysql_config, first, rogue=0)
    _save(mysql_config, second)
    conn = store.get_connection(mysql_config)
    try:
        here = _sector_id(conn, first)
        conn.execute("UPDATE rogue_planets SET sector_id = ?", (here,))
        conn.commit()
        assert edits.restricting_rows(conn, here) == {"rogue_planets": 1}
        with pytest.raises(Exception) as caught:
            conn.execute("DELETE FROM sectors WHERE id = ?", (here,))
        assert "foreign key" in str(caught.value).lower()
        conn.rollback()
    finally:
        conn.close()


def test_a_sector_emptied_with_a_resident_keeps_its_row_and_the_fill_reuses_it(mysql_config):
    first, second = _addresses()
    _save(mysql_config, first)
    _save(mysql_config, second)
    conn = store.get_connection(mysql_config)
    try:
        here = _sector_id(conn, first)
        guest = [r for r in _rogues(conn) if r["sector_id"] == _sector_id(conn, second)][0]
        conn.execute("UPDATE rogue_planets SET sector_id = ? WHERE id = ?", (here, guest["id"]))
        conn.commit()
        counts = edits.empty_sector(conn, here)
        conn.commit()
        assert counts["kept"] is True
        assert _sector_id(conn, first) == here
        assert [r["id"] for r in _rogues(conn)] == [guest["id"]]       # the one born here went, the guest stays
        assert chart_ledger.state(conn, first, 1) == chart_ledger.UNCHARTED
        assert store.get_sector_id_at(conn, *first) is None            # the cell is open to be generated again
    finally:
        conn.close()
    _save(mysql_config, first)
    conn = store.get_connection(mysql_config)
    try:
        assert _sector_id(conn, first) == here and store.get_sector_id_at(conn, *first) == here
        assert chart_ledger.state(conn, first, 1) == chart_ledger.CHARTED
        assert guest["id"] in [r["id"] for r in _rogues(conn)]
        assert {r["sector_id"] for r in _rogues(conn)} == {here}
    finally:
        conn.close()


def test_a_sector_with_no_resident_is_deleted_by_emptying_it(mysql_config):
    first, _second = _addresses()
    _save(mysql_config, first)
    conn = store.get_connection(mysql_config)
    try:
        counts = edits.empty_sector(conn, _sector_id(conn, first))
        conn.commit()
        assert counts["kept"] is False and _sector_id(conn, first) is None
    finally:
        conn.close()
