"""
GEN.202: one read gives every nebula kind of a box with no object twice. The stored nebulae, the field clouds no sector has stored and the
scattered planetary nebulae and remnants no sector has built are listed together, and an object that two of them hold (the same Object ID) is
listed once, the stored one first, whatever flag the other row carries.
"""

import pytest

from planetgen import tuning
from planetgen.db import field_clouds, nebula_index, query, store
from planetgen.galaxy.spherical import sector_address_at

from tests.test_nebula_blocks_db import _clouds, _link_to_field, _plant_field, _store_nebula

pytestmark = pytest.mark.db

EDGE_PC = float(tuning.DEFAULT_SECTOR_EDGE_PC)
BOX = ((600.0, 0.0, -250.0), (1000.0, 400.0, 350.0))


def _uid(number):
    return number.to_bytes(12, "big")


def _scatter(conn, kind, position, uid, built=None):
    shell, band, wedge = sector_address_at(position, EDGE_PC)
    conn.execute(
        "INSERT INTO phenomenon_scatter (shell_index, band_index, wedge_index, kind, subtype, position_x_mpc, position_y_mpc, position_z_mpc, seed,"
        " built_at, uid) VALUES (?, ?, ?, ?, NULL, ?, ?, ?, 1, ?, ?)", (shell, band, wedge, kind, *(round(v * 1000) for v in position), built, uid))


@pytest.fixture
def galaxy(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("DELETE FROM phenomenon_scatter")
        _plant_field(conn, _clouds(5, 60))
        cloud = field_clouds.field_clouds_in_box(conn, *BOX, 200)[0]
        linked = _store_nebula(conn, "N", cloud["radius_pc"] / 0.3066, (cloud["x"], cloud["y"], cloud["z"]))
        _link_to_field(conn, linked)         # a field cloud a sector stored: its row is stored and shares the nebula's Object ID
        other = _store_nebula(conn, "M", 40.0, (760.0, 120.0, 10.0))
        conn.execute("UPDATE nebulae SET uid = ? WHERE id = ?", (_uid(other + 1000), other))
        with conn:
            # a scattered nebula that already is the stored `other` (the same Object ID) but whose row was never marked built, a planetary nebula and a
            # remnant that are nobody's yet, and one that is built
            _scatter(conn, "planetary-nebula", (780.0, 130.0, 12.0), _uid(other + 1000))
            _scatter(conn, "planetary-nebula", (700.0, 50.0, -20.0), _uid(7001))
            _scatter(conn, "supernova-remnant", (900.0, 300.0, 100.0), _uid(7002))
            _scatter(conn, "supernova-remnant", (850.0, 200.0, 40.0), _uid(7003), built="2026-01-01 00:00:00")
        conn.commit()
        return {"config": mysql_config, "linked": linked, "other": other, "cloud": cloud}
    finally:
        conn.close()


def test_one_read_gives_every_nebula_kind_of_a_box_with_no_object_twice(galaxy):
    conn = store.get_connection(galaxy["config"])
    try:
        found = nebula_index.nebulae_in_box(conn, *BOX, EDGE_PC, 500)
        by_source = {source: [entry for entry in found if entry["source"] == source] for source in (nebula_index.STORED, nebula_index.FIELD, nebula_index.SCATTERED)}
        ids = [entry["object_id"] for entry in found if entry["object_id"] is not None]
        assert len(ids) == len(set(ids)), "no Object ID is listed twice"
        # the stored nebulae come first, the field clouds then, the scattered marks last
        assert [entry["source"] for entry in found] == sorted((entry["source"] for entry in found), key=(
            nebula_index.STORED, nebula_index.FIELD, nebula_index.SCATTERED).index)
        assert {entry["id"] for entry in by_source["stored"]} == {galaxy["linked"], galaxy["other"]}
        assert galaxy["cloud"]["id"] not in {entry["id"] for entry in by_source["field"]}
        assert len(by_source["field"]) == len(_field_only(conn)) and by_source["field"]
        # the scattered row that is the stored `other` is not a mark of its own; the two unbuilt ones that are nobody's yet are, the built one is not
        assert [(entry["type"], entry["object_id"]) for entry in by_source["scattered"]] == [
            ("planetary_nebula", _uid(7001).hex()), ("supernova_remnant", _uid(7002).hex())]
        assert {entry["type"] for entry in found} == {"nebula", "planetary_nebula", "supernova_remnant"}
    finally:
        conn.close()


def _field_only(conn):
    return field_clouds.field_clouds_in_box(conn, *BOX, 500)


def test_the_readers_alone_do_list_the_repeat_the_index_drops(galaxy):
    """The scattered reader leaves out only a built row, so the row that holds the stored nebula's Object ID is in its answer: the index is what drops it."""
    conn = store.get_connection(galaxy["config"])
    try:
        alone = query.galaxy_scattered_nebulae_in_box(conn, *BOX, EDGE_PC, 500)
        assert len(alone) == 3
        assert len([entry for entry in nebula_index.nebulae_in_box(conn, *BOX, EDGE_PC, 500) if entry["source"] == nebula_index.SCATTERED]) == 2
    finally:
        conn.close()


def test_a_limit_caps_each_source_and_a_box_with_nothing_lists_nothing(galaxy):
    conn = store.get_connection(galaxy["config"])
    try:
        few = nebula_index.nebulae_in_box(conn, *BOX, EDGE_PC, 1)
        assert max(len([e for e in few if e["source"] == s]) for s in ("stored", "field", "scattered")) == 1
        assert nebula_index.nebulae_in_box(conn, (5000.0, 5000.0, 5000.0), (5100.0, 5100.0, 5100.0), EDGE_PC, 50) == []
    finally:
        conn.close()
