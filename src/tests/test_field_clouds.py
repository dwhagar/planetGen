"""
MAP.178: a node of the Galaxy Map lists the field clouds (`nebula_field`) no sector has stored yet beside the stored nebulae, so a cloud is drawn
before its sector is generated, and once a sector stores one the stored nebula is listed and the field row is not.
"""

import numpy as np
import pytest

from planetgen import tuning
from planetgen.db import field_clouds, nodes, store
from planetgen.galaxy.drill import DrillBlock, format_drill_key
from planetgen.galaxy import nebula_blocks

from tests.test_nebula_blocks_db import _clouds, _link_to_field, _plant_field, _store_nebula

pytestmark = pytest.mark.db

EDGE_PC = float(tuning.DEFAULT_SECTOR_EDGE_PC)
BOX = ((600.0, 0.0, -250.0), (1000.0, 400.0, 350.0))


@pytest.fixture
def planted(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        _plant_field(conn, _clouds(5, 60))
    finally:
        conn.close()
    return mysql_config


def test_a_box_lists_the_field_clouds_that_reach_it_largest_first_with_the_key_of_their_row(planted):
    conn = store.get_connection(planted)
    try:
        lo, hi = BOX
        found = field_clouds.field_clouds_in_box(conn, lo, hi, 200)
        assert found and found == sorted(found, key=lambda c: (-c["radius_pc"], c["id"]))
        for cloud in found:
            assert cloud["type"] == "nebula" and cloud["field"] is True and cloud["name"] is None
            assert cloud["descriptor"] == "dark" and cloud["class"] in "MN" and cloud["id"].startswith("field-")
            gap = max(max(lo[i] - cloud[a], 0.0, cloud[a] - hi[i]) for i, a in enumerate("xyz"))
            assert gap <= cloud["radius_pc"]
        assert len({c["id"] for c in found}) == len(found)
        assert field_clouds.field_clouds_in_box(conn, lo, hi, 3) == found[:3]
    finally:
        conn.close()


def test_a_stored_field_cloud_is_listed_once_as_the_stored_nebula(planted):
    conn = store.get_connection(planted)
    try:
        lo, hi = BOX
        before = field_clouds.field_clouds_in_box(conn, lo, hi, 200)
        row = before[0]
        stored = _store_nebula(conn, "N", row["radius_pc"] / 0.3066, (row["x"], row["y"], row["z"]))
        _link_to_field(conn, stored)
        after = field_clouds.field_clouds_in_box(conn, lo, hi, 200)
        assert row["id"] not in {c["id"] for c in after}
        assert len(after) == len(before) - 1
    finally:
        conn.close()


def test_a_box_wider_than_the_reach_lists_none(planted):
    conn = store.get_connection(planted)
    try:
        assert field_clouds.field_clouds_in_box(conn, (-2000.0, -2000.0, -500.0), (2000.0, 2000.0, 500.0), 200) == []
        assert field_clouds.MAX_REACH_PC > 0
    finally:
        conn.close()


def test_a_node_lists_its_field_clouds_among_its_clouds_largest_first(planted, monkeypatch):
    monkeypatch.setattr(nodes, "_outer_shell", lambda conn: 1000)    # a galaxy whose outline holds the planted clouds
    conn = store.get_connection(planted)
    try:
        lo, hi = BOX
        first = field_clouds.field_clouds_in_box(conn, lo, hi, 1)[0]
        shell, band, wedge = nebula_blocks.locate(np.array([first["x"]]), np.array([first["y"]]), np.array([first["z"]]), 27, EDGE_PC)
        block = DrillBlock(27, int(shell[0]), int(band[0]), int(wedge[0]))
        node = nodes.galaxy_nodes(conn, [format_drill_key(block)])["nodes"][format_drill_key(block)]
        clouds = node["clouds"]
        assert clouds and first["id"] in {c["id"] for c in clouds} and any(c.get("field") for c in clouds)
        radii = [c["radius_pc"] for c in clouds]
        assert radii == sorted(radii, reverse=True) and len(clouds) <= nodes.NODE_CLOUDS
    finally:
        conn.close()
