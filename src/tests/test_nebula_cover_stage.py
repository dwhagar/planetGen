"""
MAP.174: the nebula cover in the Galaxy Map's stage and window statistics. The share of a cell inside a nebula is read
from the stored aggregates (MAP.172), which are made from the seed's field and the stored nebulae and need no sector,
so a cell that holds generated sectors and one that holds none give the same figure for the same nebula, and filling a
cell (the sectors saved, its nebulae stored, the aggregates caught up) does not change it.
"""

import numpy as np
import pytest

from planetgen import tuning
from planetgen.db import block_stats, nebula_blocks, query, store
from planetgen.galaxy import nebula_blocks as rules
from planetgen.galaxy.drill import DRILL_STEP, STAGE_TOP, DrillBlock, drill_block_bounds, drill_chain_of, drill_children, format_drill_key
from planetgen.web import galaxy_views
from planetgen.web.lib import apiclient

from tests.test_block_stats import _addresses, galaxy  # noqa: F401 -- fixture
from tests.test_nebula_blocks_db import _arrays, _link_to_field, _plant_field, _store_nebula

pytestmark = pytest.mark.db

EDGE_PC = float(tuning.DEFAULT_SECTOR_EDGE_PC)
CLOUD = ("N", 8.0)
"""tuple: The class and radius in light years of the cloud put in the middle of each block."""


def _middle(block):
    """The point in the middle of `block` (parsecs): the middle radius, colatitude and longitude of its bounds."""
    (inner, outer), (north, south), (start, end) = drill_block_bounds(block, EDGE_PC)
    radius, colatitude, longitude = (inner + outer) / 2, (north + south) / 2, (start + end) / 2
    return (radius * np.sin(colatitude) * np.cos(longitude), radius * np.sin(colatitude) * np.sin(longitude),
            radius * np.cos(colatitude))


def _blocks():
    """`(filled, unfilled, parent)`: a level-3 block that holds a generated sector, a sibling that holds none, and the
    level-9 block above both."""
    held = {key[1:] for key in _filled_keys() if key[0] == DRILL_STEP}
    for address in _addresses():
        chain = drill_chain_of(*address)
        filled, parent = chain[-2], chain[-3]
        for _shell, group in drill_children(parent):
            for sibling in group:
                if (sibling.shell, sibling.band, sibling.wedge) not in held:
                    return filled, sibling, parent
    raise AssertionError("no level-3 block without a sector beside one with")


def _filled_keys():
    keys = set()
    for address in _addresses():
        keys.update(tuple(block) for block in drill_chain_of(*address)[:-1] if block.m <= STAGE_TOP)
    return keys


def _expected(*blocks):
    """The cover the rules give each of `blocks` for one cloud in the middle of every one of them."""
    clouds = [(CLOUD[0], CLOUD[1], tuple(_middle(block))) for block in blocks]
    expected = {}
    codes, x, y, z, radius = _arrays(clouds)
    for m, level in rules.all_rows(x, y, z, radius, codes, EDGE_PC).items():
        for row in level.blocks:
            expected[(m, row["shell_index"], row["band_index"], row["wedge_index"])] = row["cover"]
    return expected, clouds


@pytest.fixture
def seeded(galaxy):  # noqa: F811 -- the fixture of test_block_stats
    """The galaxy of generated sectors with a cloud in the middle of a filled and of an unfilled level-3 block, and the
    aggregates made from the field."""
    filled, unfilled, parent = _blocks()
    expected, clouds = _expected(filled, unfilled)
    conn = store.get_connection(galaxy)
    try:
        _plant_field(conn, clouds)
        nebula_blocks.rebuild(conn, EDGE_PC)
        conn.commit()
    finally:
        conn.close()
    return galaxy, filled, unfilled, parent, expected, clouds


def _cells(conn, at):
    return {tuple(cell[:3]): cell[3] for cell in query.galaxy_nebula_cover(conn, at)["cells"]}


def _through_web(conn, at, monkeypatch):
    """The stage of `at` as the web layer serves it: the API's stage with the cover added."""
    monkeypatch.setattr(apiclient, "get_galaxy_nebula_cover",
                        lambda db, where=None, g=None: query.galaxy_nebula_cover(conn, where, g, EDGE_PC))
    return galaxy_views._with_nebula_cover(query.galaxy_stage(conn, at), "db", at)


def test_filled_and_unfilled_cells_give_the_cover_the_rules_give_the_field(seeded, monkeypatch):
    config, filled, unfilled, parent, expected, _clouds = seeded
    conn = store.get_connection(config)
    try:
        at = format_drill_key(parent)
        stage = _through_web(conn, at, monkeypatch)
        covers = {tuple(cell[:3]): cell[3] for cell in stage["nebula_cover"]}
        for block in (filled, unfilled):
            key = (block.shell, block.band, block.wedge)
            assert covers[key] == pytest.approx(expected[tuple(block)], rel=1e-9), block
        children = {(c["shell"], c["band"], c["wedge"]): c for c in stage["children"]}
        assert (unfilled.shell, unfilled.band, unfilled.wedge) not in children     # no sector: only in the list
        assert children[(filled.shell, filled.band, filled.wedge)]["stats"]["nebula_cover"] == pytest.approx(
            expected[tuple(filled)], rel=1e-9)
        assert all(0 < cover <= 1 for cover in covers.values())
        assert stage["child_m"] == DRILL_STEP
    finally:
        conn.close()


def test_filling_an_unfilled_cell_does_not_change_its_cover(seeded, monkeypatch):
    """Sectors saved in the block, the cloud it holds stored as a nebula the field knows, the aggregates caught up."""
    config, filled, unfilled, parent, expected, clouds = seeded
    at = format_drill_key(parent)
    conn = store.get_connection(config)
    try:
        before = _cells(conn, at)[(unfilled.shell, unfilled.band, unfilled.wedge)]
        _shell, group = drill_children(unfilled)[0]
        sector = group[0]
        conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES ('New', 4000, ?, ?, ?)",
                     (sector.shell, sector.band, sector.wedge))
        conn.execute("INSERT INTO sector_stats (shell_index, band_index, wedge_index, expected_systems, actual_systems,"
                     " actual_stars, mean_age_gy, total_luminosity_sol, bright_level_sol) VALUES (?, ?, ?, 2, 2, 3, 1.0, 10.0, 0)",
                     (sector.shell, sector.band, sector.wedge))
        letter, radius_ly, centre = next(c for c in clouds if tuple(c[2]) == tuple(_middle(unfilled)))
        _link_to_field(conn, _store_nebula(conn, letter, radius_ly, centre))
        block_stats.rebuild(conn)
        nebula_blocks.catch_up(conn, EDGE_PC)
        conn.commit()
        assert _cells(conn, at)[(unfilled.shell, unfilled.band, unfilled.wedge)] == pytest.approx(before, rel=1e-12)
        stage = _through_web(conn, at, monkeypatch)
        children = {(c["shell"], c["band"], c["wedge"]): c for c in stage["children"]}
        now = children[(unfilled.shell, unfilled.band, unfilled.wedge)]["stats"]["nebula_cover"]
        assert now == pytest.approx(before, rel=1e-12) and now == pytest.approx(expected[tuple(unfilled)], rel=1e-9)
    finally:
        conn.close()


def test_the_galaxy_stage_a_window_and_the_cells_of_a_sector(seeded):
    config, filled, unfilled, parent, expected, _clouds = seeded
    conn = store.get_connection(config)
    try:
        top = query.galaxy_nebula_cover(conn)
        assert top["at"] is None and top["child_m"] == 243
        assert top["cells"] and all(0 < cell[3] <= 1 for cell in top["cells"])
        window = query.galaxy_nebula_cover(conn, format_drill_key(filled), 3, EDGE_PC)
        assert window["child_m"] == DRILL_STEP and window["at"] == format_drill_key(filled)
        found = {(c[0], c[1], c[2]): c[3] for c in window["cells"]}
        assert found[(filled.shell, filled.band, filled.wedge)] == pytest.approx(expected[tuple(filled)], rel=1e-9)
        # The aggregates are of the blocks above a sector: a stage of sectors and a window of sectors have none.
        assert query.galaxy_nebula_cover(conn, format_drill_key(filled))["cells"] == []
        sector = DrillBlock(1, filled.shell, filled.band, filled.wedge)
        assert query.galaxy_nebula_cover(conn, format_drill_key(sector), 3, EDGE_PC)["cells"] == []
    finally:
        conn.close()


def test_a_malformed_key_or_width_is_refused(seeded):
    conn = store.get_connection(seeded[0])
    try:
        with pytest.raises(ValueError):
            query.galaxy_nebula_cover(conn, "nonsense")
        with pytest.raises(ValueError):
            query.galaxy_nebula_cover(conn, "27.3.5.0", 4, EDGE_PC)
    finally:
        conn.close()


def test_the_api_route_answers_the_same_cells(seeded, mysql_config):
    from planetgen.api.config import Config
    from planetgen.web.app import create_app

    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    _config, filled, unfilled, parent, expected, _clouds = seeded
    client = create_app(TestConfig).test_client()
    answer = client.get("/api/galaxy/nebula-cover", query_string={"at": format_drill_key(parent)}).get_json()
    assert answer["at"] == format_drill_key(parent) and answer["child_m"] == DRILL_STEP
    found = {(c[0], c[1], c[2]): c[3] for c in answer["cells"]}
    assert found[(unfilled.shell, unfilled.band, unfilled.wedge)] == pytest.approx(expected[tuple(unfilled)], rel=1e-9)
    assert client.get("/api/galaxy/nebula-cover", query_string={"at": "nonsense"}).status_code == 400
