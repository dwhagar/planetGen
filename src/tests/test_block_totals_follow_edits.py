# tests/test_block_totals_follow_edits.py

"""
DB.83: the totals of the Galaxy Map's blocks (`galaxy_block_stats`) follow what the generated sectors hold. A sector's own figures
(`sector_stats`) and the blocks above it are right after a sector is generated, after an edit that changes the systems, stars or
worlds in it (`store.refresh_edited_sectors`, called by `db.edits.save_system_edits` and the system delete), and after a sector is
deleted. "Right" is what a rebuild from the sectors gives.
"""

import pytest

from planetgen.api import routes
from planetgen.db import block_stats, edits, store

pytestmark = pytest.mark.db

FIGURES = ("actual_systems", "actual_stars", "mean_temperature_k", "mean_luminosity_sol", "mean_age_gy", "total_luminosity_sol",
           "habitable_worlds")


def _block_rows(conn):
    rows = conn.execute("SELECT * FROM galaxy_block_stats ORDER BY m, shell_index, band_index, wedge_index").fetchall()
    conn.commit()
    return [tuple(row.values()) if hasattr(row, "values") else tuple(row) for row in rows]


def _assert_figures_follow_the_rows(conn):
    """Every generated sector's stats row states what its rows hold now."""
    rows = conn.execute("SELECT s.id AS id, st.* FROM sectors s JOIN sector_stats st ON st.shell_index = s.shell_index"
                        " AND st.band_index = s.band_index AND st.wedge_index = s.wedge_index"
                        " WHERE st.bright_level_sol = 0").fetchall()
    conn.commit()
    assert rows
    for row in rows:
        held = dict(zip(FIGURES, store._content_figures(conn, row["id"])))
        for name in FIGURES:
            assert row[name] == pytest.approx(held[name]), (row["id"], name)


def _assert_blocks_follow_the_figures(conn):
    """The stored block rows are the ones a rebuild from the sectors' stats gives."""
    stored = _block_rows(conn)
    assert stored
    block_stats.rebuild(conn)
    conn.commit()
    rebuilt = _block_rows(conn)
    assert len(stored) == len(rebuilt)
    for have, want in zip(stored, rebuilt):
        assert have == pytest.approx(want)


ADDRESSES = ((0, 0, 0), (2, 1, 1))


@pytest.fixture
def galaxy(mysql_config):
    """Two generated sectors at grid addresses, each holding a system and stating its figures (`record_sector_stats`), and the
    blocks above them."""
    from tests.test_storage_checksum import _system
    conn = store.get_connection(mysql_config)
    sectors = []
    with conn:
        for n, address in enumerate(ADDRESSES):
            system, config = _system()
            sector_id = conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) "
                                     "VALUES (?, 4000, ?, ?, ?)", (f"Sector {n}", *address)).lastrowid
            store.insert_star_system(conn, system, config, sector_id=sector_id)
            store.record_sector_stats(conn, sector_id, address, (0.0, 0.0, 0.0))
            sectors.append(sector_id)
        block_stats.rebuild(conn)
    yield conn, sectors[0]
    conn.close()


def _single_star_system(conn, sector_id):
    for row in conn.execute("SELECT id FROM star_systems WHERE sector_id = ? ORDER BY id", (sector_id,)).fetchall():
        system = store.load_star_system(conn, row["id"])
        if system.binary_type is None and not hasattr(system.star, "event_horizon_radius_km"):
            return row["id"], system
    pytest.skip("the generated sector holds no single-star system")


def test_a_generated_sector_and_the_blocks_above_it_are_right_to_start_with(galaxy):
    conn, sector_id = galaxy
    block_stats.refresh_chain(conn, (0, 0, 0))
    conn.commit()
    _assert_figures_follow_the_rows(conn)
    _assert_blocks_follow_the_figures(conn)


def test_an_edit_of_a_star_changes_the_sector_and_the_blocks_above_it(galaxy):
    conn, sector_id = galaxy
    block_stats.refresh_chain(conn, (0, 0, 0))
    conn.commit()
    system_id, system = _single_star_system(conn, sector_id)
    before = _block_rows(conn)
    star = system.star
    star.age = star.age * 2 + 1
    star.luminosity = star.luminosity * 3
    with conn:
        edits.save_system_edits(conn, system_id, system, stars=[star])
    _assert_figures_follow_the_rows(conn)
    _assert_blocks_follow_the_figures(conn)
    assert _block_rows(conn) != before                          # the edit reached the blocks


def test_deleting_a_system_takes_it_out_of_the_sector_and_the_blocks(galaxy, mysql_config):
    conn, sector_id = galaxy
    block_stats.refresh_chain(conn, (0, 0, 0))
    conn.commit()
    systems = conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["n"]
    system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ? ORDER BY id LIMIT 1", (sector_id,)).fetchone()["id"]
    conn.commit()
    assert systems >= 1
    with conn:
        conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,))
        store.touch_sector(conn, sector_id)
        store.refresh_edited_sectors(conn, [sector_id])
    row = conn.execute("SELECT actual_systems FROM sector_stats WHERE shell_index = 0 AND band_index = 0 AND wedge_index = 0").fetchone()
    assert row["actual_systems"] == systems - 1
    _assert_figures_follow_the_rows(conn)
    _assert_blocks_follow_the_figures(conn)


def test_without_the_refresh_an_edit_leaves_the_figures_behind(galaxy, monkeypatch):
    """What DB.83 found: nothing kept the sector's figures, and so the blocks above it, after an edit."""
    conn, sector_id = galaxy
    monkeypatch.setattr(store, "refresh_edited_sectors", lambda *args, **kwargs: [])
    system_id, system = _single_star_system(conn, sector_id)
    system.star.luminosity = system.star.luminosity * 3
    with conn:
        edits.save_system_edits(conn, system_id, system, stars=[system.star])
    held = store._content_figures(conn, sector_id)[5]
    stated = conn.execute("SELECT total_luminosity_sol FROM sector_stats WHERE shell_index = 0 AND band_index = 0"
                          " AND wedge_index = 0").fetchone()["total_luminosity_sol"]
    assert stated != pytest.approx(held)


def test_the_delete_system_route_refreshes_the_figures(galaxy, monkeypatch):
    """The route's own transaction calls the refresh; the helper it calls does the work tested above."""
    import inspect
    assert "refresh_edited_sectors" in inspect.getsource(routes.delete_system)


def test_a_sector_that_is_not_a_generated_grid_cell_is_left_alone(galaxy):
    conn, _sector_id = galaxy
    with conn:
        loose = conn.execute("INSERT INTO sectors (name, edge_mpc) VALUES ('Loose', 1000)").lastrowid
    assert store.refresh_edited_sectors(conn, [loose, None, 10 ** 9]) == []


def test_deleting_the_sector_takes_it_out_of_the_blocks(galaxy):
    conn, sector_id = galaxy
    block_stats.refresh_chain(conn, (0, 0, 0))
    conn.commit()
    before = _block_rows(conn)
    with conn:
        edits.delete_sector_with_contents(conn, sector_id)
    after = _block_rows(conn)
    assert [row[4] for row in before][-1] == 2 and [row[4] for row in after][-1] == 1       # the root block's sector count
    _assert_blocks_follow_the_figures(conn)
