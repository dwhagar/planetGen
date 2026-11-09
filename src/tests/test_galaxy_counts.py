# tests/test_galaxy_counts.py

"""
DB.121: the galaxy's counts are kept as they change (`db.galaxy_counts`). A sector's rows (`sector_facets`) and the one totals row
(`galaxy_totals`) are right after a sector is generated, after an edit, after the cell is emptied, and a read that uses them gives what
counting the tables gives. "Right" is what `galaxy_counts.problems` and a rebuild from the tables say.
"""

import pytest

from planetgen.db import edits, galaxy_counts, query, store
from tests.chart_support import chart_all
from tests.test_block_totals_follow_edits import ADDRESSES, _single_star_system

pytestmark = pytest.mark.db


@pytest.fixture
def galaxy(mysql_config):
    """Two generated sectors, each holding a system and stating its figures (`record_sector_stats`)."""
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
        chart_all(conn)
    yield conn, sectors
    conn.close()


def _row(conn):
    row = conn.execute("SELECT * FROM galaxy_totals WHERE id = 1").fetchone()
    conn.commit()
    return dict(row)


def _tables(conn, table):
    n = conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]
    conn.commit()
    return int(n)


def test_a_generated_sector_has_its_rows_and_the_totals_follow(galaxy):
    conn, sectors = galaxy
    assert galaxy_counts.usable(conn) and galaxy_counts.problems(conn) == []
    row = _row(conn)
    assert row["systems"] == _tables(conn, "star_systems") == 2
    assert row["stars"] == _tables(conn, "stars") and row["planets"] == _tables(conn, "planets") and row["moons"] == _tables(conn, "moons")
    assert galaxy_counts.totals(conn)["systems"] == 2
    kept = conn.execute("SELECT COUNT(*) AS n FROM sector_facets WHERE facet = 'total' AND value = 'systems'").fetchone()["n"]
    assert kept == 2


def test_the_totals_row_is_the_sum_of_the_sectors_rows_and_the_stats_rows(galaxy):
    conn, _sectors = galaxy
    row = _row(conn)
    stats = conn.execute("SELECT SUM(total_mass_kg) AS m, SUM(total_luminosity_sol) AS l FROM sector_stats WHERE bright_level_sol = 0").fetchone()
    assert row["mass_kg"] == pytest.approx(float(stats["m"] or 0)) and row["luminosity_sol"] == pytest.approx(float(stats["l"] or 0))
    assert galaxy_counts.problems(conn) == []


def test_a_rebuild_gives_the_rows_the_save_kept(galaxy):
    conn, _sectors = galaxy
    kept = {name: galaxy_counts.facet(conn, name) for name in galaxy_counts.FACETS}
    totals = _row(conn)
    totals.pop("updated_at")
    assert galaxy_counts.rebuild(conn) == 2
    again = _row(conn)
    again.pop("updated_at")
    assert again == pytest.approx(totals)
    assert {name: galaxy_counts.facet(conn, name) for name in galaxy_counts.FACETS} == kept


def test_the_reads_give_what_counting_the_tables_gives(galaxy):
    conn, _sectors = galaxy
    stored = query._search_facet_type(conn), query._search_facet_spectral(conn), query._search_facet_class(conn), query._search_facet_body(conn)
    assert query.count_systems(conn) == 2
    conn.execute("DELETE FROM galaxy_totals")                       # not whole any more: the reads count the tables
    conn.execute("DELETE FROM sector_facets")
    conn.commit()
    assert not galaxy_counts.usable(conn)
    assert (query._search_facet_type(conn), query._search_facet_spectral(conn), query._search_facet_class(conn),
            query._search_facet_body(conn)) == stored
    assert query.systems_facets(conn)["binary"] and query.count_systems(conn) == 2


def test_the_unfiltered_system_menus_come_from_the_stored_rows(galaxy):
    conn, _sectors = galaxy
    stored = query._stored_systems_facets(conn)
    assert stored is not None
    conn.execute("DELETE FROM galaxy_totals")
    conn.execute("DELETE FROM sector_facets")
    conn.commit()
    counted = query.systems_facets(conn)
    assert stored == counted


def test_an_edit_of_a_star_moves_the_totals(galaxy):
    conn, sectors = galaxy
    system_id, system = _single_star_system(conn, sectors[0])
    star = system.star
    star.luminosity = star.luminosity * 3
    before = _row(conn)["luminosity_sol"]
    with conn:
        edits.save_system_edits(conn, system_id, system, stars=[star])
    assert _row(conn)["luminosity_sol"] > before
    assert galaxy_counts.problems(conn) == []


def test_deleting_a_system_takes_it_out_of_the_totals(galaxy):
    conn, sectors = galaxy
    system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ? ORDER BY id LIMIT 1", (sectors[0],)).fetchone()["id"]
    conn.commit()
    with conn:
        conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,))
        store.touch_sector(conn, sectors[0])
        store.refresh_edited_sectors(conn, [sectors[0]])
    assert _row(conn)["systems"] == 1
    assert galaxy_counts.problems(conn) == []


def test_emptying_a_cell_takes_its_rows_out(galaxy):
    conn, sectors = galaxy
    with conn:
        edits.delete_sector_with_contents(conn, sectors[0])
    assert _row(conn)["systems"] == 1
    assert galaxy_counts.usable(conn) and galaxy_counts.problems(conn) == []
    left = conn.execute("SELECT COUNT(*) AS n FROM sector_facets WHERE shell_index = 0 AND band_index = 0 AND wedge_index = 0").fetchone()["n"]
    assert left == 0


def test_a_system_in_no_sector_is_counted_when_asked(galaxy):
    """A standalone system is in no sector's rows: `outside` counts it, so the totals still equal the tables."""
    from tests.test_storage_checksum import _system
    conn, _sectors = galaxy
    system, config = _system()
    with conn:
        store.insert_star_system(conn, system, config, sector_id=None)
    assert galaxy_counts.totals(conn)["systems"] == 3 == _tables(conn, "star_systems")
    assert galaxy_counts.problems(conn) == []
    assert query._stored_systems_facets(conn)["placement"] == [{"value": "sector", "count": 2}, {"value": "standalone", "count": 1}]


def test_rows_that_are_not_whole_are_reported_and_not_used(galaxy):
    conn, _sectors = galaxy
    with conn:
        conn.execute("DELETE FROM sector_facets WHERE shell_index = 0 AND band_index = 0 AND wedge_index = 0")
    assert not galaxy_counts.usable(conn) and galaxy_counts.totals(conn) is None and galaxy_counts.facet(conn, "class") is None
    assert galaxy_counts.problems(conn)
    galaxy_counts.rebuild(conn)
    assert galaxy_counts.usable(conn) and galaxy_counts.problems(conn) == []


def test_a_total_that_drifted_is_reported(galaxy):
    conn, _sectors = galaxy
    with conn:
        conn.execute("UPDATE galaxy_totals SET systems = systems + 5 WHERE id = 1")
    assert any("systems" in problem for problem in galaxy_counts.problems(conn))
