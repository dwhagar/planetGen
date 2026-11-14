# tests/test_migration_habitable_worlds.py

"""
MAP.143: migration 0139 gives `sector_stats` and `galaxy_block_stats` a `habitable_worlds` count. A sector filled after it
counts the planets and moons of its systems that are Blue or Green overall (PHI-6), the block totals add them up and the stage API reports them; a
galaxy made at v134 gets the counts from the scores it already holds.
"""

import pytest

from planetgen.db import block_stats, chart_ledger, query, store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import habitability
from planetgen.util import draw, generator_epoch as generatorEpoch

pytestmark = pytest.mark.fresh_database

ADDRESS = (1, 0, 0)
CENTER_PC = (1000.0, 0.0, 0.0)


def _fill_one_sector(conn):
    """A sector at `ADDRESS` holding one system with planets, its stats recorded, the stored scores set so that exactly one world
    (the system's first planet) is Green; returns (sector id, the number of bodies the stored scores put at Blue or Green overall).
    The scores are set here, not drawn: how often a generated world scores Green moves with the generator."""
    conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES ('Habit', 13046, ?, ?, ?)",
                 ADDRESS)
    sector_id = conn.execute("SELECT id FROM sectors WHERE name = 'Habit'").fetchone()["id"]
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    with draw.bound(0):
        system = StarSystem(system_config=cfg)
    system_id = store.insert_star_system(conn, system, cfg, sector_id=sector_id)
    conn.execute("UPDATE planets SET phi4 = 0.1 WHERE star_system_id = ?", (system_id,))
    conn.execute("UPDATE moons SET phi4 = 0.1 WHERE star_system_id = ?", (system_id,))
    first = conn.execute("SELECT MIN(id) AS id FROM planets WHERE star_system_id = ?", (system_id,)).fetchone()["id"]
    assert first is not None, "the system has no planet to score"
    conn.execute("UPDATE planets SET phi4 = ? WHERE id = ?", (0.9, first))
    expected = _scored_green_or_better(conn, system_id)
    assert expected == 1
    chart_ledger.record(conn, ADDRESS, None, generatorEpoch.generator_epoch())      # a filled sector is a charted one (MAP.402)
    store.record_sector_stats(conn, sector_id, ADDRESS, CENTER_PC)
    block_stats.refresh_chain(conn, ADDRESS)
    conn.commit()
    return sector_id, expected


def _scored_green_or_better(conn, system_id):
    """The planets and moons of one system with an overall score that `habitability.tier` calls Blue or Green."""
    scores = [row["phi4"] for table in ("planets", "moons") for row in conn.execute(
        f"SELECT phi4 FROM {table} WHERE star_system_id = ? AND phi4 IS NOT NULL", (system_id,)).fetchall()]
    return sum(1 for score in scores if habitability.tier(score) in habitability.POTENTIALLY_HABITABLE_TIERS)


def test_a_filled_sector_counts_its_habitable_worlds_and_the_blocks_above_it_add_them_up(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, expected = _fill_one_sector(conn)
        assert expected > 0
        assert store.sector_habitable_worlds(conn, sector_id) == expected
        stats = query.sector_stats(conn, *ADDRESS)
        assert stats["habitable_worlds"] == expected
        totals = conn.execute("SELECT m, habitable_worlds FROM galaxy_block_stats").fetchall()
        assert totals and all(row["habitable_worlds"] == expected for row in totals)
        stage = query.galaxy_stage(conn, None)
        assert [child["stats"]["habitable_worlds"] for child in stage["children"]] == [expected]
        store.forget_sector_fill(conn, ADDRESS)
        assert conn.execute("SELECT habitable_worlds FROM sector_stats").fetchone()["habitable_worlds"] is None
        # a forgotten fill is an uncharted cell again (MAP.402), so the blocks above hold no total for it: no row, or a zero
        assert {row["habitable_worlds"] for row in conn.execute("SELECT habitable_worlds FROM galaxy_block_stats").fetchall()} <= {0}
    finally:
        conn.close()


def test_only_blue_and_green_worlds_are_counted(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, expected = _fill_one_sector(conn)
        assert store.sector_habitable_worlds(conn, sector_id) == expected
        conn.execute("UPDATE planets SET phi4 = ? WHERE phi4 >= ?", (0.5, habitability.POTENTIALLY_HABITABLE_MIN_SCORE))
        conn.execute("UPDATE moons SET phi4 = NULL WHERE phi4 >= ?", (habitability.POTENTIALLY_HABITABLE_MIN_SCORE,))
        assert store.sector_habitable_worlds(conn, sector_id) == 0
    finally:
        conn.rollback()
        conn.close()


def test_a_re_score_that_moves_a_world_off_green_moves_the_sector_and_the_blocks_above_it(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, expected = _fill_one_sector(conn)
        system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
        assert store.recount_habitable_worlds(conn, [system_id]) == []    # nothing moved yet
        conn.execute("UPDATE planets SET phi4 = 0.1 WHERE star_system_id = ?", (system_id,))
        conn.execute("UPDATE moons SET phi4 = 0.1 WHERE star_system_id = ?", (system_id,))
        assert store.recount_habitable_worlds(conn, [system_id]) == [ADDRESS]
        assert block_stats.refresh_sectors(conn, [ADDRESS]) > 0
        assert conn.execute("SELECT habitable_worlds FROM sector_stats").fetchone()["habitable_worlds"] == 0
        assert {row["habitable_worlds"] for row in conn.execute("SELECT habitable_worlds FROM galaxy_block_stats").fetchall()} == {0}
    finally:
        conn.rollback()
        conn.close()


def test_the_line_is_where_green_begins():
    floor = habitability.POTENTIALLY_HABITABLE_MIN_SCORE
    assert habitability.tier(floor) in habitability.POTENTIALLY_HABITABLE_TIERS
    assert habitability.tier(floor - 1e-9) not in habitability.POTENTIALLY_HABITABLE_TIERS
