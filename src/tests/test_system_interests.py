# tests/test_system_interests.py

"""
MAP.305, MAP.303: a star system worth a mark on the Galaxy Map has a row in `system_interests` for each reason it has (migration 0165:
a planet or moon rated Yellow or better overall on the PHI-6 score, a gas giant), the map's nodes list those systems for the blocks of 81
shells and smaller (`interests`), and a re-score keeps the rows right.
"""

import pytest

from planetgen.db import interests, query, store
from planetgen.db import nodes
from planetgen.galaxy import spherical
from planetgen.galaxy.drill import drill_chain_of, format_drill_key
from planetgen.physics.units import ly_to_pc
from planetgen.physics import habitability
from planetgen.db import block_stats
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.test_migration_habitable_worlds import ADDRESS, CENTER_PC

pytestmark = pytest.mark.fresh_database

POSITION_PC = (1000.0, 2.0, 3.0)
SCORES = (0.82, 0.55)
"""The overall scores the fixture system's first two planets are given: Green and Yellow."""


def _fill_one_sector(conn):
    """A sector at `ADDRESS` holding one system with at least two planets, the first scored Green and the second Yellow (set
    here, so the test does not hang on what the generator draws), its stats recorded; returns (sector id, None)."""
    conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES ('Interest', 13046, ?, ?, ?)",
                 ADDRESS)
    sector_id = conn.execute("SELECT id FROM sectors WHERE name = 'Interest'").fetchone()["id"]
    for seed in range(120):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        system_id = store.insert_star_system(conn, system, cfg, sector_id=sector_id)
        planets = [row["id"] for row in conn.execute("SELECT id FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
        if len(planets) >= len(SCORES):
            break
        conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,))
    else:
        pytest.fail("no system with two planets in 120 seeds")
    for planet_id, score in zip(planets, SCORES):
        conn.execute("UPDATE planets SET phi4 = ? WHERE id = ?", (score, planet_id))
    store.record_sector_stats(conn, sector_id, ADDRESS, CENTER_PC)
    block_stats.refresh_chain(conn, ADDRESS)
    conn.commit()
    return sector_id, None


def _place(conn, sector_id):
    """Puts the sector's one system at `POSITION_PC`; returns its id."""
    system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
    conn.execute("UPDATE star_systems SET galaxy_x_pc = ?, galaxy_y_pc = ?, galaxy_z_pc = ? WHERE id = ?",
                 (*POSITION_PC, system_id))
    store.recount_habitable_worlds(conn, [system_id])
    conn.commit()
    return system_id


def _box(around, half=1.0):
    return (tuple(v - half for v in around), tuple(v + half for v in around))


def _rows(conn, system_id):
    return {row["reason"]: row for row in conn.execute(
        "SELECT reason, n, best_score FROM system_interests WHERE star_system_id = ?", (system_id,)).fetchall()}


def _expected_habitable(conn, system_id):
    """How many planets and moons of the system are rated Yellow or better, and the best score among them."""
    floor = habitability.YELLOW_FLOOR
    scores = [row["phi4"] for table in ("planets", "moons") for row in conn.execute(
        f"SELECT phi4 FROM {table} WHERE star_system_id = ? AND phi4 >= ?", (system_id, floor)).fetchall()]
    return len(scores), max(scores) if scores else None


def _expected_giants(conn, system_id):
    return conn.execute("SELECT COUNT(*) AS n FROM planets WHERE star_system_id = ? AND body_type = 'g'", (system_id,)).fetchone()["n"]


def test_a_filled_sector_gives_each_system_a_row_for_each_reason(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        assert conn.execute("SELECT COUNT(*) AS n FROM system_interests").fetchone()["n"] == 0    # nowhere to draw it yet
        system_id = _place(conn, sector_id)
        rows = _rows(conn, system_id)
        count, best = _expected_habitable(conn, system_id)
        giants = _expected_giants(conn, system_id)
        assert count > 0, "the fixture system holds a world rated Yellow or better"
        assert rows[interests.HABITABLE]["n"] == count and rows[interests.HABITABLE]["best_score"] == pytest.approx(best)
        if giants:
            assert rows[interests.GAS_GIANT]["n"] == giants and rows[interests.GAS_GIANT]["best_score"] is None
        else:
            assert interests.GAS_GIANT not in rows
    finally:
        conn.close()


def test_the_box_query_lists_a_system_with_its_reasons_and_only_inside_the_box(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        lo, hi = _box(POSITION_PC)
        found = query.galaxy_system_interests_in_box(conn, lo, hi)
        assert [s["id"] for s in found] == [system_id]
        assert (found[0]["x"], found[0]["y"], found[0]["z"]) == POSITION_PC
        assert found[0]["name"]
        count, best = _expected_habitable(conn, system_id)
        habitable = found[0]["reasons"][interests.HABITABLE]
        assert habitable["n"] == count and habitable["score"] == pytest.approx(best, abs=1e-3)
        assert habitable["tier"] == habitability.tier(best)
        far_lo, far_hi = _box((-5000.0, 0.0, 0.0))
        assert query.galaxy_system_interests_in_box(conn, far_lo, far_hi) == []
        conn.execute("DELETE FROM system_interests")
        assert query.galaxy_system_interests_in_box(conn, lo, hi) == []
    finally:
        conn.rollback()
        conn.close()


def test_the_limit_cuts_the_list(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        first = _place(conn, sector_id)
        lo, hi = _box(POSITION_PC)
        assert [s["id"] for s in query.galaxy_system_interests_in_box(conn, lo, hi, limit=1)] == [first]
        assert query.galaxy_system_interests_in_box(conn, lo, hi, limit=0) == []
    finally:
        conn.rollback()
        conn.close()


def test_a_small_node_carries_the_interests_and_a_big_one_does_not(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        # A node's box stops at the outline (the furthest stored sector's shell), so the system sits in its own sector.
        edge_pc = ly_to_pc(query.DEFAULT_SECTOR_EDGE_LY)
        position = spherical.sector_position_pc(*ADDRESS, edge_pc)
        conn.execute("UPDATE star_systems SET galaxy_x_pc = ?, galaxy_y_pc = ?, galaxy_z_pc = ? WHERE id = ?",
                     (*position, system_id))
        store.recount_habitable_worlds(conn, [system_id])
        chain = drill_chain_of(*spherical.sector_address_at(position, edge_pc))
        keys = [format_drill_key(chain[-1]), format_drill_key(chain[0])]
        found = nodes.galaxy_nodes(conn, keys)["nodes"]
        assert [s["id"] for s in found[keys[0]]["interests"]] == [system_id]
        assert interests.HABITABLE in found[keys[0]]["interests"][0]["reasons"]
        assert found[keys[1]]["interests"] == []
        assert "worlds" not in found[keys[0]]
    finally:
        conn.close()


def test_a_re_score_moves_the_systems_rows(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        conn.execute("UPDATE planets SET phi4 = 0.1 WHERE star_system_id = ?", (system_id,))
        conn.execute("UPDATE moons SET phi4 = 0.1 WHERE star_system_id = ?", (system_id,))
        store.recount_habitable_worlds(conn, [system_id])
        assert interests.HABITABLE not in _rows(conn, system_id)
        # Yellow is enough for the map (the habitable rule, Blue or Green, is the sector count's).
        conn.execute("UPDATE planets SET phi4 = ? WHERE star_system_id = ? ORDER BY id LIMIT 1",
                     (habitability.YELLOW_FLOOR, system_id))
        store.recount_habitable_worlds(conn, [system_id])
        row = _rows(conn, system_id)[interests.HABITABLE]
        assert row["n"] == 1 and row["best_score"] == pytest.approx(habitability.YELLOW_FLOOR)
    finally:
        conn.rollback()
        conn.close()


def test_a_gas_giant_is_a_reason_of_its_own(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        conn.execute("UPDATE planets SET body_type = 't' WHERE star_system_id = ?", (system_id,))
        store.recount_habitable_worlds(conn, [system_id])
        assert interests.GAS_GIANT not in _rows(conn, system_id)
        conn.execute("UPDATE planets SET body_type = 'g' WHERE star_system_id = ? ORDER BY id LIMIT 2", (system_id,))
        store.recount_habitable_worlds(conn, [system_id])
        row = _rows(conn, system_id)[interests.GAS_GIANT]
        assert row["n"] == min(2, conn.execute("SELECT COUNT(*) AS n FROM planets WHERE star_system_id = ?", (system_id,)).fetchone()["n"])
        lo, hi = _box(POSITION_PC)
        found = query.galaxy_system_interests_in_box(conn, lo, hi)
        assert found[0]["reasons"][interests.GAS_GIANT] == {"n": row["n"]}
    finally:
        conn.rollback()
        conn.close()


def test_a_v164_galaxy_gets_its_rows_from_the_planets_and_moons_it_holds(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        want = {reason: (row["n"], row["best_score"]) for reason, row in _rows(conn, system_id).items()}
    finally:
        conn.close()
    # make it a v164 galaxy again: the table gone, the old one back, the rows kept
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DROP TABLE system_interests")
        conn.execute("CREATE TABLE habitable_systems (star_system_id BIGINT UNSIGNED NOT NULL PRIMARY KEY,"
                     " habitable_worlds SMALLINT UNSIGNED NOT NULL, galaxy_x_pc DOUBLE NOT NULL, galaxy_y_pc DOUBLE NOT NULL,"
                     " galaxy_z_pc DOUBLE NOT NULL)")
        conn.execute("UPDATE alembic_version SET version_num = '0164'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (164)")
        conn.commit()
        assert store._stored_version(conn, "schema_migrations") == 164
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        got = {reason: (row["n"], row["best_score"]) for reason, row in _rows(conn, system_id).items()}
        assert got == want
        assert conn.execute("SHOW TABLES LIKE 'habitable_systems'").fetchone() is None
    finally:
        conn.close()
