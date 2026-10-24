# tests/test_habitable_system_markers.py

"""
MAP.220: a star system with a habitable world has a row in `habitable_systems` (migration 0148), the
Galaxy Map's tiles list the systems that have one from the finest levels (`worlds`), and a re-score keeps the count right.
"""

import pytest

from planetgen.db import query, store
from planetgen.galaxy import viewport
from planetgen.physics import habitability
from tests.test_migration_habitable_worlds import ADDRESS, _fill_one_sector

pytestmark = pytest.mark.fresh_database

POSITION_PC = (1000.0, 2.0, 3.0)


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


def test_a_filled_sector_gives_each_system_its_count(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, expected = _fill_one_sector(conn)
        assert conn.execute("SELECT COUNT(*) AS n FROM habitable_systems").fetchone()["n"] == 0    # nowhere to draw it yet
        _place(conn, sector_id)
        counts = [row["habitable_worlds"] for row in conn.execute(
            "SELECT h.habitable_worlds FROM habitable_systems h JOIN star_systems ss ON ss.id = h.star_system_id"
            " WHERE ss.sector_id = ?", (sector_id,)).fetchall()]
        assert counts == [expected] and expected > 0
    finally:
        conn.close()


def test_the_box_query_lists_a_system_with_habitable_worlds_and_only_inside_the_box(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        lo, hi = _box(POSITION_PC)
        found = query.galaxy_habitable_systems_in_box(conn, lo, hi)
        assert [(s["id"], s["worlds"]) for s in found] == [(system_id, expected)]
        assert (found[0]["x"], found[0]["y"], found[0]["z"]) == POSITION_PC
        assert found[0]["name"]
        far_lo, far_hi = _box((-5000.0, 0.0, 0.0))
        assert query.galaxy_habitable_systems_in_box(conn, far_lo, far_hi) == []
        conn.execute("DELETE FROM habitable_systems")
        assert query.galaxy_habitable_systems_in_box(conn, lo, hi) == []
    finally:
        conn.rollback()
        conn.close()


def test_the_limit_cuts_the_list(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        first = _place(conn, sector_id)
        lo, hi = _box(POSITION_PC)
        assert [s["id"] for s in query.galaxy_habitable_systems_in_box(conn, lo, hi, limit=1)] == [first]
        assert query.galaxy_habitable_systems_in_box(conn, lo, hi, limit=0) == []
    finally:
        conn.rollback()
        conn.close()


def test_a_fine_tile_carries_the_worlds_and_a_coarse_one_does_not(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        keys = viewport.tile_keys_containing(POSITION_PC)
        fine = query.galaxy_tiles(conn, [keys[-1], keys[0]])["tiles"]
        assert [(s["id"], s["worlds"]) for s in fine[keys[-1]]["worlds"]] == [(system_id, expected)]
        assert fine[keys[0]]["worlds"] == []
    finally:
        conn.close()


def test_a_re_score_moves_the_systems_count(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        floor = habitability.POTENTIALLY_HABITABLE_MIN_SCORE
        conn.execute("UPDATE planets SET phi4 = 0.1 WHERE star_system_id = ?", (system_id,))
        conn.execute("UPDATE moons SET phi4 = 0.1 WHERE star_system_id = ?", (system_id,))
        store.recount_habitable_worlds(conn, [system_id])
        assert conn.execute("SELECT habitable_worlds FROM habitable_systems WHERE star_system_id = ?",
                            (system_id,)).fetchone() is None
        conn.execute("UPDATE planets SET phi4 = ? WHERE star_system_id = ? ORDER BY id LIMIT 1", (floor, system_id))
        store.recount_habitable_worlds(conn, [system_id])
        assert conn.execute("SELECT habitable_worlds FROM habitable_systems WHERE star_system_id = ?",
                            (system_id,)).fetchone()["habitable_worlds"] == 1
    finally:
        conn.rollback()
        conn.close()


def test_a_v147_galaxy_gets_its_counts_from_the_planets_and_moons_it_holds(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, expected = _fill_one_sector(conn)
        _place(conn, sector_id)
    finally:
        conn.close()
    # make it a v147 galaxy again: the table gone, the rows kept
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DROP TABLE habitable_systems")
        conn.execute("UPDATE alembic_version SET version_num = '0147'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (147)")
        conn.commit()
        assert store._stored_version(conn, "schema_migrations") == 147
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        counts = [row["habitable_worlds"] for row in conn.execute(
            "SELECT h.habitable_worlds FROM habitable_systems h JOIN star_systems ss ON ss.id = h.star_system_id"
            " WHERE ss.sector_id = ?", (sector_id,)).fetchall()]
        assert counts == [expected]
    finally:
        conn.close()
