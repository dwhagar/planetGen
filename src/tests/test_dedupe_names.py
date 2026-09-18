# tests/test_dedupe_names.py

"""
Tests for `src/dedupeNames.py` -- the one-off backfill script that cleans
up sector/system/planet/moon name collisions in an existing database
(v22, `stellarObjects/nameUniqueness.py`). Live generation (via
`stellarObjects._db.py`'s `reserve_*_name`/`confirm_*_name`) already
prevents any *new* duplicate, so these tests simulate a "legacy" database
that predates that guarantee: insert normally (which the reservation
system already resolves), then manually revert the row(s)/registry to
look as if no resolution had ever happened, and confirm the script
recovers a duplicate-free state -- and that a second run is a no-op.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable.

Run with: pytest tests/test_dedupe_names.py
"""

from dedupeNames import dedupe_names

from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem


def test_dedupe_fixes_a_legacy_sector_duplicate_and_is_idempotent(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            id1 = _db.insert_sector(conn, SpaceSector(name="Sol"))
        with conn:
            id2 = _db.insert_sector(conn, SpaceSector(name="Sol"))

        # Revert to a literal duplicate, as if resolution had never run.
        with conn:
            conn.execute("UPDATE sectors SET name = ? WHERE id IN (?, ?)", ("Sol", id1, id2))
            conn.execute("DELETE FROM sector_name_registry WHERE base_name = ?", ("Sol",))

        before = {r["id"]: r["name"] for r in conn.execute("SELECT id, name FROM sectors").fetchall()}
        assert before[id1] == before[id2] == "Sol"
    finally:
        conn.close()

    counts = dedupe_names(mysql_config)
    assert counts["sectors"] == 1
    assert counts["star_systems"] == 0
    assert counts["planets_and_moons"] == 0

    conn = _db.get_connection(mysql_config)
    try:
        after = {r["id"]: r["name"] for r in conn.execute("SELECT id, name FROM sectors").fetchall()}
        assert len(set(after.values())) == 2
        assert {after[id1], after[id2]} == {"Alpha Sol", "Beta Sol"}
    finally:
        conn.close()

    # A second run against an already-clean database renames nothing.
    counts_again = dedupe_names(mysql_config)
    assert sum(counts_again.values()) == 0


def test_dedupe_fixes_a_legacy_system_vs_sector_collision(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            sector_id = _db.insert_sector(conn, SpaceSector(name="Venus"))
        with conn:
            system = StarSystem(system_config=SystemConfig())
            system.star.name = "Venus"
            system_id = _db.insert_star_system(conn, system, system.system_config)

        # Revert the system's own cross-level resolution, simulating a
        # database from before that check existed.
        with conn:
            conn.execute("UPDATE star_systems SET name = ? WHERE id = ?", ("Venus", system_id))
            conn.execute("DELETE FROM system_name_registry WHERE base_name = ?", ("Venus",))
    finally:
        conn.close()

    counts = dedupe_names(mysql_config)
    assert counts["star_systems"] == 1
    assert counts["sectors"] == 0

    conn = _db.get_connection(mysql_config)
    try:
        sector_name = conn.execute("SELECT name FROM sectors WHERE id = ?", (sector_id,)).fetchone()["name"]
        system_name = conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"]
    finally:
        conn.close()

    # The sector is never decorated; only the system is.
    assert sector_name == "Venus"
    assert system_name != "Venus"
    assert system_name.endswith("Venus")

    counts_again = dedupe_names(mysql_config)
    assert sum(counts_again.values()) == 0


def test_dedupe_is_a_true_no_op_on_a_database_with_no_collisions(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_sector(conn, SpaceSector(name="Mercury"))
        with conn:
            system = StarSystem(system_config=SystemConfig())
            system.star.name = "Uniquesys"
            _db.insert_star_system(conn, system, system.system_config)
    finally:
        conn.close()

    counts = dedupe_names(mysql_config)
    assert sum(counts.values()) == 0
