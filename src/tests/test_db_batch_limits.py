# tests/test_db_batch_limits.py

"""
TEST.17: batched writes at the limits. `Connection.flush` writes a batch
whose rows together pass the server's `max_allowed_packet` (read at run
time) as several statements instead of one the server refuses, and
refuses a single row too big to send with a clear error instead of a
dropped connection. Rows held back child-table-first are still written
parents first: an ordinary parent table, a table that refers to itself
(`nebulae.inside_nebula_id`), and the FK cycle `_table_ranks` mentions
(star systems, stars, black holes, supernova remnants), whose same-rank
rows keep their insertion order.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import pymysql
import pytest

from stellarObjects import _db
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.generation.system import StarSystem

TEXT_MAX = 65535
"""int: Bytes a TEXT column holds."""


def _seeded(config):
    """Saves a sector with a planet-bearing system, a nebula and a
    supernova remnant with a black hole core, and returns one row of
    each table the tests copy (`{table: row}`)."""
    sector = SpaceSector("Template Sector", edge_ly=11.5)
    for _ in range(30):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        system = StarSystem(system_config=cfg)
        if any(p.body_type != "a" for p in system.planets):
            break
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
    sector.add_phenomenon(Nebula(SystemConfig(), name="Template Nebula"), "nebula", position=(1.0, 1.0, 1.0))
    for _ in range(200):
        remnant = SupernovaRemnant(SystemConfig(), name="Template Remnant")
        if remnant.compact_remnant is not None and type(remnant.compact_remnant).__name__ == "BlackHole":
            break
    else:
        pytest.fail("no supernova remnant with a black hole core")
    sector.add_phenomenon(remnant, "supernova-remnant", position=(2.0, 2.0, 2.0))
    _db.save_sector(sector, config=config)
    conn = _db.get_connection(config)
    try:
        return {table: conn.execute(f"SELECT * FROM {table} ORDER BY id LIMIT 1").fetchone()
                for table in ("star_systems", "stars", "planets", "nebulae", "black_holes", "supernova_remnants")}
    finally:
        conn.close()


def _insert(conn, template, table, **changes):
    """A plain INSERT of a copy of `template` with `changes` -- held back
    by `batched`, its id from `id_blocks`."""
    values = {column: value for column, value in template[table].items() if column != "id"}
    values.update(changes)
    columns = list(values)
    return conn.execute(
        f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({', '.join('?' * len(columns))})",
        tuple(values.values()),
    ).lastrowid


def _max_allowed_packet(conn):
    return conn.execute("SELECT @@max_allowed_packet AS n").fetchone()["n"]


def test_a_batch_bigger_than_max_allowed_packet_is_split(mysql_config):
    template = _seeded(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        packet = _max_allowed_packet(conn)
        text = "x" * TEXT_MAX
        # Four full TEXT columns a row: comfortably past the packet within
        # one `_BATCH_ROWS` statement on a 16 MB (MariaDB) or 64 MB (MySQL 8) server.
        count = packet // (4 * TEXT_MAX) + 2
        if count > _db._BATCH_ROWS:
            pytest.skip(f"max_allowed_packet {packet} is too large to pass in one {_db._BATCH_ROWS}-row statement")
        system_id = template["star_systems"]["id"]
        with conn:
            with conn.batched():
                ids = [
                    _insert(conn, template, "planets", name=f"Bulk {n}", orbital_index=100 + n, description=text,
                            atmosphere=text, composition=text, flavor_text=text)
                    for n in range(count)
                ]
        rows = conn.execute(
            "SELECT id, star_system_id, LENGTH(description) + LENGTH(atmosphere) + LENGTH(composition)"
            " + LENGTH(flavor_text) AS size FROM planets WHERE name LIKE ? ORDER BY id", ("Bulk %",)
        ).fetchall()
        assert [r["id"] for r in rows] == ids
        assert {r["size"] for r in rows} == {4 * TEXT_MAX}
        assert {r["star_system_id"] for r in rows} == {system_id}
        assert sum(r["size"] for r in rows) > packet
    finally:
        conn.close()


def test_a_single_row_too_big_to_send_is_a_clear_error(mysql_config, monkeypatch):
    template = _seeded(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        # Pretend the server's limit is tiny; a real row that big can't
        # fit in a TEXT column.
        monkeypatch.setattr(_db, "_max_packet", lambda connection: 64 * 1024)
        with pytest.raises(pymysql.err.OperationalError) as caught:
            with conn:
                with conn.batched():
                    _insert(conn, template, "planets", name="Small", orbital_index=200)
                    _insert(conn, template, "planets", name="Huge", orbital_index=201, description="y" * TEXT_MAX)
        assert caught.value.args[0] == 1153
        assert "planets" in caught.value.args[1] and "max_allowed_packet" in caught.value.args[1]
        # The connection is still there, and nothing of the batch was kept.
        assert conn.execute("SELECT COUNT(*) AS n FROM planets WHERE name IN ('Small', 'Huge')").fetchone()["n"] == 0
    finally:
        conn.close()


def test_child_rows_held_first_are_written_after_their_parents(mysql_config, monkeypatch):
    template = _seeded(mysql_config)
    monkeypatch.setattr(_db, "_BATCH_ROWS", 2)
    conn = _db.get_connection(mysql_config)
    try:
        old_system = template["star_systems"]["id"]
        with conn:
            with conn.batched():
                # The planets statement is the batch's first, but its later
                # rows point at systems held after it.
                planets = [_insert(conn, template, "planets", name="Early Planet", orbital_index=300,
                                   star_system_id=old_system, star_id=None)]
                for n in range(3):
                    system = _insert(conn, template, "star_systems", name=f"Late System {n}")
                    planets.append(_insert(conn, template, "planets", name=f"Late Planet {n}", orbital_index=301 + n,
                                           star_system_id=system, star_id=None))
        rows = conn.execute(
            "SELECT p.id, s.name FROM planets p JOIN star_systems s ON s.id = p.star_system_id"
            " WHERE p.name LIKE ? ORDER BY p.id", ("% Planet%",)
        ).fetchall()
        assert [r["id"] for r in rows] == planets
        assert [r["name"] for r in rows] == [template["star_systems"]["name"]] + [f"Late System {n}" for n in range(3)]
    finally:
        conn.close()


def test_self_referencing_rows_keep_their_order_across_statements(mysql_config, monkeypatch):
    template = _seeded(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        self_referencing = {r["t"].lower() for r in conn.execute(
            "SELECT TABLE_NAME AS t FROM information_schema.KEY_COLUMN_USAGE"
            " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = REFERENCED_TABLE_NAME").fetchall()}
        assert self_referencing == {"nebulae"}
        monkeypatch.setattr(_db, "_BATCH_ROWS", 2)
        with conn:
            with conn.batched():
                chain = [_insert(conn, template, "nebulae", name="Nest 0", inside_nebula_id=None)]
                for n in range(1, 5):
                    chain.append(_insert(conn, template, "nebulae", name=f"Nest {n}", inside_nebula_id=chain[-1]))
        rows = conn.execute("SELECT id, inside_nebula_id FROM nebulae WHERE name LIKE ? ORDER BY id",
                            ("Nest %",)).fetchall()
        assert [(r["id"], r["inside_nebula_id"]) for r in rows] == list(zip(chain, [None] + chain[:-1]))
    finally:
        conn.close()


def test_the_containment_cycle_tables_share_a_rank():
    ranks = _db._condensed_ranks({
        "star_systems": {"supernova_remnants"}, "stars": {"star_systems"}, "black_holes": {"stars"},
        "supernova_remnants": {"black_holes"}, "planets": {"star_systems", "stars"},
    })
    assert len({ranks[t] for t in ("star_systems", "stars", "black_holes", "supernova_remnants")}) == 1
    assert ranks["planets"] > ranks["stars"]


def test_cycle_rows_are_written_in_insertion_order(mysql_config):
    template = _seeded(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        ranks = _db._table_ranks(conn._conn, mysql_config._key())
        cycle = ("star_systems", "stars", "black_holes", "supernova_remnants")
        assert len({ranks[t] for t in cycle}) == 1
        old_system = template["star_systems"]["id"]
        old_remnant = template["supernova_remnants"]["id"]
        with conn:
            with conn.batched():
                # The stars statement is held first, for an existing
                # system; the next rows each point at the one before.
                stars = [_insert(conn, template, "stars", star_system_id=old_system, name="Old Companion",
                                 role="secondary")]
                system = _insert(conn, template, "star_systems", name="Cycle System", inside_remnant_id=old_remnant)
                stars.append(_insert(conn, template, "stars", star_system_id=system, name="Cycle System"))
                hole = _insert(conn, template, "black_holes", star_id=stars[-1], name="Cycle Hole")
                remnant = _insert(conn, template, "supernova_remnants", name="Cycle Remnant",
                                  compact_remnant_black_hole_id=hole)
                inner = _insert(conn, template, "star_systems", name="Inner System", inside_remnant_id=remnant)
        assert conn.execute("SELECT star_system_id FROM stars WHERE id = ?", (stars[1],)).fetchone()["star_system_id"] \
            == system
        assert conn.execute("SELECT star_id FROM black_holes WHERE id = ?", (hole,)).fetchone()["star_id"] == stars[1]
        assert conn.execute("SELECT compact_remnant_black_hole_id AS h FROM supernova_remnants WHERE id = ?",
                            (remnant,)).fetchone()["h"] == hole
        assert conn.execute("SELECT inside_remnant_id AS r FROM star_systems WHERE id = ?",
                            (inner,)).fetchone()["r"] == remnant
        assert conn.execute("SELECT COUNT(*) AS n FROM stars WHERE id IN (?, ?)", tuple(stars)).fetchone()["n"] == 2
    finally:
        conn.close()
