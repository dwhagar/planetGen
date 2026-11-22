# tests/test_fingerprint.py

"""
GEN.58: `planetgen fingerprint` (`planetgen.db.fingerprint`) gives two
builds of one galaxy seed the same digests whatever their row ids, a
different seed different ones, and a change to one sector's content
changes only that sector's digest; clocks, location text and nearest
links don't count.
"""

import contextlib
import io
import re
import sys


from planetgen.cli import generate as generate_cli
from planetgen.db import fingerprint, store
from planetgen.generation import run_common

from tests.test_galaxy_gen import _mysql_argv, _plan_wide_galaxy
from tests.test_parallel_generation import GALAXY_SEED, make_database  # noqa: F401 -- a fixture

_SCHEMA = store.os.path.join(store.os.path.dirname(store.__file__), "schema.sql")


def test_every_table_is_compared_or_left_out_on_purpose():
    with open(_SCHEMA, encoding="utf-8") as file:
        tables = set(re.findall(r"CREATE TABLE IF NOT EXISTS (\w+)", file.read()))
    assert not fingerprint.CONTENT_TABLES & fingerprint.LEFT_OUT_TABLES
    assert tables == fingerprint.CONTENT_TABLES | fingerprint.LEFT_OUT_TABLES


def test_every_content_table_is_in_one_section():
    assert set(fingerprint.SECTIONS["plan"]) == set(fingerprint.PLAN_TABLES)
    tables = [table for tables in fingerprint.SECTIONS.values() for table in tables]
    assert len(tables) == len(set(tables)) and set(tables) == fingerprint.CONTENT_TABLES


RAW_ID_NOT_A_ROW = {(table, column) for table in ("asteroid_belts", "black_holes", "comets", "interstellar_comets", "moons", "neutron_stars", "planets", "rogue_planets", "stars") for column in ("parent_id",)} | {
    ("asteroid_belts", "host_node_id"), ("planets", "host_node_id"), ("stars", "parent_node_id"), ("asteroid_fields", "origin_event_id"), ("sectors", "max_mass_id"), ("sectors", "override_id")}
"""set: The `*_id` columns of the content tables that are no foreign key and need no key: `parent_id` is keyed by `parent_type`, the node numbers are of the row's own system, the rest are not row ids of a content table."""


def test_a_column_that_holds_a_row_id_is_keyed_not_read_raw(make_database):
    """GEN.1827: a column that names another row of a content table (no declared foreign key) is in `ROW_ID_COLUMNS`, so the fingerprint reads that row's key and not the id it got in this run."""
    conn = store.get_connection(make_database())
    try:
        reader = fingerprint._Reader(conn)
        raw_ids = {(table, column) for table in fingerprint.CONTENT_TABLES for column in reader.columns(table)
                   if column.endswith("_id") and column not in reader.foreign_keys(table)}
        assert fingerprint.ROW_ID_COLUMNS.keys() <= raw_ids
        assert raw_ids - fingerprint.ROW_ID_COLUMNS.keys() == RAW_ID_NOT_A_ROW, "a new column holds a row id: add it to ROW_ID_COLUMNS, or here if it is no row id"
        assert ("planets", "lagrange_host_id") in fingerprint.ROW_ID_COLUMNS
        assert all(parent in fingerprint.CONTENT_TABLES for parent in fingerprint.ROW_ID_COLUMNS.values())
    finally:
        conn.close()


def _generate(config, monkeypatch, seed=GALAXY_SEED, shift_ids=False, limit=None):
    _plan_wide_galaxy(config, galaxy_seed=seed)
    if shift_ids:
        conn = store.get_connection(config)
        try:
            for table in sorted(store.ID_BLOCK_TABLES):
                conn.execute("INSERT INTO id_blocks (table_name, next_id) VALUES (?, 5000)"
                             " ON DUPLICATE KEY UPDATE next_id = 5000", (table,))
                conn.execute(f"ALTER TABLE {table} AUTO_INCREMENT = 5000")
            conn.commit()
        finally:
            conn.close()
        store.forget_id_blocks(config._key())
    monkeypatch.setenv("PLANETGEN_CONTROL_DATABASE", config.database)
    for counter in run_common.RUN_COUNTS:
        run_common.RUN_COUNTS[counter] = 0
    _cli(["sectors", "--shell", "1", "--num-systems", "4", "--workers", "1", *(["--limit", str(limit)] if limit else [])], config)


def _cli(argv, config):
    old_argv = sys.argv
    out = io.StringIO()
    try:
        sys.argv = ["planetgen"] + argv + _mysql_argv(config)
        with contextlib.redirect_stdout(out):
            generate_cli.main()
    finally:
        sys.argv = old_argv
    return out.getvalue()


def _print(config, *options):
    return [line for line in _cli(["fingerprint", *options], config).splitlines()
            if re.match(r"^(\d+ -?\d+ \d+|plan|region|unplaced) ", line)]


def _region(config, **kwargs):
    conn = store.get_connection(config)
    try:
        return fingerprint.region_fingerprint(conn, **kwargs)
    finally:
        conn.close()


SEED_SECTORS = 8
"""int: The sectors of the shell each galaxy of the another-seed test generates (it asserts at least eight sectors in common)."""

ROW_ID_SECTORS = 2
"""int: The sectors of the shell the row-id test generates (a sector takes about a minute, with its scatter cells; the second galaxy starts its ids at 5000 so none of them can match the first's)."""


def test_one_seed_gives_one_fingerprint_whatever_the_row_ids(make_database, monkeypatch):
    first, second = make_database(), make_database()
    _generate(first, monkeypatch, limit=ROW_ID_SECTORS)        # TEST.751: two sectors of the shell are enough to carry every table's ids; the whole shell took 12 minutes a galaxy
    _generate(second, monkeypatch, shift_ids=True, limit=ROW_ID_SECTORS)
    conn = store.get_connection(second)
    try:
        assert conn.execute("SELECT MIN(id) AS n FROM star_systems").fetchone()["n"] >= 5000
    finally:
        conn.close()
    a, b = _region(first), _region(second)
    assert len(a.sectors) >= ROW_ID_SECTORS and a.plan is not None
    assert a == b
    assert _print(first) == _print(second) == a.lines()


def test_another_seed_gives_another_fingerprint(make_database, monkeypatch):
    first, second = make_database(), make_database()
    _generate(first, monkeypatch, limit=SEED_SECTORS)        # TEST.417: the first SEED_SECTORS sectors of the shell are enough for the eight common ones this asserts; the whole shell took 12 minutes a galaxy
    _generate(second, monkeypatch, seed=bytes([GALAXY_SEED[0] ^ 0x80]) + GALAXY_SEED[1:], limit=SEED_SECTORS)
    assert _region(first).region != _region(second).region
    a, b = _region(first, shells=[1]), _region(second, shells=[1])
    common = dict(a.sectors).keys() & dict(b.sectors).keys()
    assert len(common) >= 8  # the shell's sectors; the scatter's cells off them differ too
    assert all(dict(a.sectors)[label] != dict(b.sectors)[label] for label in common)


def test_a_change_moves_only_its_own_sectors_digest(make_database, monkeypatch):
    config = make_database()
    _generate(config, monkeypatch, limit=ROW_ID_SECTORS)        # TEST.417: two sectors show that one changed sector moves only its own digest; the whole shell took 12 minutes
    before = _region(config)
    conn = store.get_connection(config)
    try:
        # Clocks, the location text and the nearest links don't count.
        conn.execute("UPDATE planets SET epoch_unix = 12345, next_update_due = 67890")
        conn.execute("UPDATE star_systems SET location = 'somewhere else'")
        conn.execute("DELETE FROM nearest_systems")
        conn.commit()
        assert _region(config) == before
        planet = conn.execute(
            "SELECT p.id, s.shell_index, s.band_index, s.wedge_index FROM planets p"
            " JOIN star_systems ss ON ss.id = p.star_system_id JOIN sectors s ON s.id = ss.sector_id"
            " ORDER BY p.id LIMIT 1").fetchone()
        conn.execute("UPDATE planets SET mass_kg = mass_kg * 1.000001 WHERE id = ?", (planet["id"],))
        conn.commit()
    finally:
        conn.close()
    after = _region(config)
    changed = [label for (label, old), (_label, new) in zip(before.sectors, after.sectors) if old != new]
    assert changed == [f"{planet['shell_index']} {planet['band_index']} {planet['wedge_index']}"]
    assert after.region != before.region and after.plan == before.plan
    one = _region(config, addresses=[(planet["shell_index"], planet["band_index"], planet["wedge_index"])])
    assert one.plan is None and one.sectors == [(changed[0], dict(after.sectors)[changed[0]])]
    assert _print(config, "--sector", *changed[0].split())[0] == f"{changed[0]} {dict(after.sectors)[changed[0]]}"


def test_a_changed_number_moves_only_its_sections_float_digest(make_database, monkeypatch):
    """GEN.304: the structural digest of a section holds when a number changes; the float digest moves."""
    config = make_database()
    _generate(config, monkeypatch, limit=1)        # TEST.417: one planet in one sector is all this changes
    conn = store.get_connection(config)
    try:
        planet = conn.execute(
            "SELECT p.id, s.id AS sector_id FROM planets p JOIN star_systems ss ON ss.id = p.star_system_id"
            " JOIN sectors s ON s.id = ss.sector_id ORDER BY p.id LIMIT 1").fetchone()
        before = fingerprint.sector_sections(conn, sector_id=planet["sector_id"])
        conn.execute("UPDATE planets SET mass_kg = mass_kg * 1.000001 WHERE id = ?", (planet["id"],))
        conn.commit()
        after = fingerprint.sector_sections(conn, sector_id=planet["sector_id"])
        assert {name for name in before if before[name] != after[name]} == {"planet"}
        assert before["planet"].structural == after["planet"].structural
        assert before["planet"].floats != after["planet"].floats and before["planet"].full != after["planet"].full
        conn.execute("UPDATE planets SET name = 'Renamed' WHERE id = ?", (planet["id"],))
        conn.commit()
        renamed = fingerprint.sector_sections(conn, sector_id=planet["sector_id"])
        assert renamed["planet"].structural != after["planet"].structural
    finally:
        conn.close()
