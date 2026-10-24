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


def _generate(config, monkeypatch, seed=GALAXY_SEED, shift_ids=False):
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
    _cli(["sectors", "--shell", "1", "--num-systems", "4", "--workers", "1"], config)


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


def test_one_seed_gives_one_fingerprint_whatever_the_row_ids(make_database, monkeypatch):
    first, second = make_database(), make_database()
    _generate(first, monkeypatch)
    _generate(second, monkeypatch, shift_ids=True)
    conn = store.get_connection(second)
    try:
        assert conn.execute("SELECT MIN(id) AS n FROM star_systems").fetchone()["n"] >= 5000
    finally:
        conn.close()
    a, b = _region(first), _region(second)
    assert len(a.sectors) >= 8 and a.plan is not None
    assert a == b
    assert _print(first) == _print(second) == a.lines()


def test_another_seed_gives_another_fingerprint(make_database, monkeypatch):
    first, second = make_database(), make_database()
    _generate(first, monkeypatch)
    _generate(second, monkeypatch, seed=bytes([GALAXY_SEED[0] ^ 0x80]) + GALAXY_SEED[1:])
    assert _region(first).region != _region(second).region
    a, b = _region(first, shells=[1]), _region(second, shells=[1])
    common = dict(a.sectors).keys() & dict(b.sectors).keys()
    assert len(common) >= 8  # the shell's sectors; the scatter's cells off them differ too
    assert all(dict(a.sectors)[label] != dict(b.sectors)[label] for label in common)


def test_a_change_moves_only_its_own_sectors_digest(make_database, monkeypatch):
    config = make_database()
    _generate(config, monkeypatch)
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
    _generate(config, monkeypatch)
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
