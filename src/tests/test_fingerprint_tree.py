# tests/test_fingerprint_tree.py

"""
GEN.305: a sector's leaf digest is taken from the sector in memory when it is saved, and the band and shell nodes and
the galaxy root above it (`planetgen.db.fingerprint_tree`) are a pure function of the leaves: changing one sector
changes exactly its band node, its shell node and the root; a region's digest changes when a sector is added; and
none of it depends on the order sectors were saved in or on the number of workers.
"""

import contextlib
import io
import re
import sys

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import fingerprint_tree as tree, store

from tests.test_galaxy_gen import _mysql_argv, _plan_wide_galaxy
from tests.test_parallel_generation import GALAXY_SEED, make_database  # noqa: F401 -- a fixture


def _digest(number):
    return bytes([number]) * 32


def test_a_node_depends_on_every_leaf_and_their_order():
    leaves = [(0, _digest(1)), (1, _digest(2)), (2, _digest(3))]
    band, count = tree._band(1, 0, leaves)
    assert count == 3
    assert tree._band(1, 0, leaves)[0] == band
    assert tree._band(1, 0, leaves[:2])[0] != band
    assert tree._band(1, 0, [(0, _digest(1)), (1, _digest(2)), (2, _digest(4))])[0] != band
    assert tree._band(1, 0, [(0, _digest(2)), (1, _digest(1)), (2, _digest(3))])[0] != band
    assert tree._band(1, 1, leaves)[0] != band  # the band is part of the node
    shell, shell_count = tree._shell(1, [(0, band, 3), (1, band, 3)])
    assert shell_count == 6 and shell != tree._shell(1, [(0, band, 3)])[0]
    assert tree._root([(1, shell, 6)])[0] != tree._root([(1, shell, 6), (2, shell, 6)])[0]


def test_the_collision_decoration_does_not_count():
    assert tree.plain({"name": "Beta Vestara", "xs": ["Little Omicron Kaan", 1.5]}) == {
        "name": "Vestara", "xs": ["Kaan", 1.5]}


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


def _generate(config, monkeypatch, workers=1, limit=None):
    _plan_wide_galaxy(config, galaxy_seed=GALAXY_SEED)
    monkeypatch.setenv("PLANETGEN_CONTROL_DATABASE", config.database)
    argv = ["galaxy", "--shell", "1", "--num-systems", "4", "--workers", str(workers)]
    if limit:
        argv += ["--limit", str(limit)]
    _cli(argv, config)


def _connection(config):
    return store.get_connection(config)


def test_every_saved_sector_has_a_leaf_and_the_stored_nodes_are_the_tree(make_database, monkeypatch):
    config = make_database()
    _generate(config, monkeypatch)
    conn = _connection(config)
    try:
        sectors = conn.execute("SELECT COUNT(*) AS n FROM sectors WHERE shell_index IS NOT NULL").fetchone()["n"]
        assert sectors >= 8
        assert conn.execute("SELECT COUNT(*) AS n FROM sectors WHERE leaf_digest IS NULL"
                            " AND shell_index IS NOT NULL").fetchone()["n"] == 0
        assert conn.execute("SELECT COUNT(*) AS n FROM sector_tree_nodes WHERE dirty = 1").fetchone()["n"] > 0
        root, count = tree.galaxy_root(conn)
        conn.commit()
        assert count == sectors and root
        assert conn.execute("SELECT COUNT(*) AS n FROM sector_tree_nodes WHERE dirty = 1").fetchone()["n"] == 0
        found = tree.tree(conn)
        assert found["root"] == (root, count)
        for (shell, band), (digest, n) in found["bands"].items():
            row = conn.execute("SELECT digest, child_count FROM sector_tree_nodes WHERE level = 1 AND shell_index = ?"
                               " AND band_index = ?", (shell, band)).fetchone()
            assert (bytes(row["digest"]).hex(), row["child_count"]) == (digest, n)
        for shell, (digest, n) in found["shells"].items():
            row = conn.execute("SELECT digest, child_count FROM sector_tree_nodes WHERE level = 2 AND shell_index = ?",
                               (shell,)).fetchone()
            assert (bytes(row["digest"]).hex(), row["child_count"]) == (digest, n)
    finally:
        conn.close()
    printed = [line for line in _cli(["fingerprint", "--stored"], config).splitlines()
               if re.match(r"^(band|shell|root|stored-region) ", line)]
    assert printed[-2].startswith(f"root {root} ({count} sector(s))") and printed[-1].startswith("stored-region ")


def test_one_seed_gives_one_tree_whatever_the_workers(make_database, monkeypatch):
    results = []
    for workers in (1, 4):
        config = make_database()
        _generate(config, monkeypatch, workers=workers)
        conn = _connection(config)
        try:
            results.append(tree.tree(conn))
        finally:
            conn.close()
    assert results[0]["root"] is not None
    assert results[0] == results[1]


def test_changing_one_sector_changes_its_band_its_shell_and_the_root(make_database, monkeypatch):
    config = make_database()
    _generate(config, monkeypatch)
    conn = _connection(config)
    try:
        tree.refresh(conn)
        before = tree.tree(conn)
        two_bands = [address for address in before["bands"]]
        assert len(two_bands) >= 2
        row = conn.execute("SELECT id, shell_index, band_index, wedge_index FROM sectors WHERE shell_index IS NOT NULL"
                           " ORDER BY shell_index, band_index, wedge_index LIMIT 1").fetchone()
        address = (row["shell_index"], row["band_index"], row["wedge_index"])
        stored_before = conn.execute("SELECT level, band_index, digest FROM sector_tree_nodes ORDER BY level, band_index"
                                     ).fetchall()
        tree.store_leaf(conn, row["id"], address, _digest(9))
        conn.commit()
        assert tree.refresh(conn) == 2  # its band node and its shell node
        after = tree.tree(conn)
        changed_bands = [key for key in before["bands"] if before["bands"][key] != after["bands"][key]]
        changed_shells = [key for key in before["shells"] if before["shells"][key] != after["shells"][key]]
        assert changed_bands == [address[:2]] and changed_shells == [address[0]]
        assert before["root"] != after["root"] and before["root"][1] == after["root"][1]
        stored_after = conn.execute("SELECT level, band_index, digest FROM sector_tree_nodes ORDER BY level, band_index"
                                    ).fetchall()
        assert sum(a["digest"] != b["digest"] for a, b in zip(stored_before, stored_after)) == 2
        assert tree.galaxy_root(conn)[0] == after["root"][0]
    finally:
        conn.close()


def test_a_region_digest_changes_when_a_sector_is_added(make_database, monkeypatch):
    config = make_database()
    _generate(config, monkeypatch)
    conn = _connection(config)
    try:
        whole, count = tree.region_digest(conn)
        row = conn.execute("SELECT id, shell_index, band_index, wedge_index, leaf_digest FROM sectors"
                           " WHERE shell_index IS NOT NULL ORDER BY id LIMIT 1").fetchone()
        conn.execute("UPDATE sectors SET leaf_digest = NULL WHERE id = ?", (row["id"],))
        smaller, fewer = tree.region_digest(conn)
        assert fewer == count - 1 and smaller != whole
        conn.execute("UPDATE sectors SET leaf_digest = ? WHERE id = ?", (bytes(row["leaf_digest"]), row["id"]))
        assert tree.region_digest(conn) == (whole, count)
        one = (row["shell_index"], row["band_index"], row["wedge_index"])
        assert tree.region_digest(conn, addresses=[one])[1] == 1
    finally:
        conn.close()
