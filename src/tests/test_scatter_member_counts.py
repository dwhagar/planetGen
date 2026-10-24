# tests/test_scatter_member_counts.py

"""
PERF.158: `scatter_systems.store_member_trees` sets the member counts and system masses of a batch of rows with one
statement, and the rows hold the same bits as the ones `store_member_tree` sets a row at a time.
"""

import pytest

from planetgen.db import scatter_systems, store
from planetgen.generation import multistar
from planetgen.util import draw
from tests.test_galactic_motion import SCATTER_RADIUS_PC, _insert_scatter_point, _seed_scatter_grid

pytestmark = pytest.mark.fresh_database

EPOCH = 1_700_000_000.0


def _trees(count):
    return [multistar.draw_member_tree(multistar.Placed("star", 1.0 + 0.37 * index), 5e9, draw.Stream(index + 1),
                                       n_stars=2 + index % 3) for index in range(count)]


def test_a_batch_of_counts_is_one_statement_and_the_same_bits_as_one_row_at_a_time(mysql_config, monkeypatch):
    _seed_scatter_grid(mysql_config)
    trees = _trees(7)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            single = [_insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, float(index), 0.0))
                      for index in range(len(trees))]
            batch = [_insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, 100.0 + index, 0.0))
                     for index in range(len(trees))]
            for row_id, tree in zip(single, trees):
                scatter_systems.store_member_tree(conn, "bright_stars", row_id, tree, EPOCH)
            literal = []
            real = conn.execute_literal
            monkeypatch.setattr(conn, "execute_literal", lambda sql: literal.append(sql) or real(sql))
            assert scatter_systems.store_member_trees(conn, "bright_stars", list(zip(batch, trees)), EPOCH) == len(trees)
            read = {row["id"]: (row["member_count"], row["system_mass_solar"]) for row in conn.execute(
                "SELECT id, member_count, system_mass_solar FROM bright_stars WHERE id IN (%s)" % ",".join(
                    map(str, single + batch))).fetchall()}
    finally:
        conn.close()
    assert len(literal) == 1 and literal[0].startswith("UPDATE bright_stars")
    assert [read[row_id] for row_id in batch] == [read[row_id] for row_id in single]
    assert all(count and mass for count, mass in read.values())
