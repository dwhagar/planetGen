"""PERF.147: the pieces of the population pass that work through many rows (`civilizations`, `technology`) say where they are, with the real counts, and write the rows they write without it."""

import pytest

from planetgen.db import store
from planetgen.population import model
from tests.population_support import every_window_has_a_civilization
from tests.test_engine_kill_population import _stored_galaxy

pytestmark = [pytest.mark.fresh_database]


class Counted:
    """What the engine's context gives a piece to say where it is: `step(name, total)` and `advance(amount)`."""

    def __init__(self):
        self.steps, self.done = [], 0

    def step(self, name, total=None):
        self.steps.append((name, total))

    def advance(self, amount=1):
        self.done += amount


def _rows(conn):
    species = [tuple(row.values()) for row in conn.execute("SELECT * FROM species ORDER BY id").fetchall()]
    polities = [tuple(value for key, value in row.items() if key != "id") for row in conn.execute("SELECT * FROM polities ORDER BY species_id").fetchall()]
    interests = [tuple(row.values()) for row in conn.execute("SELECT star_system_id, reason, best_score FROM system_interests ORDER BY star_system_id, reason").fetchall()]
    return species, polities, interests


def _populated(mysql_config):
    _stored_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    with every_window_has_a_civilization():
        for piece in model.population_pieces(conn):
            model.do_population_piece(conn, piece)
            conn.commit()
    conn.execute("INSERT INTO system_interests (star_system_id, reason, n, best_score, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc)"
                 " SELECT id, 'civilization', 1, 0, 0, 0, 0 FROM star_systems WHERE id IN (SELECT star_system_id FROM species)")
    conn.commit()
    return conn


def _spoil(conn, piece):
    """What the piece writes, as an earlier pass left it out of date, so the piece has something to write."""
    if piece == ["civilizations"]:
        conn.execute("UPDATE species SET era = 'stale', spacefaring = 0")
        conn.execute("UPDATE polities SET reach_ly = reach_ly + 1")
    else:
        conn.execute("UPDATE system_interests SET best_score = -5 WHERE reason = 'civilization'")
    conn.commit()


@pytest.mark.parametrize("piece", [["civilizations"], ["technology"]], ids=lambda piece: piece[0])
def test_a_piece_that_counts_writes_the_rows_it_writes_without_counting(mysql_config, piece, monkeypatch):
    monkeypatch.setattr(model, "COUNT_BATCH_ROWS", 1)       # every row is a batch of its own
    conn = _populated(mysql_config)
    try:
        _spoil(conn, piece)
        model.do_population_piece(conn, piece)
        conn.commit()
        plain = _rows(conn)
        _spoil(conn, piece)
        counted = Counted()
        model.do_population_piece(conn, piece, counted)
        conn.commit()
        assert _rows(conn) == plain, "the same rows"
        assert counted.steps, "the piece began a count"
        by_name = {}
        for name, total in counted.steps:
            by_name[name] = total
        assert all(total for total in by_name.values()), f"every count has its total: {counted.steps}"
        assert counted.done >= 1
    finally:
        conn.close()


def test_the_technology_piece_counts_the_system_ids_it_has_gone_through(mysql_config, monkeypatch):
    from planetgen.db import interests
    monkeypatch.setattr(interests, "TECH_BATCH_SYSTEMS", 1)
    conn = _populated(mysql_config)
    try:
        counted = Counted()
        model.do_population_piece(conn, ["technology"], counted)
        [(name, total)] = counted.steps
        assert total > 0 and counted.done == total, f"{name}: the ids covered are the ids there are"
    finally:
        conn.close()
