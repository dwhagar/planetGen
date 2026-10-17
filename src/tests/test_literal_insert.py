# tests/test_literal_insert.py

"""
The bulk inserts that build their own statement text (`store._insert_literal_rows`: `insert_bright_stars`,
`insert_phenomenon_scatter`; PERF.73, PERF.48): the rows read back as written, the statements stay under their
size cap, and an empty stream writes nothing.
"""

import pytest

from planetgen.db import store


def _row(number):
    return (number % 7, number % 5, number, number * 1000, -number * 1000 + 1, 17,
            "young" if number % 2 else "old", "G2" if number % 3 else "B", "V", 2e30 + number, 7e5 + 0.125 * number,
            5772.0, 3.828e26 * (1 + number), 4.6 + number / 1e6, 10.0, 1.0 + number / 7, None if number % 4 else 12.5,
            number, None if number % 2 else -0.15)


@pytest.fixture
def conn(mysql_config):
    connection = store.get_connection(mysql_config)
    yield connection
    connection.close()


def _stored(conn):
    columns = ", ".join(store.BRIGHT_STAR_COLUMNS)
    return [tuple(row.values()) for row in conn.execute(
        f"SELECT {columns} FROM bright_stars ORDER BY wedge_index").fetchall()]


def test_bright_stars_read_back_exactly_as_written(conn):
    rows = [_row(number) for number in range(1, 301)]
    assert store.insert_bright_stars(conn, rows, batch_size=100) == 300
    conn.commit()
    stored = _stored(conn)
    assert len(stored) == 300
    for expected, got in zip(rows, stored):
        assert got == expected


def test_a_long_run_is_cut_into_statements_under_the_size_cap(conn, monkeypatch):
    sent = []
    real = conn.execute_literal

    def record(sql):
        sent.append(len(sql))
        return real(sql)

    monkeypatch.setattr(conn, "execute_literal", record)
    monkeypatch.setattr(store, "LITERAL_STATEMENT_CHARS", 3000)
    rows = [_row(number) for number in range(1, 201)]
    assert store.insert_bright_stars(conn, rows, batch_size=10_000) == 200
    conn.commit()
    assert len(sent) > 5 and max(sent) < 3000 + 600  # one row over the cap at most
    assert len(_stored(conn)) == 200


def test_no_rows_write_nothing(conn):
    assert store.insert_bright_stars(conn, []) == 0
    assert store.insert_phenomenon_scatter(conn, iter(())) == 0
