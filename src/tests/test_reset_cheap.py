"""OPS.254: a reset of a galaxy with nothing in it is cheap. A table that holds no row and whose id counter was never used has nothing to clear, so the wipe leaves it alone (Boss's New galaxy run on an
empty database spent 36 s truncating 111 tables of which one held a row); a table that is empty but whose counter has moved is still truncated, so the next galaxy numbers from 1 as after any reset."""

import pytest

from planetgen.cli import reset
from planetgen.db import store

pytestmark = pytest.mark.fresh_database


def _truncated(monkeypatch):
    seen = []
    real = reset.longsql.run

    def run(conn, sql, label, **kwargs):
        seen.append(sql)
        return real(conn, sql, label, **kwargs)

    monkeypatch.setattr(reset.longsql, "run", run)
    return seen


def test_a_wipe_of_a_galaxy_with_nothing_in_it_truncates_nothing(mysql_config, monkeypatch):
    store.get_connection(mysql_config).close()          # the schema, no rows
    seen = _truncated(monkeypatch)
    wiped = reset.wipe(mysql_config)
    assert seen == [], f"tables truncated with nothing to clear: {seen[:5]}"
    assert wiped.with_rows == 0 and len(wiped.tables) > 100


def test_a_table_with_rows_and_a_table_whose_counter_moved_are_truncated_and_no_other(mysql_config, monkeypatch):
    conn = store.get_connection(mysql_config)
    try:
        counter = next(row["name"] for row in conn.execute(
            "SELECT table_name AS name FROM information_schema.tables WHERE table_schema = ? AND auto_increment IS NOT NULL ORDER BY 1", (mysql_config.database,)).fetchall())
        columns = [row["name"] for row in conn.execute(
            "SELECT column_name AS name FROM information_schema.columns WHERE table_schema = ? AND table_name = ? AND extra LIKE '%%auto_increment%%'",
            (mysql_config.database, counter)).fetchall()]
        conn.execute(f"ALTER TABLE {counter} AUTO_INCREMENT = 500")      # a counter that has moved, with no row left
        conn.commit()
    finally:
        conn.close()
    seen = _truncated(monkeypatch)
    reset.wipe(mysql_config)
    assert seen == [f"TRUNCATE TABLE {counter}"], seen
    check = store.get_connection(mysql_config)
    try:
        assert int(check.execute("SELECT auto_increment AS n FROM information_schema.tables WHERE table_schema = ? AND table_name = ?",
                                 (mysql_config.database, counter)).fetchone()["n"] or 1) == 1 and columns
    finally:
        check.close()


class _Spy:
    """A connection that says every statement it runs."""

    def __init__(self, conn):
        self._conn, self.sql = conn, []

    def execute(self, sql, *args, **kwargs):
        self.sql.append(sql)
        return self._conn.execute(sql, *args, **kwargs)

    def __getattr__(self, name):
        return getattr(self._conn, name)


def test_clearing_the_bright_stars_of_a_new_galaxy_only_reads(mysql_config):
    from tests.test_galaxy_gen import _seed_skeleton
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1)])
    conn = store.get_connection(mysql_config)
    try:
        spy = _Spy(conn)
        store.clear_bright_stars(spy)
        written = [sql for sql in spy.sql if sql.lstrip().split(None, 1)[0].upper() in ("TRUNCATE", "DELETE", "UPDATE", "INSERT", "ALTER")]
        assert written == [], f"a new galaxy has no bright stars to clear, and the clear wrote: {written}"
    finally:
        conn.close()


def test_clearing_the_bright_stars_of_a_scatter_that_was_drawn_clears_it_all(mysql_config):
    from tests.test_galaxy_gen import _seed_skeleton
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1)])
    conn = store.get_connection(mysql_config)
    try:
        store.record_bright_star_scatter(conn, 1.5, 7, 3.0)
        conn.commit()
        assert not store.bright_stars_hold_nothing(conn)
        spy = _Spy(conn)
        store.clear_bright_stars(spy)
        assert any(sql.startswith("TRUNCATE TABLE bright_stars") for sql in spy.sql) and spy.sql[-1].startswith("TRUNCATE TABLE bright_stars"), "the TRUNCATE comes last"
        row = conn.execute("SELECT bright_star_min_luminosity_sol AS a, bright_star_seed AS b, bright_star_mass_limit_sol AS c FROM galaxy_shape").fetchone()
        assert row["a"] is None and row["b"] is None and row["c"] is None
        assert store.bright_stars_hold_nothing(conn)
    finally:
        conn.close()
