"""A connection MySQL drops mid-unit (the Boss's job-20261107-032338-3037 log: `InterfaceError: (0, '')` ended stage 13) shows its real error and a sky unit goes on."""
import pymysql
import pytest

from planetgen.db import store
from planetgen.generation import run_galaxy, sky_charting


class DeadConnection:
    """What pymysql holds after the server dropped the socket: every command raises `InterfaceError(0, '')`."""

    def __init__(self):
        self.closed = 0

    def rollback(self):
        raise pymysql.err.InterfaceError(0, "")

    def close(self):
        self.closed += 1


def a_store_connection(raw):
    conn = object.__new__(store.Connection)
    conn._conn, conn._batch, conn._flagged, conn._after_commit, conn._parity_touched = raw, {}, set(), [], set()
    conn._release_txn_locks = lambda: None
    return conn


def test_a_rollback_on_a_dropped_connection_does_not_hide_the_error_being_cleaned_up():
    conn = a_store_connection(DeadConnection())
    conn.rollback()                                   # used to raise InterfaceError(0, '') over the real error
    with pytest.raises(pymysql.err.OperationalError, match="gone away"):
        try:
            raise pymysql.err.OperationalError(2006, "MySQL server has gone away")
        except BaseException:
            conn.rollback()
            raise


def test_a_rollback_still_raises_other_errors():
    class Broken(DeadConnection):
        def rollback(self):
            raise pymysql.err.OperationalError(1205, "lock wait timeout")

    with pytest.raises(pymysql.err.OperationalError):
        a_store_connection(Broken()).rollback()


def test_which_errors_say_the_connection_is_gone():
    assert store.is_lost_connection(pymysql.err.InterfaceError(0, ""))
    assert store.is_lost_connection(pymysql.err.OperationalError(2013, "Lost connection"))
    assert not store.is_lost_connection(pymysql.err.OperationalError(1213, "deadlock"))
    assert not store.is_lost_connection(ValueError("x"))


def test_a_sky_unit_goes_on_after_losing_its_connection(monkeypatch):
    opened, calls = [], []

    def get_connection(config=None, **options):
        conn = a_store_connection(DeadConnection())
        conn.close = lambda: opened.append("closed")
        opened.append("opened")
        return conn

    def chart(conn, center, limit, on_batch=None, progress=None):
        calls.append(center)
        if len(calls) == 1:
            raise pymysql.err.OperationalError(2013, "Lost connection to MySQL server during query")
        return {"tiers": [], "stars": 7, "seconds": 1.0, "cell": (0, 0, 0), "floors": []}

    monkeypatch.setattr(store, "get_connection", get_connection)
    monkeypatch.setattr(sky_charting, "chart", chart)
    monkeypatch.setattr(run_galaxy, "_mysql_config_of", lambda args: None)
    result = run_galaxy._sky_unit({"database": {}, "sky": {"center": [0, 0, 0], "limit": 7.5}}, lambda done, total: None)
    assert len(calls) == 2 and opened == ["opened", "closed", "opened", "closed"] and result["stars"] == 7


def test_a_sky_unit_gives_up_after_its_tries_and_other_errors_fail_at_once(monkeypatch):
    monkeypatch.setattr(store, "get_connection", lambda config=None, **options: a_store_connection(DeadConnection()))
    monkeypatch.setattr(run_galaxy, "_mysql_config_of", lambda args: None)
    calls = []

    def lost(conn, center, limit, on_batch=None, progress=None):
        calls.append(1)
        raise pymysql.err.OperationalError(2006, "gone away")

    monkeypatch.setattr(sky_charting, "chart", lost)
    unit = {"database": {}, "sky": {"center": [0, 0, 0], "limit": 7.5}}
    with pytest.raises(pymysql.err.OperationalError):
        run_galaxy._sky_unit(unit, lambda done, total: None)
    assert len(calls) == run_galaxy.SKY_UNIT_CONNECTION_TRIES
    calls.clear()

    def other(conn, center, limit, on_batch=None, progress=None):
        calls.append(1)
        raise ValueError("a bug")

    monkeypatch.setattr(sky_charting, "chart", other)
    with pytest.raises(ValueError):
        run_galaxy._sky_unit(unit, lambda done, total: None)
    assert len(calls) == 1


def facts(**changes):
    return {"uptime": 5000, "aborted_clients": 3, "aborted_connects": 0, "wait_timeout": 28800, "net_read_timeout": 30, "net_write_timeout": 60, "max_allowed_packet": 67108864, **changes}


def test_the_verdict_names_a_restart_an_idle_timeout_an_abort_or_nothing_it_can_see():
    from planetgen.db import connection_loss as loss
    before = facts()
    assert "restarted" in loss.verdict(before, facts(uptime=40), 100, 1800)                                  # the uptime fell: the server came back
    assert "restarted" in loss.verdict(before, facts(uptime=5100), 100, 1800)                                # it grew by 100 s in 1800 s: it was up for only the last 100
    assert loss.verdict(before, facts(uptime=6800), 29000, 1800).startswith("idle timeout")                  # idle for more than wait_timeout
    assert "aborted the client" in loss.verdict(before, facts(uptime=6800, aborted_clients=4), 100, 1800)
    assert "something between" in loss.verdict(before, facts(uptime=6800), 100, 1800)
    assert "cannot tell" in loss.verdict(facts(uptime=None), facts(), 100, 1800)


def test_facts_reads_what_it_can_and_a_refused_status_gives_none():
    from planetgen.db import connection_loss as loss

    class Refuses:
        def execute(self, sql, params=()):
            raise pymysql.err.OperationalError(1227, "Access denied")

    assert loss.facts(Refuses()) == dict.fromkeys(("uptime", "aborted_clients", "aborted_connects", "wait_timeout", "net_read_timeout", "net_write_timeout", "max_allowed_packet"))


def test_a_sky_unit_that_loses_its_connection_says_why_on_the_new_one(monkeypatch):
    from planetgen.db import connection_loss as loss
    said, calls = [], []
    monkeypatch.setattr(loss, "facts", lambda conn: facts(uptime=40 if calls else 5000))
    monkeypatch.setattr(run_galaxy.log, "normal", said.append)
    monkeypatch.setattr(store, "get_connection", lambda config=None, **options: a_store_connection(DeadConnection()))
    monkeypatch.setattr(run_galaxy, "_mysql_config_of", lambda args: None)

    def chart(conn, center, limit, on_batch=None, progress=None):
        calls.append(1)
        if len(calls) == 1:
            raise pymysql.err.OperationalError(2013, "Lost connection to MySQL server during query")
        return {"tiers": [], "stars": 1, "seconds": 1.0, "cell": (0, 0, 0), "floors": []}

    monkeypatch.setattr(sky_charting, "chart", chart)
    run_galaxy._sky_unit({"database": {}, "sky": {"center": [0, 0, 0], "limit": 7.5}}, lambda done, total: None)
    assert any("MySQL connection lost: MySQL was restarted" in line for line in said)


def test_facts_reads_a_real_server(mysql_config):
    from planetgen.db import connection_loss as loss
    conn = store.get_connection(mysql_config)
    try:
        found = loss.facts(conn)
    finally:
        conn.close()
    assert found["uptime"] > 0 and found["wait_timeout"] > 0 and found["max_allowed_packet"] > 0 and found["aborted_clients"] >= 0
