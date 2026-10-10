# tests/test_count_cache.py

"""
Stored counts for the table pages (PERF.64, `planetgen/db/countcache.py`): a
page no longer counts a whole table inside the request, and a statement that
hits the web time limit is named in the activity log.
"""

import time

import pytest

from planetgen.admin import activity_log
from planetgen.db import countcache, query, store
from tests.test_api import _place_sector


@pytest.fixture
def counting(monkeypatch, redis_server):
    """The stored counts on, with nothing left over from another test."""
    monkeypatch.setenv(countcache.ENV_VAR, "on")
    monkeypatch.setattr(countcache, "MIN_INTERVAL_SECONDS", 0.0)
    countcache.forget_all()
    yield
    countcache.forget_all()


def _until(check, timeout=30):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        value = check()
        if value:
            return value
        time.sleep(0.1)
    raise AssertionError("never became true")


def _systems(mysql_config):
    conn = query.open_readonly(mysql_config)
    try:
        return query.count_systems(conn)
    finally:
        conn.close()


def test_the_first_count_is_exact_and_a_later_change_shows_after_the_next_background_count(counting, mysql_config):
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    assert _systems(mysql_config) == 1
    _place_sector(mysql_config, "Two", (50.0, 5.0, 5.0))
    # The old answer is served at once; the new one is being made.
    first = _systems(mysql_config)
    assert first in (1, 2)
    assert _until(lambda: _systems(mysql_config) == 2)


def test_a_count_under_an_unchanged_stamp_is_not_made_again(counting, mysql_config, monkeypatch):
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    assert _systems(mysql_config) == 1
    started = []
    real = countcache._start
    monkeypatch.setattr(countcache, "_start", lambda *a: started.append(1) or real(*a))
    assert _systems(mysql_config) == 1 and _systems(mysql_config) == 1
    assert started == []


def test_a_request_never_counts_the_table_itself_when_the_answer_is_not_there(counting, mysql_config, monkeypatch):
    """Without a stored answer the request returns an estimate rather than running a table-wide count."""
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    monkeypatch.setattr(countcache, "WAIT_SECONDS", 0.0)
    monkeypatch.setattr(countcache, "_start", lambda *a: None)  # no background count either
    conn = query.open_readonly(mysql_config)
    seen = []
    real = conn.execute
    conn.execute = lambda sql, params=(): seen.append(" ".join(sql.split())) or real(sql, params)
    try:
        assert query.count_systems(conn) >= 0
        assert query.count_sectors(conn) >= 0
        assert query.systems_facets(conn) == {"placement": [], "binary": [], "octant": []}
        assert query.sectors_facets(conn) == {"quadrant": []}
    finally:
        conn.close()
    assert not [sql for sql in seen if "COUNT(" in sql.upper() and "LIMIT" not in sql.upper()
                and "information_schema" not in sql and "SELECT COALESCE(MAX" not in sql], seen


def test_the_stored_facets_are_the_exact_ones(counting, mysql_config):
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    conn = query.open_readonly(mysql_config)
    try:
        facets = query.systems_facets(conn)
        assert [option["value"] for option in facets["placement"]] == ["sector"]
        assert facets["placement"][0]["count"] == 1
    finally:
        conn.close()


def test_a_slow_first_count_gives_the_estimate_and_the_count_arrives_later(counting, mysql_config, monkeypatch):
    monkeypatch.setattr(countcache, "WAIT_SECONDS", 0.05)
    conn = query.open_readonly(mysql_config)
    try:
        def slow(c):
            time.sleep(1.0)
            return 7

        assert countcache.cached(conn, ["slow"], slow, lambda c: -1, lambda c: "s") == -1
        assert _until(lambda: countcache.cached(conn, ["slow"], slow, lambda c: -1, lambda c: "s") == 7)
    finally:
        conn.close()


def test_the_cache_can_be_turned_off(mysql_config, monkeypatch):
    monkeypatch.setenv(countcache.ENV_VAR, "off")
    conn = query.open_readonly(mysql_config)
    try:
        assert countcache.cached(conn, ["x"], lambda c: 3, lambda c: -1, lambda c: "s") == 3
    finally:
        conn.close()


def test_a_statement_that_hits_the_time_limit_is_named_in_the_activity_log(mysql_config, monkeypatch):
    import pymysql
    events = []
    monkeypatch.setattr(activity_log, "event", lambda *args, **fields: events.append((args, fields)))
    conn = store.get_connection(mysql_config, ensure_schema=False, statement_timeout_s=0.2)
    try:
        with pytest.raises(pymysql.err.OperationalError):
            conn.execute("SELECT COUNT(*) FROM (SELECT SLEEP(5)) t")
    finally:
        conn.close()
    assert events and events[0][0][:2] == ("DB", "statement_timeout")
    assert "SLEEP(5)" in events[0][1]["sql"] and float(events[0][1]["seconds"]) < 4
