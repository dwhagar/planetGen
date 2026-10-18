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
from tests.test_api import _place_sector, _scatter_stars


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


def test_a_scatter_or_a_fill_of_scattered_phenomena_moves_the_stored_phenomena_count(counting, mysql_config):
    """UX.110: the content stamp follows sectors and systems only, so a stored phenomena count kept its old figure
    after a scatter, a fill that built scattered phenomena or a saved phenomenon -- the overview showed it for hours."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.phenomena.nebula import Nebula

    def phenomena():
        conn = query.open_readonly(mysql_config)
        try:
            return query.count_phenomena(conn)
        finally:
            conn.close()

    nebula = Nebula(SystemConfig())
    nebula.name = "Alpha Cloud"
    store.save_phenomenon(nebula, SystemConfig(), "nebula", config=mysql_config)
    assert phenomena() == 1
    conn = store.get_connection(mysql_config)
    try:
        store.insert_phenomenon_scatter(conn, [(0, 0, slot, "neutron-star", None, slot, 0, 0, None, None, None,
                                                slot + 1, None, None) for slot in range(5)])
        store.record_phenomenon_scatter_classes(conn, {("neutron-star", ""): 5})
        conn.commit()
        assert _until(lambda: phenomena() == 6)
        store.mark_phenomena_built(conn, [conn.execute("SELECT MIN(id) AS i FROM phenomenon_scatter").fetchone()["i"]])
        conn.commit()
    finally:
        conn.close()
    assert _until(lambda: phenomena() == 5)
    other = Nebula(SystemConfig())
    other.name = "Zeta Cloud"
    store.save_phenomenon(other, SystemConfig(), "nebula", config=mysql_config)
    assert _until(lambda: phenomena() == 6)


# --- The uncharted-star count (PERF.112) ------------------------------------------------------------------------


def _uncharted(mysql_config):
    conn = query.open_readonly(mysql_config)
    try:
        return query.count_uncharted_systems(conn)
    finally:
        conn.close()


def _write(mysql_config, work):
    conn = store.get_connection(mysql_config)
    try:
        work(conn)
        conn.commit()
    finally:
        conn.close()


def test_the_uncharted_count_is_exact_first_and_read_from_the_stored_answer_after(counting, mysql_config):
    _scatter_stars(mysql_config)
    assert _uncharted(mysql_config) == 3
    conn = query.open_readonly(mysql_config)
    seen = []
    real = conn.execute
    conn.execute = lambda sql, params=(): seen.append(" ".join(sql.split())) or real(sql, params)
    try:
        assert query.count_uncharted_systems(conn) == 3 and query.count_uncharted_systems(conn) == 3
    finally:
        conn.close()
    assert not [sql for sql in seen if "star_system_id IS NULL" in sql], "the scatter table was scanned again"


def test_a_request_without_a_stored_answer_scans_at_most_the_cap(counting, mysql_config, monkeypatch):
    _scatter_stars(mysql_config)
    monkeypatch.setattr(countcache, "WAIT_SECONDS", 0.0)
    monkeypatch.setattr(countcache, "_start", lambda *a: None)
    conn = query.open_readonly(mysql_config)
    seen = []
    real = conn.execute
    conn.execute = lambda sql, params=(): seen.append(" ".join(sql.split())) or real(sql, params)
    try:
        assert query.count_uncharted_systems(conn) == 3
    finally:
        conn.close()
    scans = [sql for sql in seen if "star_system_id IS NULL" in sql]
    assert scans and all("LIMIT" in sql for sql in scans)


def test_the_uncharted_count_drops_after_a_star_is_built_into_a_system(counting, mysql_config):
    ids = _scatter_stars(mysql_config)
    assert _uncharted(mysql_config) == 3
    system_id = _place_sector(mysql_config, "Built", (5.0, 5.0, 5.0))
    _write(mysql_config, lambda conn: store.mark_bright_star_filled(
        conn, ids[0], conn.execute("SELECT MAX(id) AS n FROM star_systems").fetchone()["n"]))
    assert _until(lambda: _uncharted(mysql_config) == 2)
    assert system_id


def test_the_uncharted_count_follows_a_rescatter_and_a_reset(counting, mysql_config):
    _scatter_stars(mysql_config)
    assert _uncharted(mysql_config) == 3
    _write(mysql_config, store.clear_bright_stars)
    assert _until(lambda: _uncharted(mysql_config) == 0)
    _scatter_stars(mysql_config, (800, 600, 900, 700, 500))
    _write(mysql_config, lambda conn: store.record_bright_star_scatter(conn, 1.0e4, 7))
    assert _until(lambda: _uncharted(mysql_config) == 5)
    _write(mysql_config, store.clear_bright_stars)
    assert _until(lambda: _uncharted(mysql_config) == 0)
