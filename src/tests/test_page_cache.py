# tests/test_page_cache.py

"""
The pages' in-memory API response cache (`planetgen/web/lib/pagecache.py`, PERF.2):
the cache on its own with a fake clock, then the Flask pages serving a
repeat visit without the API, and an API write clearing it.
"""

import json
import secrets

import pytest

from planetgen.web.app import create_app
from planetgen.api.config import Config

from planetgen import web
from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.web.lib import pagecache  # noqa: E402
from planetgen.db import dataversion  # noqa: E402
from planetgen.db import store  # noqa: E402
from planetgen.queue import redisclient  # noqa: E402
from planetgen.galaxy.sector import SpaceSector  # noqa: E402
from tests.publicids import pid


class _Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


@pytest.fixture
def db(redis_server):
    """A database name of this test's own, so its versions and entries in the shared Redis are its own."""
    name = f"pgtest_{secrets.token_hex(6)}"
    yield name
    client = redisclient.client(decode_responses=True)
    for key in client.scan_iter(f"planetgen:*:{name}*"):
        client.delete(key)


def _cache(stamps, clock, **settings):
    """A cache whose changes come from `stamps[db]`: a state string, an exception, or a dict (`stamp`, `state`, ...)."""
    calls = []

    def changes_for(db, since):
        calls.append(db)
        value = stamps[db]
        if isinstance(value, Exception):
            raise value
        if isinstance(value, dict):
            return dict(value, seen=since)
        return {"stamp": value, "state": f"state-{value}"}
    cache = pagecache.ResponseCache(changes_for, dict(pagecache.DEFAULTS, **settings), clock=clock)
    return cache, calls


def _miss_then_put(cache, db, target, body):
    assert cache.get(db, target) is None
    cache.put(db, target, body, cache.begin(db))


def test_hit_until_the_state_changes(db):
    clock = _Clock()
    stamps = {db: "a"}
    cache, calls = _cache(stamps, clock, stamp_seconds=15)
    _miss_then_put(cache, db, "/sectors/1?db=g", '{"id": 1}')
    assert cache.get(db, "/sectors/1?db=g") == '{"id": 1}'
    assert calls == [db]  # one check covers the next 15 s
    stamps[db] = "b"
    clock.now += 10
    assert cache.get(db, "/sectors/1?db=g") == '{"id": 1}'
    clock.now += 10
    assert cache.get(db, "/sectors/1?db=g") is None  # a change nothing places bumps `epoch`
    assert calls == [db, db]


def test_a_change_to_one_sector_rebuilds_only_what_read_it(db):
    clock = _Clock()
    stamps = {db: "a"}
    cache, _ = _cache(stamps, clock, stamp_seconds=0)
    for sector in ("AAA", "BBB"):
        _miss_then_put(cache, db, f"/sectors/{sector}", json.dumps({"id": sector, "neighbors": []}))
    _miss_then_put(cache, db, "/systems/SYS-1", json.dumps({"id": "SYS-1", "sector_id": "BBB"}))
    _miss_then_put(cache, db, "/polities", "[]")
    stamps[db] = {"stamp": "b", "state": "state-b", "full": False, "sectors": ["AAA"], "grew": False}
    assert cache.get(db, "/sectors/AAA") is None
    assert cache.get(db, "/sectors/BBB") is not None
    assert cache.get(db, "/systems/SYS-1") is not None
    assert cache.get(db, "/polities") == "[]"


def test_a_busy_fill_keeps_the_entries_and_the_state_until_they_age_out(db):
    """PERF.38: while the galaxy changes faster than pages can be rebuilt, a new state
    no longer empties the cache at every check."""
    clock = _Clock()
    stamps = {db: "a"}
    cache, calls = _cache(stamps, clock, stamp_seconds=15, max_age_seconds=300)
    _miss_then_put(cache, db, "/t", "body")
    clock.now += 20
    stamps[db] = {"stamp": "b", "state": "state-b", "full": True, "busy": True}
    assert cache.get(db, "/t") is None  # the first busy check is an ordinary change
    cache.put(db, "/t", "body", cache.begin(db))
    clock.now += 20
    stamps[db] = {"stamp": "c", "state": "state-c", "full": True, "busy": True}
    assert cache.get(db, "/t") == "body"  # still busy: a fill, so kept
    clock.now += 20
    stamps[db] = {"stamp": "d", "state": "state-d", "full": True, "busy": True}
    assert cache.get(db, "/t") == "body"
    assert cache._stamps[db][2] == "state-b"  # still comparing against what it accepted
    stamps[db] = {"stamp": "d", "state": "state-d"}
    clock.now += 20
    assert cache.get(db, "/t") is None  # settled: a plain change drops it
    assert cache._stamps[db][2] == "state-d"


def test_a_fill_that_stays_busy_past_max_age_drops_the_entries(db):
    clock = _Clock()
    stamps = {db: "a"}
    cache, _ = _cache(stamps, clock, stamp_seconds=15, max_age_seconds=300)
    _miss_then_put(cache, db, "/t", "body")
    stamps[db] = {"stamp": "b", "state": "state-b", "full": True, "busy": True}
    clock.now += 20
    assert cache.get(db, "/t") is None  # the first busy check
    cache.put(db, "/t", "body", cache.begin(db))
    stamps[db] = {"stamp": "c", "state": "state-c", "full": True, "busy": True}
    clock.now += 100
    assert cache.get(db, "/t") == "body"
    clock.now += 250  # busy since 370 s ago: the state is taken at last
    assert cache.get(db, "/t") is None
    assert cache._stamps[db][2] == "state-c"


def test_a_change_in_one_database_leaves_the_other(db):
    clock = _Clock()
    other = db + "_other"
    stamps = {db: "a", other: "x"}
    cache, _ = _cache(stamps, clock, stamp_seconds=0)
    _miss_then_put(cache, db, "/t1", "1")
    _miss_then_put(cache, other, "/t2", "2")
    stamps[db] = "b"
    assert cache.get(db, "/t1") is None
    assert cache.get(other, "/t2") == "2"


def test_entries_expire_after_max_age(db):
    clock = _Clock()
    cache, _ = _cache({db: "a"}, clock, max_age_seconds=300, stamp_seconds=1e9)
    _miss_then_put(cache, db, "/t", "body")
    clock.now += 299
    assert cache.get(db, "/t") == "body"
    clock.now += 2
    assert cache._entries.get((db, "/t")) is None


def test_unreadable_state_means_uncached(db):
    clock = _Clock()
    cache, _ = _cache({db: RuntimeError("no galaxy tables")}, clock)
    assert cache.get(db, "/t") is None
    cache.put(db, "/t", "body", cache.begin(db))
    assert len(cache) == 0


def test_a_change_beats_a_fetch_in_flight(db):
    clock = _Clock()
    cache, _ = _cache({db: "a"}, clock)
    assert cache.get(db, "/t") is None
    token = cache.begin(db)
    dataversion.record(db, "everything")  # an edit lands while the old answer is being fetched
    cache.put(db, "/t", "old", token)
    assert cache.get(db, "/t") is None


def test_clear_beats_a_fetch_in_flight(db):
    clock = _Clock()
    cache, _ = _cache({db: "a"}, clock)
    assert cache.get(db, "/t") is None
    token = cache.begin(db)
    cache.clear()
    cache.put(db, "/t", "old", token)
    assert cache.get(db, "/t") is None


def test_least_recently_used_goes_first_and_size_is_bounded(db):
    clock = _Clock()
    cache, _ = _cache({db: "a"}, clock, max_entries=2, stamp_seconds=1e9)
    for target in ("/a", "/b"):
        _miss_then_put(cache, db, target, target)
    assert cache.get(db, "/a") == "/a"
    _miss_then_put(cache, db, "/c", "/c")
    assert (db, "/b") not in cache._entries and cache.get(db, "/a") == "/a"
    small, _ = _cache({db: "a"}, clock, max_mb=0.001)
    _miss_then_put(small, db, "/big", "x" * 1000)  # over a quarter of the budget
    assert len(small) == 0


def test_settings():
    from planetgen.util import settings as settings_model
    assert settings_model.build({}, environ={}).page_cache.model_dump() == pagecache.DEFAULTS
    assert settings_model.build({"page_cache": {"max_mb": 8}}, environ={}).page_cache.max_mb == 8
    assert settings_model.build({}, environ={"PLANETGEN_PAGE_CACHE": "off"}).page_cache.enabled is False
    assert pagecache.is_cacheable("/sectors/5")
    assert not pagecache.is_cacheable("/galaxy/changes")
    assert not pagecache.is_cacheable("/galaxy/tiles")


# --- The pages --------------------------------------------------------------------------

@pytest.fixture
def db_app(mysql_config):
    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        WEB_DATABASE = ""
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    application = create_app(RealConfig)
    application.testing = True
    return application


def _count_sends(monkeypatch):
    seen = []
    real = apiclient._send

    def counting(method, target, *args, **kwargs):
        seen.append((method, target.partition("?")[0]))
        return real(method, target, *args, **kwargs)
    monkeypatch.setattr(apiclient, "_send", counting)
    return seen


def test_repeat_visit_is_served_from_the_cache(db_app, mysql_config, monkeypatch, redis_server):
    sector = SpaceSector(name="Cached Sector", edge_ly=40.0)
    sector_id = store.save_sector(sector, config=mysql_config)
    client = db_app.test_client()
    seen = _count_sends(monkeypatch)
    first = client.get(f"/sector/{pid('sector', sector_id)}").get_data(as_text=True)
    assert ("GET", f"/sectors/{pid('sector', sector_id)}") in seen
    del seen[:]
    second = client.get(f"/sector/{pid('sector', sector_id)}").get_data(as_text=True)
    assert "Cached Sector" in first and "Cached Sector" in second
    assert ("GET", f"/sectors/{pid('sector', sector_id)}") not in seen


def test_a_write_through_the_api_makes_what_it_changed_stale(db_app, mysql_config, redis_server):
    @db_app.route("/api/test-write", methods=["POST"])
    def test_write():
        return {"ok": True}

    @db_app.route("/api/test-refused", methods=["POST"])
    def test_refused():
        return {"error": "no"}, 409

    sector = SpaceSector(name="Before", edge_ly=40.0)
    sector_id = store.save_sector(sector, config=mysql_config)
    client = db_app.test_client()
    client.get(f"/sector/{pid('sector', sector_id)}")
    cache = db_app.extensions[web.PAGE_CACHE_EXTENSION]
    database = mysql_config.database
    target = f"/sectors/{pid('sector', sector_id)}?db={database}"
    assert cache.get(database, target) is not None
    # A refused write changes nothing, so the copy stands...
    assert client.post("/api/test-refused").status_code == 409
    assert cache.get(database, target) is not None
    # ...and a successful one nothing describes makes every copy stale.
    assert client.post("/api/test-write").status_code == 200
    assert cache.get(database, target) is None


def test_cache_can_be_turned_off(mysql_config, monkeypatch):
    monkeypatch.setenv("PLANETGEN_PAGE_CACHE", "off")

    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        SECRET_KEY = "test-secret"
    assert web.PAGE_CACHE_EXTENSION not in create_app(RealConfig).extensions
