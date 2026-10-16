# tests/test_cache_invalidation.py

"""
PERF.104: a stale copy is never served as current, whatever changed it, and a change rebuilds
only the caches that read it (docs/design/cache-warm-up.md).

The pages of a real app on a throwaway database are read twice around each kind of change:
a write through the API (a sector rename, a system rename, a deletion), a change the command
line records (`dataversion.record`: generation, the orbital update, the population pass), a
direct database edit nothing recorded (found by the galaxy's content stamp), and a Redis that
lost its data. A spy on `apiclient._send` says which answers were rebuilt.
"""

import json
import secrets

import pytest

from planetgen import web
from planetgen.db import dataversion, store
from planetgen.galaxy.sector import SpaceSector
from planetgen.queue import redisclient
from planetgen.web.lib import apiclient, pagecache, pagetags
from tests.publicids import pid
from tests.test_api import admin_client, client, default_admin_client, first_admin_password, seeded_sector  # noqa: F401


# --- The two tables ----------------------------------------------------------------------

def test_a_sector_page_reads_its_neighbours_and_nothing_else_by_id():
    body = json.dumps({"id": "A", "neighbors": [{"sector_id": "b-1", "exists": True}, {"sector_id": None}]})
    tags = pagetags.read_tags("/sectors/A", body)
    assert "sector:A" in tags and "sector:B-1" in tags
    assert not any(tag.startswith("system:") for tag in tags)
    assert tags[0] == dataversion.TAG_EPOCH


def test_a_system_page_reads_its_system_and_its_sector():
    tags = pagetags.read_tags("/systems/S-1", json.dumps({"id": "S-1", "sector_id": "AA"}))
    assert "system:S-1" in tags and "sector:AA" in tags
    assert dataversion.TAG_SECTOR_SET not in tags
    # A sector that is not named means every set could matter.
    assert dataversion.TAG_SECTOR_SET in pagetags.read_tags("/systems/S-1", json.dumps({"id": "S-1"}))


def test_an_unknown_path_reads_everything():
    tags = pagetags.read_tags("/something-new", "{}")
    assert set(pagetags.ALL_TAGS) <= set(tags) and dataversion.TAG_EPOCH in tags


@pytest.mark.parametrize("method, path, kind", [
    ("PATCH", "/sectors/AA", "sector-edited"),
    ("DELETE", "/sectors/AA", "sector-set"),
    ("POST", "/sectors", "sector-set"),
    ("DELETE", "/sectors/AA/contents", "sector-contents"),
    ("POST", "/sectors/AA/regenerate", "sector-contents"),
    ("POST", "/sectors/AA/generate-neighborhood", "generated"),
    ("PATCH", "/systems/S-1", "system-edited"),
    ("DELETE", "/systems/S-1", "system-set"),
    ("POST", "/systems", "system-set"),
    ("PATCH", "/planets/P-1", "system-any"),
    ("POST", "/facilities", "system-any"),
    ("DELETE", "/phenomena/nebulae/N-1", "phenomena-set"),
    ("POST", "/phenomena/nebulae/N-1/regenerate", "phenomenon-edited"),
    ("POST", "/admin/naming-key", "everything"),
])
def test_each_write_maps_to_its_change(method, path, kind):
    assert pagetags.change_for_write(method, path)[0] == kind


def test_every_change_kind_has_tags():
    for kind in dataversion.CHANGES:
        ids = {"sector_id": "AA", "system_id": "S-1", "phenomenon": ("nebulae", "N-1")}
        assert dataversion.tags_for_change(kind, **ids)
    with pytest.raises(ValueError):
        dataversion.tags_for_change("nonsense")
    with pytest.raises(ValueError):
        dataversion.tags_for_change("sector-edited")


# --- Redis ---------------------------------------------------------------------------------

@pytest.fixture
def database(redis_server):
    name = f"pgtest_{secrets.token_hex(6)}"
    yield name
    for key in redisclient.client(decode_responses=True).scan_iter(f"planetgen:*:{name}*"):
        redisclient.client().delete(key)


def test_bumps_are_seen_by_every_reader(database):
    first, found = dataversion.versions(database, ["sector:AA"])
    assert found == {"sector:AA": 0}
    dataversion.record(database, "sector-edited", sector_id="aa")
    second, found = dataversion.versions(database, ["sector:AA", dataversion.TAG_SECTOR_DATA, "sector:BB"])
    assert second == first and found == {"sector:AA": 1, dataversion.TAG_SECTOR_DATA: 1, "sector:BB": 0}


def test_a_redis_that_lost_its_data_starts_a_new_incarnation(database):
    first, _ = dataversion.versions(database, [])
    redisclient.client().delete(dataversion.key(database))
    second, _ = dataversion.versions(database, [])
    assert first != second


def test_the_first_process_to_claim_a_state_bumps_and_the_others_skip(database):
    assert dataversion.claim_state(database, "s1", None)
    assert not dataversion.claim_state(database, "s1", None)
    assert dataversion.claim_state(database, "s2", "s1")
    assert not dataversion.claim_state(database, "s3", "s1")
    assert dataversion.recorded_state(database) == "s2"


def test_a_copy_one_process_made_is_served_by_another_until_data_changes(database):
    state = {"value": "a"}

    def changes_for(db, since):
        return {"stamp": state["value"], "state": f"state-{state['value']}"}

    def make():
        return pagecache.ResponseCache(changes_for, dict(pagecache.DEFAULTS, stamp_seconds=0))

    warm, web_worker = make(), make()
    assert warm.get(database, "/sectors/AA?db=x") is None
    warm.put(database, "/sectors/AA?db=x", json.dumps({"id": "AA", "neighbors": []}), warm.begin(database))
    assert len(web_worker) == 0
    assert json.loads(web_worker.get(database, "/sectors/AA?db=x"))["id"] == "AA"  # from Redis, not rebuilt
    dataversion.record(database, "sector-edited", sector_id="AA")
    assert web_worker.get(database, "/sectors/AA?db=x") is None  # the copy it holds is stale now
    assert warm.get(database, "/sectors/AA?db=x") is None
    dataversion.record(database, "sector-edited", sector_id="ZZ")


# --- The pages -----------------------------------------------------------------------------

class _Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


@pytest.fixture
def site(admin_client, seeded_sector, monkeypatch):  # noqa: F811
    """The admin's app on the seeded database: its page cache on a fake clock, and a spy on what is rebuilt."""
    config, sector_id, system_ids = seeded_sector
    other_id = store.save_sector(SpaceSector(name="Other Sector", edge_ly=10.0), config=config)
    app = admin_client.application
    cache = app.extensions[web.PAGE_CACHE_EXTENSION]
    cache._clock = _Clock()
    seen = []
    real = apiclient._send

    def counting(method, target, *args, **kwargs):
        if method == "GET":
            seen.append(target.partition("?")[0])
        return real(method, target, *args, **kwargs)
    monkeypatch.setattr(apiclient, "_send", counting)

    class Site:
        pass

    found = Site()
    found.client, found.config, found.cache, found.seen = admin_client, config, cache, seen
    found.sector, found.other, found.systems = sector_id, other_id, system_ids
    found.sector_pid, found.other_pid = pid("sector", sector_id, config), pid("sector", other_id, config)
    found.system_pids = [pid("system", system_id, config) for system_id in system_ids]
    yield found
    for key in redisclient.client(decode_responses=True).scan_iter(f"planetgen:*:{config.database}*"):
        redisclient.client().delete(key)


def _visit(site):
    """Reads both sector pages and both system pages, and says which API answers were rebuilt."""
    del site.seen[:]
    pages = {}
    for name, path in (("sector", f"/sector/{site.sector_pid}"), ("other", f"/sector/{site.other_pid}"),
                       ("system0", f"/system/{site.system_pids[0]}"), ("system1", f"/system/{site.system_pids[1]}")):
        response = site.client.get(path)
        assert response.status_code == 200, (path, response.status_code)
        pages[name] = response.get_data(as_text=True)
    return pages, list(site.seen)


def _rebuilt(site, seen, *, sector=False, other=False, system0=False, system1=False):
    wanted = {f"/sectors/{site.sector_pid}": sector, f"/sectors/{site.other_pid}": other,
              f"/systems/{site.system_pids[0]}": system0, f"/systems/{site.system_pids[1]}": system1}
    for path, rebuilt in wanted.items():
        assert (path in seen) == rebuilt, (path, "rebuilt" if path in seen else "from the cache", seen)


def test_a_second_visit_rebuilds_nothing(site):
    _visit(site)
    _pages, seen = _visit(site)
    _rebuilt(site, seen)


def test_renaming_a_sector_through_the_api_rebuilds_only_what_read_it(site):
    _visit(site)
    response = site.client.patch(f"/api/sectors/{site.sector_pid}", json={"name": "Renamed By The Admin"})
    assert response.status_code == 200, response.get_data(as_text=True)
    pages, seen = _visit(site)
    assert "Renamed By The Admin" in pages["sector"]
    # The sector, and the systems in it (their pages name their sector); not the other sector.
    _rebuilt(site, seen, sector=True, system0=True, system1=True)


def test_renaming_a_system_rebuilds_its_page_and_its_sectors_not_another_sector(site):
    _visit(site)
    response = site.client.patch(f"/api/systems/{site.system_pids[0]}", json={"name": "Renamed System"})
    assert response.status_code == 200, response.get_data(as_text=True)
    pages, seen = _visit(site)
    assert "Renamed System" in pages["system0"]
    assert "Renamed System" in pages["sector"]
    # Its sibling names it in "other systems here", so it is rebuilt too; the other sector is not.
    _rebuilt(site, seen, sector=True, system0=True, system1=True)


def test_deleting_a_system_is_never_served_again(site):
    _visit(site)
    assert site.client.delete(f"/api/systems/{site.system_pids[0]}").status_code == 200
    assert site.client.get(f"/system/{site.system_pids[0]}").status_code == 404
    pages = site.client.get(f"/sector/{site.sector_pid}").get_data(as_text=True)
    assert site.system_pids[0] not in pages


def test_deleting_a_sector_is_never_served_again(site):
    _visit(site)
    assert site.client.delete(f"/api/sectors/{site.other_pid}").status_code == 200
    assert site.client.get(f"/sector/{site.other_pid}").status_code == 404


def test_the_population_pass_rebuilds_what_shows_population_not_the_sector_lists(site):
    _visit(site)
    dataversion.record(site.config.database, "population")
    _pages, seen = _visit(site)
    _rebuilt(site, seen, sector=True, other=True, system0=True, system1=True)


def test_the_orbital_update_rebuilds_what_shows_positions(site):
    _visit(site)
    dataversion.record(site.config.database, "orbits")
    _pages, seen = _visit(site)
    _rebuilt(site, seen, sector=True, other=True, system0=True, system1=True)


def test_a_change_to_another_sector_leaves_these_alone(site):
    _visit(site)
    dataversion.record(site.config.database, "sector-edited", sector_id=site.other_pid)
    _pages, seen = _visit(site)
    _rebuilt(site, seen, other=True)


def test_a_database_edit_nothing_recorded_is_found_by_the_stamp(site):
    pages, _seen = _visit(site)
    conn = store.get_connection(site.config)
    try:
        conn.execute("UPDATE sectors SET name = ? WHERE id = ?", ("Edited Behind The Cache", site.sector))
        conn.commit()
    finally:
        conn.close()
    pages, seen = _visit(site)
    assert "Edited Behind The Cache" not in pages["sector"]  # the stamp is not due yet
    site.cache._clock.now += site.cache._stamp_seconds + 1
    pages, seen = _visit(site)
    assert "Edited Behind The Cache" in pages["sector"]
    _rebuilt(site, seen, sector=True, system0=True, system1=True)


def test_a_redis_that_lost_its_data_serves_nothing_it_held(site):
    _visit(site)
    redisclient.client().delete(dataversion.key(site.config.database))
    for key in redisclient.client(decode_responses=True).scan_iter(f"planetgen:page:{site.config.database}:*"):
        redisclient.client().delete(key)
    _pages, seen = _visit(site)
    _rebuilt(site, seen, sector=True, other=True, system0=True, system1=True)


def test_a_redis_that_refuses_a_bump_makes_every_copy_stale_instead(database, monkeypatch, capsys):
    first, _ = dataversion.versions(database, ["sector:AA"])

    class _Refusing:
        def pipeline(self, transaction=True):
            raise dataversion.redis.ResponseError("OOM command not allowed when used memory > 'maxmemory'")

        def __getattr__(self, name):
            return getattr(redisclient.client(timeout=0.5, retries=1, decode_responses=True), name)

    monkeypatch.setattr(dataversion, "_client", lambda: _Refusing())
    assert dataversion.bump(database, "sector:AA") is True
    monkeypatch.undo()
    second, found = dataversion.versions(database, ["sector:AA"])
    assert second != first and found == {"sector:AA": 0}, "a new incarnation: every copy taken before is stale"
    assert "made stale" in capsys.readouterr().err
