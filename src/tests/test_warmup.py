# tests/test_warmup.py

"""
Tests for `planetgen/web/warmup.py` (MAP.134): the Galaxy Map's opening
view built into the tile cache ahead of the first visit. The API calls
are fakes, as in test_tilecache.py.
"""

from planetgen.web import warmup
from planetgen.web.lib import tilecache
from tests.test_tilecache import FakeApi, api  # noqa: F401 -- the fixture


def test_warming_caches_the_exact_tiles_and_stage_the_first_visit_asks_for(api, monkeypatch):  # noqa: F811
    monkeypatch.setattr(warmup.apiclient, "get_galaxy_shape", lambda db: None)
    keys = warmup.opening_request(None)

    done = warmup.warm_opening_view("mydb")
    assert done["tiles"] == len(keys) and done["cached"] == 0
    assert api.tile_calls == [keys] and api.stage_calls == [None]

    # The first visit now reads everything from disk.
    visit = tilecache.fetch_tiles("mydb", keys)
    tilecache.fetch_stage("mydb")
    assert visit["cached"] == len(keys)
    assert api.tile_calls == [keys] and api.stage_calls == [None]

    again = warmup.warm_opening_view("mydb")
    assert again["cached"] == len(keys)


def test_a_new_generation_is_warmed_again(api, monkeypatch):  # noqa: F811
    monkeypatch.setattr(warmup.apiclient, "get_galaxy_shape", lambda db: None)
    warmup.warm_opening_view("mydb")
    api.stamp = "00000000000000bb"
    api.changed = None  # a release or a re-plan: everything changes
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    assert warmup.warm_opening_view("mydb")["cached"] == 0
