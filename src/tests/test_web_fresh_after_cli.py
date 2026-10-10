# tests/test_web_fresh_after_cli.py

"""
TEST.42: the pages show what `planetgen` wrote straight to the
database, not through the API.

An API write clears the page cache at once (`test_page_cache.py`), but a
command-line run happens in another process, so the only things that
can notice it are the galaxy's content stamp (`GET /api/galaxy/changes`)
-- checked at most once per `page_cache.stamp_seconds` by the page cache
(`planetgen/web/lib/pagecache.py`) and once per `STAMP_TTL_SECONDS` by the Galaxy Map's
disk tile cache (`planetgen/web/lib/tilecache.py`) -- and the page cache's age limit.

Each test warms the pages and both caches, runs `planetgen` through
its real entry point (`main()`, via `sys.argv`, as
`test_galaxy_gen.py` does) against the same throwaway database, lets the
stamp checks fall due (a fake clock for the page cache, `stamp.json`'s
mtime moved back for the tile cache), and fetches again:

- `/galaxy`, `/galaxy/tiles` and `/api/galaxy/tiles` show sectors the
  `galaxy` subcommand placed;
- `/sector/<id>` links the neighbors generated next to it;
- `/sectors` and `/systems` list a new unplaced sector, and `/systems`
  a new standalone system;
- `/sector/<id>` and the tiles show a nebula `phenomenon --sector-id`
  added to it.

Before the stamp check falls due the old answer is still served -- that
is what proves the pages really were cached and not just read afresh.
Skipped without a MySQL server.
"""

import math
import os
import re
import sys

import pytest
from markupsafe import escape

from planetgen.web.app import create_app
from planetgen.api.config import Config

from planetgen.cli import generate as generate_cli
from planetgen import web
from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.web.lib import pagecache  # noqa: E402
from planetgen.web.lib import tilecache  # noqa: E402
from planetgen.db import store  # noqa: E402
from planetgen import tuning
from planetgen.galaxy.density import build_galaxy_shape  # noqa: E402
from planetgen.physics.units import ly_to_pc  # noqa: E402

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)

# The same small toy galaxy test_galaxy_gen.py plans: ring 0 holds three
# slots, all deep inside it, so `galaxy --ring 0` is quick.
_SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)

WHOLE_GALAXY_TILE = "0/0/0/0"


class _Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


class Site:
    """The Flask app on the fixture's database, with its page cache on a
    fake clock and its tile cache in `tile_root`."""

    def __init__(self, app, clock, tile_root, mysql_config):
        self.app = app
        self.client = app.test_client()
        self.clock = clock
        self.tile_root = tile_root
        self.mysql_config = mysql_config

    @property
    def page_cache(self):
        return self.app.extensions[web.PAGE_CACHE_EXTENSION]

    def get(self, path):
        response = self.client.get(path)
        assert response.status_code == 200, (path, response.status_code)
        return response

    def html(self, path):
        return self.get(path).get_data(as_text=True)

    def let_stamp_checks_fall_due(self):
        """Time passes: past the page cache's `stamp_seconds` (but well
        short of `max_age_seconds`, so only the stamp can drop an entry)
        and past the tile cache's `STAMP_TTL_SECONDS`."""
        self.clock.now += self.page_cache._stamp_seconds + 1
        assert self.page_cache._stamp_seconds + 1 < self.page_cache._max_age
        past = os.path.getmtime(self._stamp_file()) - tilecache.STAMP_TTL_SECONDS - 1
        os.utime(self._stamp_file(), (past, past))

    def _stamp_file(self):
        root = tilecache.cache_dir()
        assert root == str(self.tile_root)
        path = os.path.join(tilecache._db_dir(root, self.mysql_config.database), "stamp.json")
        assert os.path.exists(path), "the tile cache was never warmed"
        return path


@pytest.fixture
def site(mysql_config, tmp_path, monkeypatch):
    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        WEB_DATABASE = ""
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")

    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    monkeypatch.delenv("PLANETGEN_PAGE_CACHE", raising=False)
    monkeypatch.delenv("PLANETGEN_TILE_CACHE_MAX_MB", raising=False)
    tile_root = tmp_path / "tiles"
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tile_root))
    app = create_app(RealConfig)
    app.testing = True
    clock = _Clock()
    cache = app.extensions[web.PAGE_CACHE_EXTENSION]
    assert isinstance(cache, pagecache.ResponseCache)
    cache._clock = clock
    return Site(app, clock, tile_root, mysql_config)


def _mysql_argv(mysql_config):
    return [
        "--mysql-host", mysql_config.host,
        "--mysql-port", str(mysql_config.port),
        "--mysql-user", mysql_config.user,
        "--mysql-password", mysql_config.password,
        "--mysql-database", mysql_config.database,
    ]


def _generate(mysql_config, *argv):
    """Runs `planetgen <argv...>` in-process, as from a shell."""
    old_argv = sys.argv
    try:
        sys.argv = ["planetgen", *argv, *_mysql_argv(mysql_config), "--quiet"]
        generate_cli.main()
    finally:
        sys.argv = old_argv


def _plan_galaxy(mysql_config):
    store.save_galaxy_shape(_SHAPE, edge_pc=EDGE_PC, outer_ring_index=999,
                          expected_system_count_at_density_1=1.0, config=mysql_config)
    store.replace_galaxy_layers([(layer, 999) for layer in range(5, -6, -1)], config=mysql_config)


def _rows(mysql_config, sql):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute(sql).fetchall()
    finally:
        conn.close()


def _sectors(mysql_config):
    return _rows(mysql_config, "SELECT id, name, center_x_pc FROM sectors ORDER BY id")


def _placed_ids(tile):
    return sorted(entry["id"] for entry in tile["placed"])


def _placed_count(html):
    match = re.search(r"(\d+) placed sectors?", html)
    assert match, "no placed-sector count on /galaxy"
    return int(match.group(1))


def test_galaxy_map_and_tiles_show_sectors_the_cli_placed(site, mysql_config):
    _plan_galaxy(mysql_config)
    _generate(mysql_config, "galaxy", "--ring", "0", "--limit", "1", "--num-systems", "1", "--backfill-from", "none")
    [first] = _sectors(mysql_config)

    tiles_url = f"/galaxy/tiles?tiles={WHOLE_GALAXY_TILE}"
    api_tiles_url = f"/api/galaxy/tiles?tiles={WHOLE_GALAXY_TILE}"
    # Warm the page cache and the disk tile cache.
    assert _placed_count(site.html("/galaxy")) == 1
    warm = site.get(tiles_url).get_json()
    assert "placed" not in warm["tiles"][WHOLE_GALAXY_TILE], "MAP.157: the web tile carries only what the page reads"
    assert site.get(tiles_url).get_json()["cached"] == 1
    assert _placed_ids(site.get(api_tiles_url).get_json()["tiles"][WHOLE_GALAXY_TILE]) == [first["id"]]

    _generate(mysql_config, "galaxy", "--ring", "0", "--num-systems", "1", "--backfill-from", "none")
    sectors = _sectors(mysql_config)
    assert len(sectors) == 3
    all_ids = sorted(row["id"] for row in sectors)

    # The API itself caches nothing: fresh at once.
    assert _placed_ids(site.get(api_tiles_url).get_json()["tiles"][WHOLE_GALAXY_TILE]) == all_ids
    # The web caches still trust their last stamp check...
    assert _placed_count(site.html("/galaxy")) == 1
    stale = site.get(tiles_url).get_json()
    assert stale["cached"] == 1 and stale["stamp"] == warm["stamp"]

    # ...until it falls due, and then show the new sectors.
    site.let_stamp_checks_fall_due()
    page = site.html("/galaxy")
    assert _placed_count(page) == 3
    fresh = site.get(tiles_url).get_json()
    assert fresh["stamp"] != warm["stamp"]
    assert fresh["cached"] == 0
    # The changed tile was cached again.
    again = site.get(tiles_url).get_json()
    assert again["cached"] == 1 and again["stamp"] == fresh["stamp"]
    # The browser holding the old stamp is told which tiles to drop.
    told = site.get(f"{tiles_url}&stamp={warm['stamp']}").get_json()
    assert WHOLE_GALAXY_TILE in [key for entry in told["history"] for key in entry["tiles"]]
    assert _placed_ids(site.get(api_tiles_url).get_json()["tiles"][WHOLE_GALAXY_TILE]) == all_ids


def test_sector_map_shows_neighbors_the_cli_generated(site, mysql_config):
    """The sector page's map draws its neighbors from `/sector/<id>/scene`
    (MAP.68), cached like a page, so a sector the CLI generated shows
    once the stamp is checked again."""
    _plan_galaxy(mysql_config)
    _generate(mysql_config, "galaxy", "--ring", "0", "--limit", "1", "--num-systems", "1")
    [first] = _sectors(mysql_config)
    scene_url = f"/sector/{first['id']}/scene"
    assert escape(first["name"]) in site.html(f"/sector/{first['id']}")

    def linked():
        return {entry["sectorId"] for entry in site.get(scene_url).get_json()["neighbors"] if entry.get("sectorId")}

    assert linked() == set()

    _generate(mysql_config, "galaxy", "--ring", "0", "--num-systems", "1")
    others = [row for row in _sectors(mysql_config) if row["id"] != first["id"]]
    assert len(others) == 2
    assert linked() == set(), "still the cached answer"
    site.let_stamp_checks_fall_due()
    assert linked() == {row["id"] for row in others}


def test_sector_and_system_lists_show_cli_sectors_and_systems(site, mysql_config):
    _plan_galaxy(mysql_config)
    _generate(mysql_config, "galaxy", "--ring", "0", "--limit", "1", "--num-systems", "1")
    [first] = _sectors(mysql_config)
    old_systems = {row["name"] for row in _rows(mysql_config, "SELECT name FROM star_systems")}
    sectors_before = site.html("/sectors")
    systems_before = site.html("/systems")
    assert escape(first["name"]) in sectors_before
    assert all(escape(name) in systems_before for name in old_systems)
    site.html(f"/galaxy/tiles?tiles={WHOLE_GALAXY_TILE}")

    _generate(mysql_config, "sector", "--name", "Fresh From The Shell", "--num-systems", "2")
    _generate(mysql_config, "system", "--name", "Lonely Shell Star", "-planets")
    new_systems = {row["name"] for row in _rows(mysql_config, "SELECT name FROM star_systems")} - old_systems
    assert "Lonely Shell Star" in new_systems and len(new_systems) == 3

    assert site.html("/sectors") == sectors_before
    assert site.html("/systems") == systems_before
    site.let_stamp_checks_fall_due()
    sectors_after = site.html("/sectors")
    systems_after = site.html("/systems")
    assert "Fresh From The Shell" in sectors_after and escape(first["name"]) in sectors_after
    for name in new_systems | old_systems:
        assert escape(name) in systems_after, name



def test_sector_page_and_tiles_show_a_phenomenon_the_cli_added(site, mysql_config):
    """`planetgen phenomenon --sector-id` adds a row to no table the
    stamp reads, so it must touch its sector (`store.save_phenomenon`):
    before it did, the sector page stayed stale until the page cache's age
    limit, and a cached tile never showed the new cloud at all."""
    _plan_galaxy(mysql_config)
    _generate(mysql_config, "galaxy", "--ring", "0", "--limit", "1", "--num-systems", "1")
    [first] = _sectors(mysql_config)
    page_url = f"/sector/{first['id']}"
    tiles_url = f"/galaxy/tiles?tiles={WHOLE_GALAXY_TILE}"
    before = site.html(page_url)
    warm = site.get(tiles_url).get_json()
    assert "Shell Cloud" not in before
    assert warm["tiles"][WHOLE_GALAXY_TILE]["clouds"] == []

    _generate(mysql_config, "phenomenon", "--type", "nebula", "--sector-id", str(first["id"]),
              "--name", "Shell Cloud")
    site.let_stamp_checks_fall_due()
    assert "Shell Cloud" in site.html(page_url)
    fresh = site.get(tiles_url).get_json()
    assert fresh["stamp"] != warm["stamp"]
    assert [cloud["name"] for cloud in fresh["tiles"][WHOLE_GALAXY_TILE]["clouds"]] == ["Shell Cloud"]
