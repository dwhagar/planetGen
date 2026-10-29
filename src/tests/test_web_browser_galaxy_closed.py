"""
PERF.149: the Galaxy Map page (`/galaxy`) says a generation is running while a run changes the galaxy, and shows the map
once the run has ended. Skipped without Playwright, Chromium or a Redis test server.
"""

import pytest

pytest.importorskip("playwright.sync_api")

from planetgen.queue import redisqueue, runrecord  # noqa: E402
from planetgen.queue.runrecord import Stage  # noqa: E402
from planetgen.web.app import create_app  # noqa: E402
from planetgen.web.lib import apiclient, galaxystatus, tilecache  # noqa: E402
from tests.served import served  # noqa: E402
from tests.test_web_a11y import browser  # noqa: E402,F401 -- fixture
from tests.test_web_galaxy import FakeData, _FakeConfig  # noqa: E402


@pytest.fixture(scope="module")
def base_url(tmp_path_factory):
    data = FakeData()
    patch = pytest.MonkeyPatch()
    for name in ("get_galaxy_sectors", "get_galaxy_shape", "get_polities", "get_superbubbles", "get_containers",
                 "auth_me", "get_galaxy_nebula_cover"):
        patch.setattr(apiclient, name, getattr(data, name))
    for name in ("get_galaxy_changes", "get_galaxy_nodes", "get_galaxy_stage", "get_galaxy_window"):
        patch.setattr(tilecache, name, getattr(data, name))
    patch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    patch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path_factory.mktemp("tiles")))
    app = create_app(_FakeConfig)
    app.testing = True
    try:
        with served(app) as url:
            yield url
    finally:
        patch.undo()


@pytest.fixture
def run_conn(redis_server, key_prefix, monkeypatch):
    connection = redisqueue.connect(redis_server)
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    monkeypatch.setattr(galaxystatus, "_client", lambda: connection)
    galaxystatus.clear_cache()
    yield connection
    galaxystatus.clear_cache()


def test_the_galaxy_map_is_closed_while_a_run_changes_the_galaxy_and_opens_when_it_ends(browser, base_url, run_conn):
    page = browser.new_page(viewport={"width": 1000, "height": 700})
    try:
        run_id = runrecord.start_run(run_conn, "Galaxy", [Stage("mass", "Scatter the massive stars", units=4)],
                                     database=_FakeConfig.WEB_DATABASE)
        runrecord.begin_stage(run_conn, run_id, 0)
        galaxystatus.clear_cache()
        page.goto(base_url + "/galaxy", wait_until="load")
        assert "closed while the galaxy is generated" in page.locator("#galaxy-generating").inner_text()
        assert page.locator("#galaxymap3d-canvas").count() == 0

        runrecord.finish_run(run_conn, run_id, "done")
        galaxystatus.clear_cache()
        page.goto(base_url + "/galaxy", wait_until="load")
        assert page.locator("#galaxy-generating").count() == 0
        page.locator("#galaxymap3d-canvas").wait_for(state="attached", timeout=10000)
    finally:
        page.close()
