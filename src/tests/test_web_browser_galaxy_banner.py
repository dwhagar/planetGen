"""
UX.3: the galaxy banner (`static/galaxybanner.js`) appears on an open page when a background run starts, goes
when it ends, and needs no reload. Skipped without Playwright, Chromium or a Redis test server.
"""

import pytest

pytest.importorskip("playwright.sync_api")

from planetgen.queue import redisqueue, runrecord  # noqa: E402
from planetgen.queue.runrecord import Stage  # noqa: E402
from planetgen.web.app import create_app  # noqa: E402
from planetgen.web.lib import apiclient, galaxystatus  # noqa: E402
from tests.served import served  # noqa: E402
from tests.test_web_a11y import browser  # noqa: E402,F401 -- fixture
from tests.test_web_pages import FakeData, _FakeConfig  # noqa: E402


@pytest.fixture
def run_conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture(scope="module")
def base_url():
    data = FakeData()
    patch = pytest.MonkeyPatch()
    for name in ("get_sectors", "get_systems", "get_uncharted_systems", "get_phenomena", "auth_me"):
        patch.setattr(apiclient, name, getattr(data, name))
    app = create_app(_FakeConfig)
    app.testing = True
    try:
        with served(app) as url:
            yield url
    finally:
        patch.undo()


@pytest.fixture
def prefixed(run_conn, key_prefix, monkeypatch):
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    monkeypatch.setattr(galaxystatus, "_client", lambda: run_conn)
    galaxystatus.clear_cache()
    yield
    galaxystatus.clear_cache()


def _ask_again(page):
    """The tab coming back to the front asks at once, which is what the 30 second timer does too."""
    page.evaluate("() => document.dispatchEvent(new Event('visibilitychange'))")


def test_the_banner_appears_and_goes_without_a_reload(browser, base_url, run_conn, prefixed):
    page = browser.new_page(viewport={"width": 1000, "height": 700})
    try:
        page.goto(base_url + "/", wait_until="load")
        banner = page.locator(".galaxy-banner")
        assert banner.count() == 1 and not banner.is_visible()

        run_id = runrecord.start_run(run_conn, "Fill shell 1",
                                     [Stage("mass", "Mass scatter", units=4, cost=8.0, expected_seconds=7200)],
                                     database=_FakeConfig.WEB_DATABASE)
        runrecord.begin_stage(run_conn, run_id, 0)
        galaxystatus.clear_cache()
        _ask_again(page)
        banner.wait_for(state="visible", timeout=10000)
        assert "The galaxy is being changed" in banner.inner_text()
        assert "about 3 hours" in banner.inner_text()

        runrecord.finish_run(run_conn, run_id, "done")
        galaxystatus.clear_cache()
        _ask_again(page)
        banner.wait_for(state="hidden", timeout=10000)
    finally:
        page.close()
