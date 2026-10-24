# tests/test_web_browser_slow_tables.py

"""
UX.136: a data table's own row fetches (scroll and sort) wait on a page job when the data route answers 202, in a real
browser: the Phenomena list over the fake API of `test_web_browser_tables.py`, where the fake raises
`slowpages.Pending` for the table's JSON route while the test says the rows are still being built. The job's status is
whatever the test sets. The route's answer is in `test_slow_route.py`; the polling logic in
`tests/js/datatablestate.test.mjs`.
"""

import threading

import pytest
from flask import request
from tests.served import make_quiet_server
from tests.test_web_browser_tables import (  # noqa: F401 -- the shared fakes and fixtures
    DB, TOTAL, _SiteConfig, _fake_phenomena, _fake_sector, _fake_sectors, _fake_systems, browser,
)

sync_api = pytest.importorskip("playwright.sync_api")

from planetgen.queue import pagejobs  # noqa: E402
from planetgen.web import sector_page  # noqa: E402
from planetgen.web.app import create_app  # noqa: E402
from planetgen.web.lib import apiclient, slowpages  # noqa: E402

JOB = "0123456789abcdef"


def _job(state="running", done=300, total=1000, eta=5.0, error=None):
    return {"id": JOB, "state": state, "done": done, "total": total, "unit": "rows", "estimated": True,
            "eta_seconds": eta, "seconds": 3.0, "error": error, "cancelling": False,
            "cache_key": "planetgen:page:x:y", "target": f"/phenomena?db={DB}"}


class World:
    """What the fake API and the fake job say, changed by the tests."""

    def __init__(self):
        self.reset()

    def reset(self):
        self.slow = False       # the table route is waiting on the job
        self.job = _job()
        self.retried = []


WORLD = World()


def _phenomena(db, **kwargs):
    if request.path.startswith("/table/"):
        WORLD.retried.append(request.headers.get("X-Slow-Retry") == "1")
        if WORLD.slow:
            raise slowpages.Pending(WORLD.job, f"/phenomena?db={db}")
    return _fake_phenomena(db, **kwargs)


@pytest.fixture(scope="module")
def slow_site():
    patch = pytest.MonkeyPatch()
    try:
        patch.setattr(apiclient, "get_phenomena", _phenomena)
        patch.setattr(apiclient, "get_sectors", _fake_sectors)
        patch.setattr(apiclient, "get_systems", _fake_systems)
        patch.setattr(apiclient, "auth_me", lambda cookie_header: None)
        patch.setattr(apiclient, "get_population_status",
                      lambda db: {"generated": True, "species": True, "polities": True, "territories": False})
        patch.setattr(apiclient, "get_sector", _fake_sector)
        patch.setattr(apiclient, "get_galaxy_shape", lambda db: None)
        patch.setattr(sector_page, "fetch_nodes", lambda db, keys, known_stamp=None, ahead=False: {
            "stamp": "0" * 16, "nodes": {}, "edge_pc": 3.066, "has_shape": False})
        patch.setattr(pagejobs, "status", lambda job_id: dict(WORLD.job) if job_id == JOB else None)
        app = create_app(_SiteConfig)
        app.testing = True
        server = make_quiet_server("127.0.0.1", 0, app, threaded=True)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            yield f"http://127.0.0.1:{server.server_port}"
        finally:
            server.shutdown()
            thread.join(timeout=10)
    finally:
        patch.undo()


@pytest.fixture(params=[("dark", 1280), ("light", 1280), ("dark", 390), ("light", 390)],
                ids=["dark", "light", "dark-phone", "light-phone"])
def page(browser, request):  # noqa: F811
    scheme, width = request.param
    WORLD.reset()
    context = browser.new_context(viewport={"width": width, "height": 900}, color_scheme=scheme, reduced_motion="reduce")
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.errors = errors
    page.width = width
    try:
        yield page
    finally:
        context.close()
    assert not errors, f"uncaught exceptions: {errors}"


def _open(page, base_url):
    assert page.goto(f"{base_url}/phenomena", wait_until="load").status == 200
    page.locator("[data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)


def _sort_by_name(page):
    page.locator("th[data-col='name'] .datatable-sort").click()


def _panel(page):
    return page.locator("[data-table-job]")


def test_a_slow_sort_shows_the_jobs_progress_and_the_sorted_rows_replace_the_old_ones_when_it_is_built(page, slow_site):
    _open(page, slow_site)
    WORLD.slow = True
    _sort_by_name(page)
    _panel(page).wait_for(state="visible", timeout=10000)
    assert page.locator("[data-job-bar]").get_attribute("value") == "30"
    assert page.locator("[data-job-percent]").inner_text() == "30%"
    assert page.locator("[data-job-status]").inner_text() == "About 30% (an estimate), about 5 s left."
    assert page.locator("[data-job-stop]").is_visible() and not page.locator("[data-job-retry]").is_visible()
    assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth"), "no sideways scroll"

    WORLD.job = _job(done=700, eta=2.0)
    page.locator("[data-job-percent]", has_text="70%").wait_for(state="attached", timeout=10000)
    WORLD.slow = False
    WORLD.job = _job("succeeded")
    _panel(page).wait_for(state="hidden", timeout=10000)
    page.locator("tbody tr[data-index='0']", has_text=f"Phenomenon {TOTAL - 1:03d}").wait_for(state="attached", timeout=10000)


def test_a_slow_scroll_waits_on_the_job_too(page, slow_site):
    _open(page, slow_site)
    WORLD.slow = True
    WORLD.job = _job("queued", done=0)
    page.locator(".datatable-scroll").evaluate("el => { el.scrollTop = el.scrollHeight; }")
    _panel(page).wait_for(state="visible", timeout=10000)
    assert page.locator("[data-job-status]").inner_text() == "Waiting for a worker to start building these rows."
    assert page.locator("[data-job-bar]").get_attribute("value") is None
    WORLD.slow = False
    WORLD.job = _job("succeeded")
    _panel(page).wait_for(state="hidden", timeout=10000)
    page.locator("tbody tr[data-index]", has_text=f"Phenomenon {TOTAL - 1:03d}").wait_for(state="attached", timeout=10000)


def test_a_failed_job_says_why_and_try_again_asks_for_a_new_one(page, slow_site):
    _open(page, slow_site)
    WORLD.slow = True
    WORLD.job = _job("failed", done=0, error="Disk full")
    _sort_by_name(page)
    page.locator("[data-job-error]").wait_for(state="visible", timeout=10000)
    assert page.locator("[data-job-error]").inner_text() == "Disk full"
    assert page.locator("[data-job-retry]").is_visible() and not page.locator("[data-job-stop]").is_visible()
    WORLD.slow = False
    page.locator("[data-job-retry]").click()
    _panel(page).wait_for(state="hidden", timeout=10000)
    assert True in WORLD.retried, "Try again asked for a new job"
    page.locator("tbody tr[data-index='0']", has_text=f"Phenomenon {TOTAL - 1:03d}").wait_for(state="attached", timeout=10000)


def test_stop_waiting_ends_the_polling_and_check_again_picks_it_up(page, slow_site):
    _open(page, slow_site)
    WORLD.slow = True
    _sort_by_name(page)
    _panel(page).wait_for(state="visible", timeout=10000)
    page.locator("[data-job-stop]").click()
    page.locator("[data-job-resume]").wait_for(state="visible", timeout=10000)
    assert page.locator("[data-job-status]").inner_text().startswith("Stopped waiting.")
    assert not page.locator("[data-job-stop]").is_visible()
    WORLD.slow = False
    page.locator("[data-job-resume]").click()
    _panel(page).wait_for(state="hidden", timeout=10000)
    page.locator("tbody tr[data-index='0']", has_text=f"Phenomenon {TOTAL - 1:03d}").wait_for(state="attached", timeout=10000)
