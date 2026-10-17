"""
The loading state of a slow page in a real browser (UX.118): the Phenomena list while a queued job builds it.

The page is the real Flask view over the fake API of `test_web_browser_tables.py`; the fake API client raises
`slowpages.Pending` while the test says the answer is still being built, and the job's status is whatever the test
sets, so every state (queued, running, failed, cancelled, finished) is shown without a queue. The queue and the
job itself are in `test_slowpages.py` and `test_slow_route.py`; the polling logic in `tests/js/slowload.test.mjs`.
"""

import threading

import pytest
from flask import request
from tests.served import make_quiet_server
from tests.test_web_browser_tables import (  # noqa: F401 -- the shared fakes and fixtures
    DB, _SiteConfig, _fake_phenomena, _fake_sector, _fake_sectors, _fake_systems, browser,
)

sync_api = pytest.importorskip("playwright.sync_api")

from planetgen.queue import pagejobs  # noqa: E402
from planetgen.web import sector_page  # noqa: E402
from planetgen.web.app import create_app  # noqa: E402
from planetgen.web.lib import apiclient, slowpages  # noqa: E402

JOB_A, JOB_B = "0123456789abcdef", "fedcba9876543210"


def _job(state="running", done=300, total=1000, eta=5.0, error=None, job_id=JOB_A):
    return {"id": job_id, "state": state, "done": done, "total": total, "unit": "rows", "estimated": True,
            "eta_seconds": eta, "seconds": 3.0, "error": error, "cancelling": False,
            "cache_key": "planetgen:page:x:y", "target": f"/phenomena?db={DB}"}


class World:
    """What the fake API and the fake job say, changed by the tests."""

    def __init__(self):
        self.reset()

    def reset(self):
        self.ready = False
        self.job = _job()
        self.then_wait_for = None        # a job to wait for once the first has finished (a page with two answers)
        self.retried = []
        self.polls = 0


WORLD = World()


def _phenomena(db, **kwargs):
    WORLD.retried.append(request.headers.get("X-Slow-Retry") == "1")
    if not WORLD.ready or WORLD.then_wait_for is not None:
        job = WORLD.job if not WORLD.ready else WORLD.then_wait_for
        raise slowpages.Pending(job, f"/phenomena?db={db}")
    return _fake_phenomena(db, **kwargs)


def _status(job_id):
    WORLD.polls += 1
    for job in (WORLD.job, WORLD.then_wait_for):
        if job is not None and job["id"] == job_id:
            return dict(job)
    return None


@pytest.fixture(scope="module")
def slow_site():
    patch = pytest.MonkeyPatch()
    try:
        patch.setattr(apiclient, "get_phenomena", _phenomena)
        patch.setattr(apiclient, "get_sectors", _fake_sectors)
        patch.setattr(apiclient, "get_systems", _fake_systems)
        patch.setattr(apiclient, "get_uncharted_systems",
                      lambda db, limit=50, offset=0, **kwargs: {"items": [], "total": 0, "limit": limit, "offset": offset})
        patch.setattr(apiclient, "auth_me", lambda cookie_header: None)
        patch.setattr(apiclient, "get_population_status",
                      lambda db: {"generated": True, "species": True, "polities": True, "territories": False})
        patch.setattr(apiclient, "get_sector", _fake_sector)
        patch.setattr(apiclient, "get_galaxy_shape", lambda db: None)
        patch.setattr(sector_page, "fetch_tiles", lambda db, keys, known_stamp=None: {
            "stamp": "0" * 16, "tiles": {}, "density": {}, "edge_pc": 3.066, "has_shape": False})
        patch.setattr(pagejobs, "status", _status)
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
    context = browser.new_context(viewport={"width": width, "height": 900}, color_scheme=scheme)
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.errors = errors
    page.scheme, page.width = scheme, width
    try:
        yield page
    finally:
        context.close()
    assert not errors, f"uncaught exceptions: {errors}"


def _open(page, base_url):
    response = page.goto(f"{base_url}/phenomena", wait_until="load")
    assert response.status == 202
    page.locator("[data-slow-load]").wait_for(state="visible", timeout=15000)


def _status_line(page):
    return page.locator("[data-slow-status]").inner_text()


def test_a_running_job_shows_an_estimated_bar_that_follows_it_and_the_finished_list_replaces_it(page, slow_site):
    _open(page, slow_site)
    assert page.title().startswith("Phenomena")
    assert page.locator("h1").inner_text() == "Phenomena"
    bar = page.locator("[data-slow-bar]")
    assert bar.is_visible() and bar.get_attribute("value") == "30"
    assert page.locator("[data-slow-percent]").inner_text() == "30%"
    assert _status_line(page) == "About 30% (an estimate: 300 of 1,000 rows), about 5 s left."
    page.evaluate("window.__kept = 'same page'")

    WORLD.job = _job(done=700, eta=2.0)
    page.locator("[data-slow-percent]", has_text="70%").wait_for(state="attached", timeout=10000)
    assert "about 2 s left" in _status_line(page)

    WORLD.job, WORLD.ready = _job("succeeded", done=1000, eta=0.0), True
    page.locator("[data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)
    assert page.locator("[data-slow-load]").count() == 0
    assert page.evaluate("window.__kept") == "same page", "the list came in without a page load"
    assert page.locator(".datatable-count-line").inner_text() == "230 phenomena"
    assert page.locator("tbody tr[data-index]").first.locator("td").first.inner_text() == "Phenomenon 000"


def test_the_bar_fits_the_screen_and_the_text_is_readable(page, slow_site):
    _open(page, slow_site)
    assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
    box = page.locator("[data-slow-bar]").bounding_box()
    assert box["width"] >= 100 and box["x"] >= 0 and box["x"] + box["width"] <= page.width
    colour = page.evaluate(
        "(() => { const el = document.querySelector('[data-slow-status]'); "
        "const text = getComputedStyle(el).color; const back = getComputedStyle(document.body).backgroundColor; "
        "return [text, back]; })()")
    assert colour[0] != colour[1]


def test_a_queued_job_says_it_is_waiting(page, slow_site):
    WORLD.job = _job("queued", done=0, eta=None)
    _open(page, slow_site)
    assert page.locator("[data-slow-bar]").get_attribute("value") is None
    assert _status_line(page) == "Waiting for a worker to start building this page."
    assert page.locator("[data-slow-percent]").inner_text() == ""


def test_a_reader_who_asked_for_less_motion_gets_words_not_a_moving_bar(browser, slow_site):
    WORLD.reset()
    WORLD.job = _job("queued", done=0, eta=None)
    context = browser.new_context(viewport={"width": 1280, "height": 900}, reduced_motion="reduce")
    page = context.new_page()
    try:
        _open(page, slow_site)
        assert not page.locator("[data-slow-bar]").is_visible()
        assert _status_line(page).startswith("Waiting for a worker")
    finally:
        context.close()


def test_the_bar_is_a_progressbar_a_screen_reader_can_follow(page, slow_site):
    _open(page, slow_site)
    labelled_by = page.locator("progress[data-slow-bar]").get_attribute("aria-labelledby")
    assert page.locator(f"#{labelled_by}").inner_text() == "Building this page"
    status = page.locator("[data-slow-status]")
    assert status.get_attribute("role") == "status" and status.get_attribute("aria-live") == "polite"
    assert page.locator("[data-slow-error]").get_attribute("role") == "alert"
    assert page.get_by_role("progressbar").count() == 1


def test_a_failed_job_shows_its_error_and_try_again_starts_a_new_one(page, slow_site):
    WORLD.job = _job("failed", done=0, error="Lost connection to the database")
    _open(page, slow_site)
    error = page.locator("[data-slow-error]")
    error.wait_for(state="visible")
    assert error.inner_text() == "Lost connection to the database"
    assert page.locator("[data-slow-retry]").is_visible() and not page.locator("[data-slow-stop]").is_visible()
    WORLD.job = _job(job_id=JOB_B, done=100)
    page.locator("[data-slow-retry]").click()
    page.locator("[data-slow-percent]", has_text="10%").wait_for(state="attached", timeout=10000)
    assert WORLD.retried[-1] is True, "the retry asked for a new job"
    assert page.locator("[data-slow-load]").get_attribute("data-job") == JOB_B
    assert not error.is_visible()


def test_a_job_an_administrator_cancelled_says_so(page, slow_site):
    WORLD.job = _job("cancelled", done=0)
    _open(page, slow_site)
    page.locator("[data-slow-error]").wait_for(state="visible")
    assert page.locator("[data-slow-error]").inner_text() == "An administrator stopped building this page."
    assert page.locator("[data-slow-retry]").is_visible()


def test_stop_waiting_stops_the_polling_and_check_again_picks_it_up(page, slow_site):
    _open(page, slow_site)
    page.locator("[data-slow-stop]").click()
    assert _status_line(page).startswith("Stopped waiting.")
    polls = WORLD.polls
    page.wait_for_timeout(2500)
    assert WORLD.polls == polls, "a stopped page does not ask again"
    WORLD.job = _job(done=600)
    page.locator("[data-slow-resume]").click()
    page.locator("[data-slow-percent]", has_text="60%").wait_for(state="attached", timeout=10000)
    assert page.locator("[data-slow-stop]").is_visible()


def test_a_page_with_a_second_answer_to_wait_for_carries_on_with_that_job(page, slow_site):
    _open(page, slow_site)
    WORLD.job, WORLD.ready = _job("succeeded", done=1000, eta=0.0), True
    WORLD.then_wait_for = _job(job_id=JOB_B, done=500)
    page.locator("[data-slow-load][data-job='%s']" % JOB_B).wait_for(state="attached", timeout=15000)
    page.locator("[data-slow-percent]", has_text="50%").wait_for(state="attached", timeout=10000)
    WORLD.then_wait_for = None
    page.locator("[data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)


def test_without_scripts_the_page_refreshes_itself_until_it_is_ready(browser, slow_site):
    WORLD.reset()
    context = browser.new_context(viewport={"width": 1280, "height": 900}, java_script_enabled=False)
    page = context.new_page()
    try:
        _open_plain(page, slow_site)
        assert page.locator("h2", has_text="Building this page").is_visible()
        assert "refreshes itself every 5 seconds" in page.locator(".slow-load noscript").inner_text()
        WORLD.job, WORLD.ready = _job("succeeded", done=1000, eta=0.0), True
        page.wait_for_selector("[data-datatable]", state="attached", timeout=15000)
        assert page.locator("[data-slow-load]").count() == 0
    finally:
        context.close()


def _open_plain(page, base_url):
    assert page.goto(f"{base_url}/phenomena", wait_until="load").status == 202
