# tests/test_slow_route.py

"""
UX.118: a page whose answer is being built by a queued job (PERF.105) shows its frame with a progress bar
(`web/lib/slowroute.py`, `templates/slow_loading.html`) instead of the content, and the finished page replaces it.

The first group fakes the API client (no database): what the loading page contains and sends, the retry header,
the notice. The second runs the real chain over a small galaxy and the real queue: `/phenomena` and `/search`
answer `202` with a bar while their job is queued and `200` with the list once it has run.
"""

import pytest

from planetgen.api.config import Config
from planetgen.api.limiter import PAGE_LIMITS_OFF
from planetgen import tuning
from planetgen.db import dataversion, store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.queue import pagejobs, redisclient, workers
from planetgen.web.app import create_app
from planetgen.web import searchpage
from planetgen.web.lib import apiclient, slowpages, slowroute
from tests.test_api import seeded_sector  # noqa: F401 -- the fixture galaxy needs
from tests.test_cache_warmup import galaxy  # noqa: F401 -- the fixture
from tests.test_slowpages import _site, slow  # noqa: F401 -- the fixtures


def _job(state="running", done=300, total=1000, eta=5.0, error=None, job_id="0123456789abcdef"):
    return {"id": job_id, "state": state, "done": done, "total": total, "unit": "rows", "estimated": True,
            "eta_seconds": eta, "seconds": 3.0, "error": error, "cancelling": False,
            "cache_key": "planetgen:page:x:y", "target": "/phenomena?db=x"}


class _Site(Config):
    WEB_DATABASE = "fixture_db"
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "slow-route-secret"
    RATELIMIT_PAGES = PAGE_LIMITS_OFF


@pytest.fixture
def client(monkeypatch):
    monkeypatch.setattr(apiclient, "auth_me", lambda cookie_header: None)
    monkeypatch.setattr(apiclient, "get_population_status", lambda db: {"generated": True, "species": True, "polities": True, "territories": False})
    app = create_app(_Site)
    app.testing = True
    return app.test_client()


def _pending(job, seen=None):
    def fake(db, **kwargs):
        if seen is not None:
            seen.append(slowpages._retrying.get())
        raise slowpages.Pending(job, f"/phenomena?db={db}")
    return fake


# --- The loading page ----------------------------------------------------------------------------------

def test_a_waiting_page_is_the_pages_frame_with_a_bar_and_is_never_cached_or_indexed(client, monkeypatch):
    monkeypatch.setattr(apiclient, "get_phenomena", _pending(_job()))
    response = client.get("/phenomena")
    html = response.get_data(as_text=True)
    assert response.status_code == 202
    assert response.headers["Cache-Control"] == "no-store" and response.headers["X-Robots-Tag"] == "noindex"
    assert "<title>Phenomena" in html and '<h1 class="page-title">Phenomena</h1>' in html
    assert 'data-job="0123456789abcdef"' in html and 'data-status-url="/api/page-jobs/0123456789abcdef"' in html
    assert '<progress data-slow-bar max="100"' in html and 'value="30"' in html
    assert "About 30% (an estimate: 300 of 1,000 rows)" in html and "about 5 s left" in html
    assert '<meta name="robots" content="noindex">' in html
    assert '<noscript><meta http-equiv="refresh" content="5"></noscript>' in html
    assert "slowload.js" in html and "Stop waiting" in html


def test_a_queued_job_has_a_bar_with_no_value_and_says_it_is_waiting(client, monkeypatch):
    monkeypatch.setattr(apiclient, "get_phenomena", _pending(_job("queued", done=0)))
    html = client.get("/phenomena").get_data(as_text=True)
    assert "<progress data-slow-bar max=\"100\"" in html and "value=" not in html.split("<progress")[1].split(">")[0]
    assert "Waiting for a worker to start building this page." in html


@pytest.mark.parametrize("state,words", [("failed", "Disk full"), ("cancelled", "An administrator stopped")])
def test_a_job_that_ended_badly_shows_why_and_offers_to_try_again(client, monkeypatch, state, words):
    monkeypatch.setattr(apiclient, "get_phenomena", _pending(_job(state, done=0, error="Disk full")))
    html = client.get("/phenomena").get_data(as_text=True)
    assert words in html
    assert '<p class="error" role="alert" data-slow-error>' in html
    assert "data-slow-retry>Try again" in html and "data-slow-stop hidden" in html.replace(">Stop waiting", "")


def test_try_again_asks_for_a_new_job(client, monkeypatch):
    seen = []
    monkeypatch.setattr(apiclient, "get_phenomena", _pending(_job(), seen))
    client.get("/phenomena")
    client.get("/phenomena", headers={slowroute.RETRY_HEADER: "1"})
    assert seen == [False, True]


def test_the_search_page_waits_the_same_way(client, monkeypatch):
    def pending(db, *args, **kwargs):
        raise slowpages.Pending(_job(), "/search?db=x")
    monkeypatch.setattr(apiclient, "get_search", pending)
    response = client.get("/search?q=sol")
    assert response.status_code == 202
    html = response.get_data(as_text=True)
    assert '<h1 class="page-title">Search</h1>' in html and "data-slow-load" in html


def test_a_page_that_built_a_slow_answer_itself_says_why(client, monkeypatch):
    def fake(db, **kwargs):
        slowpages._notices.text = "This page is slow to build and the work queue is not available, so it was built while you waited."
        return {"items": [], "total": 0, "limit": 50, "offset": 0, "facets": {"type": [], "descriptor": []}}
    monkeypatch.setattr(apiclient, "get_phenomena", fake)
    html = client.get("/phenomena").get_data(as_text=True)
    assert 'class="hint slow-notice" role="status">This page is slow to build' in html
    assert "slow-notice" not in client.get("/").get_data(as_text=True)         # said once


def test_a_notice_left_on_the_thread_is_not_the_next_pages(client, monkeypatch):
    slowpages._notices.text = "left over"
    monkeypatch.setattr(apiclient, "get_phenomena", lambda db, **kwargs: {
        "items": [], "total": 0, "limit": 50, "offset": 0, "facets": {"type": [], "descriptor": []}})
    assert "left over" not in client.get("/phenomena").get_data(as_text=True)


# --- The words and the bar's value ------------------------------------------------------------------------

def test_the_line_under_the_bar_is_an_estimate_that_says_so():
    assert slowroute.progress_text(_job()) == "About 30% (an estimate: 300 of 1,000 rows), about 5 s left."
    assert slowroute.progress_text(_job(eta=None)).endswith("still working.")
    assert slowroute.progress_text(_job("queued", done=0)) == slowroute.STATE_TEXT["queued"]
    assert slowroute.progress_text(_job(done=999, total=1000)).startswith("About 99%")        # never 100 before it is built
    assert slowroute.progress_text(_job(total=None, done=None)) == slowroute.STATE_TEXT["running"]


def test_the_bar_has_no_value_until_it_can_say():
    assert slowroute.percent(_job()) == 30 and slowroute.percent(_job("succeeded", done=1000)) == 100
    assert slowroute.percent(_job("queued", done=0)) is None and slowroute.percent(_job(total=None, done=None)) is None
    assert slowroute.percent(_job(done=1000)) == 99


# --- The real chain ----------------------------------------------------------------------------------------

def test_a_waiting_page_queues_one_job_for_every_visitor_who_asks(slow):  # noqa: F811
    config = slow.config
    client = _site(config).test_client()
    waiting = client.get("/phenomena")
    assert waiting.status_code == 202 and "data-slow-load" in waiting.get_data(as_text=True)
    [job_id] = _queued_ids(config.database)
    assert f'data-job="{job_id}"' in waiting.get_data(as_text=True)
    assert _site(config).test_client().get("/phenomena").status_code == 202    # another process: the same job
    assert _queued_ids(config.database) == [job_id]
    search = client.get("/search?q=sol")
    assert search.status_code == 202 and 'class="page-title">Search<' in search.get_data(as_text=True)


def test_the_finished_page_is_served_once_the_job_has_run(slow, monkeypatch):  # noqa: F811
    config = slow.config
    monkeypatch.setattr(workers, "ensure_workers", slow.run_jobs)           # the tests' inline worker runs the job as it is queued
    real_ensure = pagejobs.ensure
    monkeypatch.setattr(pagejobs, "ensure", lambda db, target, rows=None, expected=None, config=None, retry=False:
                        real_ensure(db, target, rows, expected, slow.config, retry))
    client = _site(config).test_client()
    first = client.get("/phenomena")
    assert first.status_code == 202                                         # this request found the answer missing
    page = client.get("/phenomena")
    html = page.get_data(as_text=True)
    assert page.status_code == 200 and "data-datatable" in html and "data-slow-load" not in html


def test_the_warm_up_builds_a_slow_page_in_its_own_request_and_never_gets_the_loading_page(slow):  # noqa: F811
    """PERF.156: every request of the galaxy counts as slow here; a visitor is told to wait (202), the warm-up's visit
    is answered with the page (200) and nothing is queued for it."""
    from planetgen.web import warmup
    config = slow.config
    assert _site(config).test_client().get("/search").status_code == 202
    visitor = warmup._Visitor(warmup.build_app(config), config.database)
    assert visitor.visit("/search").status_code == 200 and visitor.visit("/phenomena").status_code == 200
    assert visitor.failed == []


def _add_an_m_class_planet(config):
    """A sector with a planet of class M in the fixture galaxy (which has no planets); returns its name."""
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.MOONS = False
    cfg.BINARY_SYSTEM = False
    for _ in range(60):
        system = StarSystem(system_config=cfg)
        if any(planet.body_type != "a" for planet in system.planets):
            break
    sector = SpaceSector("Mclass Field", edge_ly=10.0)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
    store.save_sector(sector, config=config)
    conn = store.get_connection(config)
    try:
        planet_id = conn.execute("SELECT MIN(id) AS id FROM planets WHERE body_type <> 'a'").fetchone()["id"]
        conn.execute("UPDATE planets SET planet_class = 'M', name = 'Emmaline Prime' WHERE id = ?", (planet_id,))
        conn.commit()
    finally:
        conn.close()
    return "Emmaline Prime"


def _every_build_is_slow(monkeypatch, slow):
    """The cost rule calls every request slow whatever it measured, and the tests' worker runs a job as it is queued."""
    monkeypatch.setattr(tuning, "SLOW_PAGE_SECONDS", -1)
    monkeypatch.setattr(workers, "ensure_workers", slow.run_jobs)
    real_ensure = pagejobs.ensure
    monkeypatch.setattr(pagejobs, "ensure", lambda db, target, rows=None, expected=None, config=None, retry=False:
                        real_ensure(db, target, rows, expected, slow.config, retry))


def _until_it_answers(client, address, most=12):
    """The status of each request for `address` until one is not the loading page (at most `most`), and the last page."""
    codes = []
    for _ in range(most):
        response = client.get(address)
        codes.append(response.status_code)
        if response.status_code != 202:
            break
    return codes, response.get_data(as_text=True)


def test_a_search_for_m_class_planets_asks_for_all_its_answers_before_it_waits(slow):  # noqa: F811
    """UX.234: the search page's eight answers (the facets and each result panel) are queued together, so the wait is
    the slowest build and not the eight in turn, which each went stale before the last was done."""
    _add_an_m_class_planet(slow.config)
    client = _site(slow.config).test_client()
    assert client.get("/search?class=M").status_code == 202
    assert len(_queued_ids(slow.config.database)) == 1 + len(searchpage.PANEL_NAMES)


def test_a_search_for_m_class_planets_ends_with_rows_even_while_the_data_keeps_changing(slow, monkeypatch):  # noqa: F811
    """UX.234: the data changes under every build (a fill is running), so the page cache refuses each finished answer;
    the page that waited is still answered with its rows instead of being sent to queue the same jobs again."""
    name = _add_an_m_class_planet(slow.config)
    _every_build_is_slow(monkeypatch, slow)
    real_send = apiclient._send

    def send(*args, **kwargs):
        answer = real_send(*args, **kwargs)
        if args[0] == "GET" and "/search" in args[1]:
            dataversion.bump(slow.config.database, dataversion.TAG_SECTOR_DATA)
        return answer
    monkeypatch.setattr(apiclient, "_send", send)
    codes, html = _until_it_answers(_site(slow.config).test_client(), "/search?class=M")
    assert codes[-1] == 200 and len(codes) <= 3, codes
    assert name in html and "data-slow-load" not in html


def test_a_search_for_m_class_planets_ends_with_rows_when_redis_has_no_room_for_the_cache(slow, monkeypatch):  # noqa: F811
    """UX.234: the page cache stops adding answers once Redis is nearly full; the finished jobs hand theirs over."""
    from planetgen.web.lib import pagecache
    name = _add_an_m_class_planet(slow.config)
    _every_build_is_slow(monkeypatch, slow)
    monkeypatch.setattr(pagecache.ResponseCache, "_redis_has_room", lambda self: False)
    codes, html = _until_it_answers(_site(slow.config).test_client(), "/search?class=M")
    assert codes[-1] == 200 and len(codes) <= 3, codes
    assert name in html


def test_a_job_that_the_cache_took_hands_nothing_over(slow, monkeypatch):  # noqa: F811
    """The answer is kept for the visitors only when the page cache would not take it."""
    _every_build_is_slow(monkeypatch, slow)
    client = _site(slow.config).test_client()
    _until_it_answers(client, "/phenomena")
    connection = redisclient.client(decode_responses=True)
    assert list(connection.scan_iter(f"planetgen:slowanswer:{slow.config.database}:*")) == []


def _queued_ids(db):
    connection = redisclient.client(decode_responses=True)
    found = []
    for key in connection.scan_iter("planetgen:slowmeta:*"):
        meta = connection.hgetall(key)
        if meta.get("db") == db and meta.get("queue_job"):
            found.append(key.rpartition(":")[2])
    return found
