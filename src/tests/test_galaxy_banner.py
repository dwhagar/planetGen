"""UX.3: the banner every visitor sees while a background run changes the galaxy.

The whole-hour rule (`queue/hourdisplay.py`), the public status read from the run record in Redis
(`web/lib/galaxystatus.py`), `/galaxy-status`, and the banner in `base.html`.
"""

import random
import time

import pytest

from planetgen.queue import hourdisplay, progress_rate, redisqueue, runrecord
from planetgen.queue.runrecord import Stage
from planetgen.web.app import create_app
from planetgen.web.lib import apiclient, galaxystatus
from tests.test_web_pages import FakeData, _FakeConfig

HOUR = 3600.0


# --- The hour rule ---------------------------------------------------------------

def test_a_fresh_figure_is_the_top_of_the_range_rounded_up_and_never_below_one_hour():
    assert hourdisplay.hours_left(0.0) == 1
    assert hourdisplay.hours_left(60.0) == 1
    assert hourdisplay.hours_left(HOUR) == 2                      # 1.15 h rounds up to 2
    assert hourdisplay.hours_left(10 * HOUR) == 12                # 11.5 h
    assert hourdisplay.hours_left(2.6 * HOUR) == 3                # 2.99 h


def test_the_figure_falls_freely_and_rises_only_when_even_the_low_end_is_longer():
    shown = hourdisplay.hours_left(10 * HOUR)
    assert shown == 12
    assert hourdisplay.hours_left(7 * HOUR, shown) == 9           # falls at once
    assert hourdisplay.hours_left(10 * HOUR, 9) == 9              # 0.88 x 10 h = 8.8 h is not over 9: stays
    assert hourdisplay.hours_left(11 * HOUR, 9) == 13             # 9.68 h is over 9: rises to the fresh figure


def test_a_noisy_estimate_seldom_raises_the_figure_and_rarely_under_promises():
    rng = random.Random(3)
    rises = under = steps = 0
    for _ in range(200):
        truth, shown = 10 * HOUR, None
        for step in range(100):
            left = truth * (1.0 - step / 100.0)
            guess = left * rng.uniform(progress_rate.ETA_LOW_FACTOR, progress_rate.ETA_HIGH_FACTOR)
            new = hourdisplay.hours_left(guess, shown)
            rises += shown is not None and new > shown
            under += new * HOUR < left
            steps += 1
            shown = new
    assert rises / 200 < 1.0
    assert under / steps < 0.05


# --- The status read from Redis --------------------------------------------------

@pytest.fixture
def conn(redis_server, key_prefix, monkeypatch):
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    connection = redisqueue.connect(redis_server)
    monkeypatch.setattr(galaxystatus, "_client", lambda: connection)
    galaxystatus.clear_cache()
    yield connection
    galaxystatus.clear_cache()


def _start(conn, database="g1", now=100.0):
    run_id = runrecord.start_run(conn, "Fill shell 1", [Stage("mass", "Mass scatter", units=4, cost=8.0,
                                                              expected_seconds=3 * HOUR)],
                                 who="boss", database=database, argv=["secret"], now=now)
    runrecord.begin_stage(conn, run_id, 0, now=now)
    return run_id


def test_no_run_means_no_banner(conn):
    assert galaxystatus.compute(conn, "g1") == {"active": False}


def test_a_run_shows_hours_from_earlier_runs_until_it_has_its_own_rate(conn):
    _start(conn)
    status = galaxystatus.compute(conn, "g1", now=160.0)
    assert status["active"] and status["hours"] == 4              # 3 h expected, 1.15 x, rounded up
    assert status["stage"] == "Mass scatter"
    assert "about 4 hours (a guess from earlier runs)" in status["text"]
    assert "boss" not in repr(status) and "secret" not in repr(status)   # nothing private


def test_the_figure_is_shared_and_only_falls_while_the_estimate_falls(conn):
    run_id = _start(conn)
    first = galaxystatus.compute(conn, "g1", now=160.0)["hours"]
    assert conn.get(runrecord.banner_key(run_id)) == str(first).encode()
    conn.hset(runrecord.run_key(run_id), "updated", 100.0 + HOUR)          # a unit finished: the run is alive
    later = galaxystatus.compute(conn, "g1", now=160.0 + HOUR)["hours"]
    assert later < first                                           # an hour of the 3 h expected has gone


def test_a_run_for_another_database_or_a_finished_one_shows_nothing(conn):
    run_id = _start(conn, database="other")
    assert galaxystatus.compute(conn, "g1", now=160.0) == {"active": False}
    runrecord.finish_run(conn, run_id, "done", now=200.0)
    assert galaxystatus.compute(conn, "other", now=210.0) == {"active": False}


def test_a_run_with_no_sign_of_life_shows_nothing(conn):
    _start(conn, now=100.0)
    assert galaxystatus.compute(conn, "g1", now=100.0 + 3 * HOUR) == {"active": False}


def test_a_stage_with_no_estimate_says_it_is_working_it_out(conn):
    run_id = runrecord.start_run(conn, "Fill", [Stage("mass", "Mass scatter", units=4)], database="g1", now=100.0)
    runrecord.begin_stage(conn, run_id, 0, now=100.0)
    status = galaxystatus.compute(conn, "g1", now=101.0)
    assert status["active"] and status["hours"] is None
    assert "Working out how long it will take." in status["text"]


def test_redis_down_means_no_banner_not_an_error(monkeypatch):
    import redis

    def broken():
        raise redis.ConnectionError("down")

    monkeypatch.setattr(galaxystatus, "_client", broken)
    galaxystatus.clear_cache()
    assert galaxystatus.status("g1") == {"active": False}
    galaxystatus.clear_cache()


def test_the_answer_is_kept_for_a_few_seconds(conn, monkeypatch):
    _start(conn)
    first = galaxystatus.status("g1")
    runrecord.finish_run(conn, runrecord.active_run(conn), "done")
    assert galaxystatus.status("g1") == first                      # still the kept answer
    galaxystatus.clear_cache()
    assert galaxystatus.status("g1") == {"active": False}


# --- The page and the endpoint ----------------------------------------------------

@pytest.fixture
def client(conn, monkeypatch):
    data = FakeData()
    for name in ("get_sectors", "get_systems", "get_uncharted_systems", "get_phenomena", "auth_me"):
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    app = create_app(_FakeConfig)
    app.testing = True
    return app.test_client()


def test_every_page_has_a_hidden_banner_when_idle_and_the_status_says_so(client):
    html = client.get("/").get_data(as_text=True)
    assert 'class="galaxy-banner"' in html and "galaxybanner.js" in html
    banner = html[html.index('class="galaxy-banner"'):].split("</p>")[0]
    assert " hidden" in banner
    reply = client.get("/galaxy-status")
    assert reply.status_code == 200 and reply.get_json() == {"active": False}
    assert reply.headers["Cache-Control"] == f"public, max-age={int(galaxystatus.CACHE_SECONDS)}"


def test_while_a_run_is_going_every_visitor_sees_the_banner_and_the_status(client, conn):
    _start(conn, database=_FakeConfig.WEB_DATABASE, now=time.time())
    galaxystatus.clear_cache()
    for path in ("/", "/sectors"):
        html = client.get(path).get_data(as_text=True)
        banner = html[html.index('class="galaxy-banner"'):].split("</p>")[0]
        assert " hidden" not in banner and "The galaxy is being changed" in banner and "about" in banner
    body = client.get("/galaxy-status").get_json()
    assert body["active"] and body["hours"] >= 1 and "text" in body
