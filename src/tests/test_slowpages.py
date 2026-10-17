# tests/test_slowpages.py

"""
PERF.105: a page request that is slow to build queues one job that builds the answer into the tagged page cache
(`planetgen/queue/pagejobs.py`, `planetgen/web/lib/slowpages.py`), with a status route, a node on the admin runs
page, a cancel, and a fallback to building it in the request when the queue can't take it.

The jobs are real RQ jobs on the test Redis. Most tests keep them queued (no worker runs), which is enough to show
how many were queued and what the status says; two let the tests' own inline worker (conftest `_run_queue_inline`)
run them as they are queued, end to end.
"""

import threading
import time

import pytest

from planetgen import tuning
from planetgen.admin import auth as adminAuth
from planetgen.db import dataversion, query, store
from planetgen.queue import api_jobs, pagejobs, redisclient, work as workQueue, workers
from planetgen.web.app import create_app
from planetgen.web.lib import apiclient, pagecache, slowpages
from tests.test_api import seeded_sector  # noqa: F401 -- the fixture galaxy needs
from tests.test_cache_warmup import _Visitor, galaxy  # noqa: F401 -- the fixture


def _site(config):
    """A web process of its own on `config`: its own in-process copies, the shared Redis."""

    class Site(_Visitor):
        MYSQL_CONFIG = WRITE_MYSQL_CONFIG = CONTROL_MYSQL_CONFIG = config
        WEB_DATABASE = ""

    app = create_app(Site)
    app.testing = True
    return app


def _list(app, db):
    """What the phenomena page asks for, as a page that opted in."""
    with app.test_request_context("/"), slowpages.allow_pending():
        return apiclient.get_phenomena(db)


def _slow_jobs(db):
    """The page jobs recorded for `db` that were queued: their ids."""
    connection = redisclient.client(decode_responses=True)
    found = []
    for key in connection.scan_iter("planetgen:slowmeta:*"):
        meta = connection.hgetall(key)
        if meta.get("db") == db and meta.get("queue_job"):
            found.append(key.rpartition(":")[2])
    return found


@pytest.fixture
def slow(galaxy, monkeypatch):  # noqa: F811
    """Every request on the phenomena list counts as slow, and no worker starts on its own."""
    monkeypatch.setattr(tuning, "SLOW_PAGE_ROW_BUDGET", -1)
    galaxy.run_jobs = workers.ensure_workers       # the tests' own (conftest `_run_queue_inline`): runs waiting jobs here
    monkeypatch.setattr(workers, "ensure_workers", lambda *args, **kwargs: 0)
    slowpages.clear_estimates()
    yield galaxy
    connection = redisclient.client(decode_responses=True)
    for key in list(connection.scan_iter("planetgen:slow*")):
        meta = connection.hgetall(key) if key.startswith("planetgen:slowmeta:") else {}
        if meta.get("db") == galaxy.config.database and meta.get("queue_job"):
            api_jobs.cancel(meta["queue_job"])
        if key.startswith(("planetgen:slowmeta:", f"planetgen:slowjob:{galaxy.config.database}:",
                           f"planetgen:slowtime:{galaxy.config.database}:",
                           f"planetgen:slowpath:{galaxy.config.database}:")):
            connection.delete(key)
    slowpages.clear_estimates()


# --- The cost rule -----------------------------------------------------------------------------------

def test_a_request_over_the_row_budget_is_slow_and_one_under_it_is_not(slow, monkeypatch):
    db = slow.config.database
    app = _site(slow.config)
    target = f"/phenomena?db={db}"
    with app.test_request_context("/"):
        assert slowpages.cost(db, "/phenomena", target)["why"] == "rows"
        monkeypatch.setattr(tuning, "SLOW_PAGE_ROW_BUDGET", 10 ** 12)
        assert slowpages.cost(db, "/phenomena", target) is None
        assert slowpages.cost(db, "/sectors/5", "/sectors/5") is None      # not a route that reads a table in full


def test_a_measured_time_beats_the_row_estimate(slow, monkeypatch):
    db = slow.config.database
    app = _site(slow.config)
    target = f"/phenomena?db={db}"
    with app.test_request_context("/"):
        pagejobs.record_duration(db, target, tuning.SLOW_PAGE_SECONDS / 10)
        assert slowpages.cost(db, "/phenomena", target) is None           # over budget by rows, but quick in fact
        pagejobs.record_duration(db, target, tuning.SLOW_PAGE_SECONDS * 3)
        assert slowpages.cost(db, "/phenomena", target) == {
            "rows": slowpages.estimated_rows(db, "/phenomena"), "expected": tuning.SLOW_PAGE_SECONDS * 3,
            "why": "timed"}
        monkeypatch.setattr(tuning, "SLOW_PAGE_ROW_BUDGET", 10 ** 12)
        assert slowpages.cost(db, "/phenomena", target)["why"] == "timed"  # slow in fact, small by rows


def test_a_request_not_seen_before_takes_its_route_s_recent_slowest_build(slow, monkeypatch):
    db = slow.config.database
    app = _site(slow.config)
    other = f"/phenomena?db={db}&limit=5"
    with app.test_request_context("/"):
        assert slowpages.cost(db, "/phenomena", other)["why"] == "rows"      # nothing on the route timed yet
        pagejobs.record_duration(db, f"/phenomena?db={db}", tuning.SLOW_PAGE_SECONDS * 3)
        assert slowpages.cost(db, "/phenomena", other)["why"] == "route"     # a wide one was slow a moment ago
        redisclient.client().delete(pagejobs.path_key(db, "/phenomena"))
        pagejobs.record_duration(db, f"/phenomena?db={db}", tuning.SLOW_PAGE_SECONDS / 10)
        assert slowpages.cost(db, "/phenomena", other) is None                # the route has been quick lately


def test_the_phenomena_route_is_costed_on_every_table_it_reads():
    listed = set(slowpages.PATH_TABLES["/phenomena"])
    assert {table for table, *_rest in query._PHENOMENON_TABLES} | {"phenomenon_scatter"} <= listed


# --- One request, one job ------------------------------------------------------------------------------

def test_a_page_that_did_not_opt_in_builds_the_answer_itself(slow):
    app = _site(slow.config)
    with app.test_request_context("/"):
        assert apiclient.get_phenomena(slow.config.database)["total"] >= 1
    assert _slow_jobs(slow.config.database) == []


def test_a_slow_request_queues_a_job_and_a_repeat_attaches_to_it(slow):
    db = slow.config.database
    app = _site(slow.config)
    with pytest.raises(slowpages.Pending) as first:
        _list(app, db)
    assert first.value.job["state"] == "queued" and first.value.target == f"/phenomena?db={db}"
    with pytest.raises(slowpages.Pending) as again:
        _list(_site(slow.config), db)                                      # another web process, same request
    assert again.value.job["id"] == first.value.job["id"]
    assert _slow_jobs(db) == [first.value.job["id"]]


def test_visitors_asking_at_the_same_moment_share_one_job(slow):
    db = slow.config.database
    apps = [_site(slow.config) for _ in range(4)]
    ids, start = [], threading.Barrier(len(apps))

    def visit(app):
        start.wait()
        try:
            _list(app, db)
        except slowpages.Pending as pending:
            ids.append(pending.job["id"])

    threads = [threading.Thread(target=visit, args=(app,)) for app in apps]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(30)
    assert len(ids) == 4 and len(set(ids)) == 1
    assert _slow_jobs(db) == ids[:1]


def test_a_different_request_gets_a_job_of_its_own(slow):
    db = slow.config.database
    app = _site(slow.config)
    with app.test_request_context("/"), slowpages.allow_pending():
        with pytest.raises(slowpages.Pending) as one:
            apiclient.get_phenomena(db)
        with pytest.raises(slowpages.Pending) as other:
            apiclient.get_phenomena(db, types=["black_hole"])
    assert one.value.job["id"] != other.value.job["id"]


# --- The answer, the status, the invalidation ----------------------------------------------------------

def _wait_for(job_id, app, states=("succeeded", "failed", "cancelled"), seconds=180):
    client = app.test_client()
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        response = client.get(f"/api/page-jobs/{job_id}")
        assert response.status_code == 200, response.get_data(as_text=True)
        body = response.get_json()
        if body["state"] in states:
            return body
        time.sleep(0.3)
    raise AssertionError(f"job {job_id} did not finish: {body}")


def test_a_worker_builds_the_answer_into_the_cache_and_an_edit_invalidates_it(slow, monkeypatch):
    monkeypatch.setattr(workers, "ensure_workers", slow.run_jobs)         # the job runs as it is queued
    config, db = slow.config, slow.config.database
    real_ensure = pagejobs.ensure                                          # the site's own database is the default; here it is this test's
    monkeypatch.setattr(pagejobs, "ensure", lambda db, target, rows=None, expected=None, config=None, retry=False:
                        real_ensure(db, target, rows, expected, slow.config, retry))
    adminAuth.bootstrap_control_schema(config)
    app = _site(config)
    with pytest.raises(slowpages.Pending) as pending:
        _list(app, db)
    job = pending.value.job
    done = _wait_for(job["id"], app)
    assert done["state"] == "succeeded", done
    assert done["done"] == done["total"] and done["eta_seconds"] == 0.0 and done["estimated"] is True
    assert redisclient.client().exists(done["cache_key"])

    # The page now finds the answer where the job stored it, in a process that never saw it, with no API call.
    seen = []
    real = apiclient._send
    monkeypatch.setattr(apiclient, "_send", lambda method, target, *a, **k: (seen.append(target), real(method, target, *a, **k))[1])
    answer = _list(_site(config), db)
    assert answer["total"] >= 1
    assert not [target for target in seen if target.startswith("/phenomena")]

    # It is on the admin runs page, beside the warm-up's jobs.
    conn = store.get_connection(config)
    try:
        roots, _total = workQueue.list_roots(conn)
        [node] = [root for root in roots if root["kind"] == "slow-page"]
        assert node["title"].startswith("Slow page: /api/phenomena") and node["status"] == "done"
        tree = workQueue.load_tree(conn, node["id"])
        assert any(task["kind"] == "slow-page-line" and "bytes in" in task["task_key"] for task in tree["tasks"])
    finally:
        conn.close()

    # An edit that bumps what the answer read throws it away: the next request, still slow, queues a new job (this
    # build was quick on a tiny test galaxy, and the cost rule trusts a measured time, so say it was not).
    pagejobs.record_duration(db, f"/phenomena?db={db}", tuning.SLOW_PAGE_SECONDS * 3)
    dataversion.bump(db, dataversion.TAG_EPOCH)
    with pytest.raises(slowpages.Pending) as rebuilt:
        _list(_site(config), db)
    assert rebuilt.value.job["id"] != job["id"]


def test_a_job_that_fails_says_why_and_is_replaced_by_a_retry(slow, monkeypatch):
    monkeypatch.setattr(workers, "ensure_workers", slow.run_jobs)
    config, db = slow.config, slow.config.database
    app = _site(config)
    target = f"/phenomena?db={db}&sort=nonsense"
    job = pagejobs.ensure(db, target, 10, 1.0, config)
    failed = _wait_for(job["id"], app)
    assert failed["state"] == "failed" and failed["error"]
    assert pagejobs.ensure(db, target, 10, 1.0, config)["id"] == job["id"]          # shown, not hammered
    assert pagejobs.ensure(db, target, 10, 1.0, config, retry=True)["id"] != job["id"]


# --- The status route ----------------------------------------------------------------------------------

def test_the_status_route_reports_a_queued_job_and_refuses_what_is_not_one(slow):
    db = slow.config.database
    app = _site(slow.config)
    with pytest.raises(slowpages.Pending) as pending:
        _list(app, db)
    client = app.test_client()
    body = client.get(f"/api/page-jobs/{pending.value.job['id']}").get_json()
    assert body["state"] == "queued" and body["done"] == 0 and body["unit"] == "rows"
    assert body["target"] == f"/phenomena?db={db}" and body["cache_key"].startswith("planetgen:page:")
    assert client.get("/api/page-jobs/0123456789abcdef").status_code == 404
    assert client.get("/api/page-jobs/not-an-id").status_code == 404
    other = api_jobs.submit(api_jobs.run_command, ["true"], {}, "/", 5)       # an admin's job is not a page job
    try:
        assert client.get(f"/api/page-jobs/{other}").status_code == 404
    finally:
        api_jobs.cancel(other)


def test_the_progress_of_a_running_job_is_an_estimate_that_never_reaches_the_end(slow, monkeypatch):
    db = slow.config.database
    job = pagejobs.ensure(db, f"/phenomena?db={db}", 1000, 10.0, slow.config)
    connection = redisclient.client(decode_responses=True)
    now = time.time()
    monkeypatch.setattr(api_jobs, "status", lambda queue_job: {"state": "running", "result": None, "error": None})
    connection.hset(pagejobs.meta_key(job["id"]), "started", f"{now - 4:.3f}")
    running = pagejobs.status(job["id"])
    assert running["state"] == "running" and 350 <= running["done"] <= 450 and 5.5 < running["eta_seconds"] < 6.5
    connection.hset(pagejobs.meta_key(job["id"]), "started", f"{now - 100:.3f}")
    late = pagejobs.status(job["id"])
    assert late["done"] == 950 and late["eta_seconds"] is None            # past its expected time: still working


# --- Cancel --------------------------------------------------------------------------------------------

def test_a_cancel_before_a_worker_takes_the_job_removes_it(slow):
    db = slow.config.database
    target = f"/phenomena?db={db}"
    job = pagejobs.ensure(db, target, 10, 1.0, slow.config)
    assert pagejobs.cancel(job["id"])["state"] == "cancelled"
    assert pagejobs.ensure(db, target, 10, 1.0, slow.config)["state"] == "cancelled"        # shown for a moment
    assert pagejobs.ensure(db, target, 10, 1.0, slow.config, retry=True)["id"] != job["id"]
    assert pagejobs.cancel("0123456789abcdef") is None


def test_a_job_whose_cancel_arrives_before_its_query_builds_nothing(slow):
    db = slow.config.database
    target = f"/phenomena?db={db}"
    job = pagejobs.ensure(db, target, 10, 1.0, slow.config)
    redisclient.client(decode_responses=True).hset(pagejobs.meta_key(job["id"]), "cancelled", "1")
    assert pagejobs.build_answer(db, target, job["id"], slow.config) == {"cancelled": True}
    assert not redisclient.client().exists(pagecache.shared_key(db, target))
    assert slowpages.take_notice() is None


def test_only_an_admin_cancels_the_job_every_waiting_visitor_shares(slow):
    app = _site(slow.config)
    with pytest.raises(slowpages.Pending) as pending:
        _list(app, slow.config.database)
    response = app.test_client().post(f"/api/page-jobs/{pending.value.job['id']}/cancel")
    assert response.status_code in (401, 403)
    assert pagejobs.status(pending.value.job["id"])["state"] == "queued"


# --- No queue ------------------------------------------------------------------------------------------

def test_without_a_queue_the_answer_is_built_in_the_request_and_the_page_is_told(slow, monkeypatch):
    def refuse(function, *args):
        raise api_jobs.NoQueue("no Redis server answers")

    monkeypatch.setattr(api_jobs, "submit", refuse)
    db = slow.config.database
    app = _site(slow.config)
    slowpages.take_notice()
    answer = _list(app, db)
    assert answer["total"] >= 1
    assert "work queue is not available" in slowpages.take_notice()
    assert slowpages.take_notice() is None
    assert _slow_jobs(db) == []                                             # nothing left half queued


def test_an_update_in_progress_builds_in_the_request_too(slow, monkeypatch):
    from planetgen.admin import updateguard
    monkeypatch.setattr(updateguard, "is_updating", lambda: True)
    monkeypatch.setattr(updateguard, "inside_update", lambda: False)
    app = _site(slow.config)
    assert _list(app, slow.config.database)["total"] >= 1
    assert "work queue is not available" in slowpages.take_notice()
