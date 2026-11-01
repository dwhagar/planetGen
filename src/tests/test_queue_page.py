# tests/test_queue_page.py

"""
The admin Queue page and its API read and control the runs in Redis (ADM.10, PERF.99): the list, one run with its stages
and log, cancel, the stale-run Clear, delete, retry, the log download, the whole-queue pause and the worker count. The
controls themselves (`runcontrol`) are tested in `tests/test_run_control.py`.

Tests that take the `mysql_config` fixture (see `conftest.py`) are skipped, not failed, when no MySQL test server is
configured/reachable; those that take `redis_server` when there is no Redis.
"""

import pytest
from html import escape as html_escape

from planetgen.admin import auth as adminAuth
from planetgen.api.config import Config
from planetgen.db import store
from planetgen.queue import commands, progress_rate, redisqueue, runcontrol, runrecord, units
from planetgen.queue.runrecord import Stage
from planetgen.web import csrf, jobs, queue_page
from planetgen.web.app import create_app
from planetgen.web.lib import apiclient

KIND = "page-test"


@pytest.fixture
def control_config(mysql_config, monkeypatch):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    return mysql_config


@pytest.fixture
def web_app(control_config, redis_server, monkeypatch, tmp_path):
    class RealConfig(Config):
        MYSQL_CONFIG = control_config
        WRITE_MYSQL_CONFIG = control_config
        CONTROL_MYSQL_CONFIG = control_config
        WEB_DATABASE = control_config.database
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))
    monkeypatch.setitem(commands.REGISTRY, KIND, commands.CommandSpec(KIND, "A test", unit=(__name__, "_noop"),
                                                                     expand=(__name__, "_noop")))

    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")
    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    application = create_app(RealConfig)
    application.testing = True
    return application


@pytest.fixture
def admin_client(web_app, control_config):
    client = web_app.test_client()
    _username, password = adminAuth.bootstrap_control_schema(control_config)
    assert client.post("/api/auth/login", json={"username": "admin", "password": password}).status_code == 200
    assert client.post("/api/auth/change-credentials", json={
        "current_password": password, "new_username": "boss", "new_password": "a-long-new-password-1",
    }).status_code == 200
    return client


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


def _post(app, client, **form):
    nonce = "n" * 43
    client.set_cookie(csrf.COOKIE_NAME, nonce)
    session = client.get_cookie("pg_admin_session")
    with app.app_context():
        form[csrf.FIELD_NAME] = csrf._sign(nonce, session.value if session else "")
    return client.post("/admin/queue/action", data=form)


def _audit(config):
    conn = store.get_control_connection(config)
    try:
        return [row["action"] for row in conn.execute("SELECT action FROM admin_audit_log ORDER BY id").fetchall()]
    finally:
        conn.close()


def _run(conn, title="Fill shell 1", units_count=4, done=1, argv=("galaxy", "--shell", "1"), database=None, queued=True,
         **kwargs):
    """A live run with a started first stage: `units_count` units, `done` of them finished, the rest waiting."""
    run_id = runrecord.start_run(conn, title, [Stage("sectors", "Sectors", units=units_count), Stage("settle", "Settle")],
                                 who="boss", argv=list(argv), database=database, **kwargs)
    runrecord.begin_stage(conn, run_id, 0)
    if queued:
        units.enqueue_units(conn, run_id, 0, KIND, [({"n": n}, 1.0) for n in range(units_count)])
    for unit in range(done):
        runrecord.unit_done(conn, run_id, 0, unit, seconds=2.0)
    runrecord.log(conn, run_id, "sector 1,0,0 saved")
    return run_id


def _noop(*args):
    return {}


def test_the_api_needs_an_admin(web_app):
    client = web_app.test_client()
    assert client.get("/api/admin/work").status_code == 401
    assert client.get("/api/admin/work/20260101-000000-deadbeef").status_code == 401
    assert client.post("/api/admin/work/queue", json={"action": "pause"}).status_code == 401
    assert client.post("/api/admin/work/20260101-000000-deadbeef/control", json={"action": "cancel"}).status_code == 401


def test_the_api_lists_the_runs_shows_one_and_cancels_it(admin_client, control_config, conn):
    run_id = _run(conn, queued=False)
    body = admin_client.get("/api/admin/work").get_json()
    assert body["available"] and body["total"] == 1
    item = body["items"][0]
    assert item["id"] == run_id and item["live"] and item["state"] == "running" and item["units"] == 5
    assert item["units_done"] == 1 and item["stage"]["label"] == "Sectors" and item["eta"] is not None
    assert body["status"]["runs_active"] == 1 and body["status"]["holder"] == run_id and body["status"]["load"]["text"]
    detail = admin_client.get(f"/api/admin/work/{run_id}").get_json()
    assert [stage["label"] for stage in detail["run"]["stages"]] == ["Sectors", "Settle"]
    assert detail["run"]["stages"][0]["state"] == "running" and detail["run"]["stages"][1]["state"] == "waiting"
    assert any(event["message"] == "sector 1,0,0 saved" for event in detail["events"])
    response = admin_client.post(f"/api/admin/work/{run_id}/control", json={"action": "cancel"})
    assert response.get_json()["ok"] is True
    assert runrecord.read_run(conn, run_id)["state"] == "cancelled"
    assert admin_client.post(f"/api/admin/work/{run_id}/control", json={"action": "cancel"}).get_json()["ok"] is False
    assert admin_client.post(f"/api/admin/work/{run_id}/control", json={"action": "bogus"}).status_code == 400
    assert admin_client.post(f"/api/admin/work/{run_id}/control", json={"action": "pause"}).status_code == 400
    assert admin_client.get("/api/admin/work/not-an-id").status_code == 404
    assert admin_client.get("/api/admin/work/20260101-000000-00000000").status_code == 404
    assert admin_client.post("/api/admin/work/20260101-000000-00000000/control", json={"action": "cancel"}).status_code == 404
    assert admin_client.post("/api/admin/work/queue", json={"action": "pause"}).get_json() == {"ok": True}
    assert admin_client.post("/api/admin/work/queue", json={"action": "resume"}).get_json() == {"ok": True}
    assert admin_client.post(f"/api/admin/work/{run_id}/delete").get_json() == {"ok": True}
    assert _audit(control_config)[-4:] == ["work.cancel", "work.queue.pause", "work.queue.resume", "work.delete"]


def test_a_live_run_cannot_be_deleted_and_clear_finished_keeps_it(admin_client, conn):
    live = _run(conn, title="Still going", queued=False)
    assert admin_client.post(f"/api/admin/work/{live}/delete").get_json() == {"ok": False}
    assert admin_client.post("/api/admin/work/clear-finished").get_json() == {"cleared": 0}
    assert admin_client.get("/api/admin/work").get_json()["total"] == 1


def test_the_queue_page_shows_workers_load_and_runs(admin_client, conn):
    run_id = _run(conn, queued=False)
    html = admin_client.get("/admin/queue").get_data(as_text=True)
    assert "Workers active" in html and "Load average" in html and "Pause the queue" in html
    assert "Fill shell 1" in html and f"/admin/queue/{run_id}" in html
    assert "The run that owns the machine is" in html


def test_the_page_says_when_redis_is_not_answering(admin_client, monkeypatch):
    from planetgen.api import workqueue
    monkeypatch.setattr(workqueue, "_connection", lambda: None)
    body = admin_client.get("/api/admin/work").get_json()
    assert body["available"] is False and body["items"] == []
    assert "Redis isn't answering" in admin_client.get("/admin/queue").get_data(as_text=True)


def test_the_runs_table_route_lists_runs_with_their_controls(web_app, admin_client, conn):
    assert web_app.test_client().get("/table/queue-jobs").status_code == 403
    run_id = _run(conn, queued=False)
    data = admin_client.get("/table/queue-jobs").get_json()
    assert data["total"] == 1
    title, status, started, _duration, left, progress, actions = data["rows"][0]
    assert title["text"] == "Fill shell 1" and title["href"] == f"/admin/queue/{run_id}"
    assert status["text"] == "Running" and "stage 1 of 2: Sectors" in progress["text"] and "1 of 5 units" in progress["text"]
    labels = {part["text"]: part["href"] for part in actions["parts"] if isinstance(part, dict)}
    assert labels["Cancel"] == f"/admin/queue/confirm/cancel?run={run_id}"
    assert "Clear" not in labels, "a run that answers is cancelled, not cleared"


def test_clear_finished_deletes_every_finished_run_and_keeps_the_live_one(web_app, admin_client, control_config, conn):
    # ADM.41: one action on the Queue page clears the past-runs list; a running run is never touched.
    for number in (1, 2):
        old = _run(conn, title=f"Old run {number}", queued=False)
        runrecord.finish_run(conn, old, "done")
    _run(conn, title="Still going", queued=False)
    page = admin_client.get("/admin/queue/confirm/clear-finished")
    assert page.status_code == 200 and "Clear finished runs" in page.get_data(as_text=True)
    assert _post(web_app, admin_client, action="clear-finished", run="").status_code == 303
    body = admin_client.get("/api/admin/work").get_json()
    assert body["total"] == 1 and body["items"][0]["title"] == "Still going"
    assert "work.clear" in _audit(control_config)


def test_the_time_left_and_the_progress_read_as_words():
    estimate = {"seconds": 125.0, "stalled": False}
    low, high = progress_rate.eta_range(125.0)
    assert queue_page._eta_text(estimate, True) == f"{queue_page._duration(low)} to {queue_page._duration(high)}"
    assert queue_page._eta_text(estimate, False) == ""
    assert queue_page._eta_text({"seconds": None, "stalled": False}, True) == "estimating"
    assert queue_page._eta_text(None, True) == "estimating"
    assert queue_page._eta_text({"seconds": 90.0, "stalled": True}, True) == "waiting on progress"
    run = {"stage": {"index": 1, "count": 3, "label": "Settle"}, "live": True, "units": 10, "units_done": 4}
    assert queue_page._progress_text(run) == "stage 2 of 3: Settle, 4 of 10 units"
    assert queue_page._progress_text({**run, "live": False}) == "3 stages, 4 of 10 units"


def test_the_run_page_shows_its_stages_and_log(admin_client, conn):
    run_id = _run(conn, queued=False)
    html = admin_client.get(f"/admin/queue/{run_id}").get_data(as_text=True)
    assert 'id="stage-1"' in html and 'id="stage-2"' in html and "1 of 4" in html
    assert "sector 1,0,0 saved" in html and "Download the log" in html
    assert f"/admin/queue/confirm/cancel?run={run_id}" in html
    assert admin_client.get("/admin/queue/20260101-000000-00000000").status_code == 404
    assert admin_client.get("/admin/queue/not-a-run").status_code == 404


def test_the_run_page_log_lines_are_the_display_lines_of_the_records_the_download_holds(admin_client, conn):
    """UX.123: the page draws each event with the logging layer's display line, the download writes the same record as a file line."""
    from planetgen.util import logrecord
    run_id = _run(conn, queued=False)
    records = runrecord.log_records(conn, run_id)
    html = admin_client.get(f"/admin/queue/{run_id}").get_data(as_text=True)
    assert records
    for record in records:
        assert html_escape(logrecord.display_line(record)) in html
    assert " INFO     run: sector 1,0,0 saved" in html
    downloaded = [logrecord.parse_file_line(line) for line in runrecord.log_text(conn, run_id).splitlines()]
    assert [logrecord.display_line(r) for r in downloaded] == [logrecord.display_line(r) for r in records]


def test_the_api_events_carry_the_fields_the_record_is_made_of(admin_client, conn):
    run_id = _run(conn, queued=False)
    events = admin_client.get(f"/api/admin/work/{run_id}").get_json()["events"]
    rebuilt = [runrecord.event_record(event) for event in events]
    assert [r.message for r in rebuilt] == [r.message for r in runrecord.log_records(conn, run_id)]
    assert any(r.fields.get("stage") for r in rebuilt)


def test_the_log_downloads_from_redis_and_then_from_the_file(admin_client, conn, tmp_path):
    run_id = _run(conn, queued=False)
    runcontrol.cancel_run(conn, run_id, who="boss")         # ends it and writes the file
    response = admin_client.get(f"/admin/queue/{run_id}/log")
    assert response.status_code == 200 and response.mimetype == "text/plain"
    assert f"run-{run_id}.log" in response.headers["Content-Disposition"]
    text = response.get_data(as_text=True)
    assert "sector 1,0,0 saved" in text and "cancelled by boss" in text
    runrecord.delete_run(conn, run_id)
    assert admin_client.get(f"/admin/queue/{run_id}").status_code == 404
    kept = admin_client.get(f"/admin/queue/{run_id}/log")
    assert kept.status_code == 200 and "sector 1,0,0 saved" in kept.get_data(as_text=True)
    assert admin_client.get("/admin/queue/20260101-000000-00000000/log").status_code == 404


def test_every_control_is_confirmed_first(web_app, admin_client, control_config, conn):
    run_id = _run(conn)
    page = admin_client.get(f"/admin/queue/confirm/cancel?run={run_id}").get_data(as_text=True)
    assert "Cancel “Fill shell 1”?" in page and 'name="action" value="cancel"' in page
    assert runrecord.read_run(conn, run_id)["state"] == "running"
    response = _post(web_app, admin_client, action="cancel", run=run_id)
    assert response.status_code == 303
    assert runrecord.read_run(conn, run_id)["state"] == "cancelled"
    assert runrecord.active_run(conn) is None
    assert "Cancelled" in admin_client.get(response.headers["Location"]).get_data(as_text=True)
    assert admin_client.get("/admin/queue/confirm/nonsense").status_code == 404
    assert _post(web_app, admin_client, action="queue-pause").status_code == 303
    assert "The queue is paused" in admin_client.get("/admin/queue").get_data(as_text=True)
    _post(web_app, admin_client, action="queue-resume")
    assert "The queue is paused" not in admin_client.get("/admin/queue").get_data(as_text=True)


def test_cancel_on_a_run_page_returns_to_the_page_it_was_opened_from(web_app, admin_client, control_config, conn):
    """UX.256: Cancel on a run's page goes back to the page it was opened from (same site only), else the run's page."""
    run_id = _run(conn)
    page = admin_client.get(f"/admin/queue/{run_id}", headers={"Referer": "http://localhost/admin/generate/jobs/x"}).get_data(as_text=True)
    assert f"/admin/queue/confirm/cancel?run={run_id}&amp;next=/admin/generate/jobs/x" in page.replace("%2F", "/")
    page = admin_client.get(f"/admin/queue/{run_id}", headers={"Referer": "http://localhost/admin/queue/confirm/cancel"}).get_data(as_text=True)
    assert "next=" not in page
    confirm = admin_client.get(f"/admin/queue/confirm/cancel?run={run_id}&next=/admin/queue").get_data(as_text=True)
    assert 'name="next" value="/admin/queue"' in confirm
    hostile = admin_client.get(f"/admin/queue/confirm/cancel?run={run_id}&next=//evil.example").get_data(as_text=True)
    assert 'name="next"' not in hostile
    response = _post(web_app, admin_client, action="cancel", run=run_id, next="/admin/queue")
    assert response.status_code == 303 and response.headers["Location"].endswith("/admin/queue")
    other = _run(conn, title="Another")
    response = _post(web_app, admin_client, action="cancel", run=other, next="https://evil.example/")
    assert response.headers["Location"].endswith(f"/admin/queue/{other}")


def test_a_run_that_stopped_answering_is_offered_a_clear(web_app, admin_client, conn):
    run_id = _run(conn, queued=False)
    conn.hset(runrecord.run_key(run_id), "updated", 1.0)                  # not heard from for a long time
    assert admin_client.get("/api/admin/work").get_json()["status"]["stale"] is True
    html = admin_client.get("/admin/queue").get_data(as_text=True)
    assert "Clear the stale run" in html and f"/admin/queue/confirm/clear?run={run_id}" in html
    assert "Clear" in admin_client.get(f"/admin/queue/confirm/clear?run={run_id}").get_data(as_text=True)
    response = _post(web_app, admin_client, action="clear", run=run_id)
    assert response.status_code == 303
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "failed" and run["error"].startswith("cleared by boss")
    assert runrecord.active_run(conn) is None


def test_a_run_that_answers_is_not_cleared(web_app, admin_client, conn):
    run_id = _run(conn, queued=False)
    runrecord.touch(conn, run_id)
    response = _post(web_app, admin_client, action="clear", run=run_id)
    assert "still answering" in admin_client.get(response.headers["Location"]).get_data(as_text=True)
    assert runrecord.read_run(conn, run_id)["state"] == "running"


def test_the_action_needs_csrf(admin_client):
    assert admin_client.post("/admin/queue/action", data={"action": "queue-pause"}).status_code == 400


def test_retry_reruns_the_command_a_run_was_started_by(web_app, admin_client, control_config, conn, monkeypatch):
    started = []
    monkeypatch.setattr(jobs, "start_job", lambda kind, title, steps, **kw: started.append((kind, title, steps, kw))
                        or "20261001-000000-abcd")
    run_id = _run(conn, title="planetgen galaxy --shell 2", argv=("galaxy", "--shell", "2"),
                  database=control_config.database, queued=False)
    assert "/admin/queue/confirm/retry" not in admin_client.get(f"/admin/queue/{run_id}").get_data(as_text=True), \
        "a live run is cancelled, not retried"
    runrecord.finish_run(conn, run_id, "failed", error="unit 3 failed")
    page = admin_client.get(f"/admin/queue/{run_id}").get_data(as_text=True)
    assert f"/admin/queue/confirm/retry?run={run_id}" in page and "unit 3 failed" in page
    _post(web_app, admin_client, action="retry", run=run_id)
    kind, title, steps, kw = started[-1]
    assert kind == "galaxy" and steps[0]["argv"][-3:] == ["galaxy", "--shell", "2"]
    assert kw["database"] == control_config.database


def test_the_queue_page_offers_the_two_fills_and_each_confirms_before_it_starts(admin_client):
    """ADM.129: the buttons are on the Queue page, one place, each with what it does."""
    page = admin_client.get("/admin/queue").get_data(as_text=True)
    assert "/admin/queue/confirm/fill-habitability" in page and "/admin/queue/confirm/fill-map-luminosity" in page
    assert "does not re-draw a planet" in page
    confirm = admin_client.get("/admin/queue/confirm/fill-habitability").get_data(as_text=True)
    assert "Regenerate habitability?" in confirm and "does not re-draw a planet" in confirm
    assert "Fill map luminosity?" in admin_client.get("/admin/queue/confirm/fill-map-luminosity").get_data(as_text=True)


@pytest.mark.parametrize("action, command, title", [
    ("fill-habitability", "rescore-habitability", "Regenerate habitability"),
    ("fill-map-luminosity", "fill-map-luminosity", "Fill map luminosity")])
def test_a_fill_button_starts_its_command_as_a_job_and_writes_the_audit_entry(
        web_app, admin_client, control_config, monkeypatch, action, command, title):
    started, events = [], []
    monkeypatch.setattr(jobs, "start_job", lambda kind, title, steps, **kw: started.append((kind, title, steps, kw))
                        or "20261001-000000-abcd")
    monkeypatch.setattr(queue_page.activity_log, "event", lambda *args, **kw: events.append((args, kw)))
    response = _post(web_app, admin_client, action=action)
    assert response.status_code == 303 and "/generate" in response.headers["Location"]
    kind, job_title, steps, kw = started[-1]
    assert kind == command and job_title == title and steps[0]["argv"][-1] == command
    assert kw["database"] == control_config.database and kw["admin"]
    assert events and events[-1][0][:2] == ("GEN", "job.start") and events[-1][1]["kind"] == command
    assert events[-1][1]["started_from"] == "queue page"


def test_a_fill_behind_a_live_run_or_a_running_job_waits_its_turn(web_app, admin_client, conn, monkeypatch):
    """PERF.154: nothing is refused for a busy machine; a fill is a job and waits in the list for the job before it
    (its own run then waits for the machine)."""
    started = []
    monkeypatch.setattr(jobs, "start_job", lambda *a, **k: started.append(a[1]) or "20261001-000000-aaaa")
    live = _run(conn, title="Fill shell 7", queued=False)
    runrecord.touch(conn, live)
    response = _post(web_app, admin_client, action="fill-habitability")
    assert response.status_code == 303 and started, "a live run does not refuse a fill"
    runcontrol.cancel_run(conn, live)

    def busy(*args, **kwargs):
        raise jobs.JobBusy({"id": "20261001-000000-abcd", "title": "Generate sectors"})
    monkeypatch.setattr(jobs, "start_job", busy)
    monkeypatch.setattr(jobs, "promote_queued", lambda *a, **k: None)
    response = _post(web_app, admin_client, action="fill-map-luminosity")
    assert response.status_code == 303
    assert [job["title"] for job in jobs.queued_jobs()] == ["Fill map luminosity"]


def test_retry_refuses_another_database_and_a_run_with_no_command(web_app, admin_client, conn, monkeypatch):
    monkeypatch.setattr(jobs, "start_job", lambda *a, **k: pytest.fail("should not start"))
    elsewhere = _run(conn, title="Elsewhere", argv=("galaxy", "--shell", "2"), database="some_other_galaxy", queued=False)
    runrecord.finish_run(conn, elsewhere, "failed")
    response = _post(web_app, admin_client, action="retry", run=elsewhere)
    assert "some_other_galaxy" in admin_client.get(response.headers["Location"]).get_data(as_text=True)
    bare = _run(conn, title="No command", argv=(), queued=False)
    runrecord.finish_run(conn, bare, "failed")
    assert "/admin/queue/confirm/retry" not in admin_client.get(f"/admin/queue/{bare}").get_data(as_text=True)
    maintenance = _run(conn, title="Warm", kind=runrecord.MAINTENANCE, argv=("warm",), queued=False)
    runrecord.finish_run(conn, maintenance, "failed")
    assert "/admin/queue/confirm/retry" not in admin_client.get(f"/admin/queue/{maintenance}").get_data(as_text=True)


def test_retry_plan_for_a_run_a_web_job_started_begins_at_the_step_that_failed(tmp_path, monkeypatch):
    root_dir = tmp_path / "jobs"
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(root_dir))
    job_id = jobs.start_job("new_galaxy", "New galaxy", [
        {"label": "Reset", "argv": ["x"]}, {"label": "Plan", "argv": ["y"]}, {"label": "Fill", "argv": ["z"]},
    ], spawn=False)
    path = root_dir / job_id
    (path / "state.json").write_text('{"status": "failed", "step": 2, "started_at": 1, "finished_at": 2}')
    jobs.release_lock(job_id)
    job, steps = jobs.remaining_steps(job_id)
    assert job["status"] == "failed" and [step["label"] for step in steps] == ["Plan", "Fill"]
    run = {"kind": "generation", "job": job_id, "title": "Fill", "argv": ["galaxy"], "database": None}
    plan = queue_page._retry_plan(run)
    assert plan["describe"] == "Plan, then Fill" and plan["title"] == "Retry: New galaxy"
    (path / "state.json").write_text('{"status": "succeeded", "step": 3, "started_at": 1, "finished_at": 2}')
    assert jobs.remaining_steps(job_id) is None
    assert queue_page._retry_plan(run)["describe"] == "planetgen galaxy"        # the job is done: the command alone


def test_the_worker_count_is_read_and_set_through_the_api(admin_client, redis_server):
    """ADM.15: GET shows the bulk worker count; POST sets a count or "Ludicrous Speed", in the audit log."""
    from planetgen.queue import workers
    connection = redisqueue.connect()
    connection.delete(workers.WANTED_KEY)
    try:
        state = admin_client.get("/api/admin/work/workers").get_json()
        assert set(state) == {"wanted", "live", "ludicrous", "cpus", "host_cap", "host_running", "host_refused"}
        assert state["wanted"] >= 1
        done = admin_client.post("/api/admin/work/workers", json={"count": 3})
        assert done.status_code == 200 and done.get_json()["wanted"] == 3 and workers.wanted(connection) == 3
        done = admin_client.post("/api/admin/work/workers", json={"ludicrous": True})
        assert done.get_json()["wanted"] == workers.ludicrous_count()
        for bad in ({}, {"count": 0}, {"count": "4"}, {"count": 2, "ludicrous": True}, {"ludicrous": False}):
            assert admin_client.post("/api/admin/work/workers", json=bad).status_code == 400, bad
    finally:
        connection.delete(workers.WANTED_KEY)


def test_the_worker_count_needs_an_admin(web_app):
    client = web_app.test_client()
    assert client.get("/api/admin/work/workers").status_code == 401
    assert client.post("/api/admin/work/workers", json={"count": 2}).status_code == 401


def test_the_queue_page_sets_the_worker_count(web_app, admin_client, redis_server):
    from planetgen.queue import workers
    connection = redisqueue.connect()
    connection.delete(workers.WANTED_KEY)
    try:
        page = admin_client.get("/admin/queue").get_data(as_text=True)
        assert 'name="count"' in page and "Ludicrous Speed" in page
        response = _post(web_app, admin_client, action="workers-set", count="2")
        assert response.status_code == 303 and workers.wanted(connection) == 2
        assert "now asks for 2 workers" in admin_client.get(response.headers["Location"]).get_data(as_text=True)
        _post(web_app, admin_client, action="workers-set", count="0")
        assert workers.wanted(connection) == 2                                  # refused: still 2
        _post(web_app, admin_client, action="workers-ludicrous")
        assert workers.wanted(connection) == workers.ludicrous_count()
    finally:
        connection.delete(workers.WANTED_KEY)


def test_stop_everything_is_confirmed_holds_everything_and_resume_lifts_it(web_app, admin_client, conn, monkeypatch):
    """PERF.186: the Queue page's Stop everything runs the shared routine; the banner shows on every admin page; Resume lifts it."""
    from planetgen.queue import pause, stopall

    monkeypatch.setattr(stopall, "kill_processes", lambda *a, **k: {"killed": 5, "survivors": []})
    run_id = _run(conn)
    html = admin_client.get("/admin/queue").get_data(as_text=True)
    assert "Stop everything" in html and "Everything is stopped" not in html
    confirm = admin_client.get("/admin/queue/confirm/queue-stop-all").get_data(as_text=True)
    assert "Stop everything?" in confirm and "ends every planetgen worker" in confirm
    try:
        response = _post(web_app, admin_client, action="queue-stop-all")
        assert response.status_code == 303
        assert runrecord.read_run(conn, run_id)["state"] == "cancelled"
        assert stopall.state()["stopped"] and stopall.state()["by"] == "boss" and pause.is_paused(conn)
        page = admin_client.get("/admin/queue").get_data(as_text=True)
        assert "Everything is stopped" in page and "5 process(es) ended" in page
        assert 'href="/admin/queue/confirm/queue-resume-all"' in page
        assert "Everything is stopped" in admin_client.get("/admin").get_data(as_text=True)     # the banner, on any admin page
        assert admin_client.get("/api/admin/work").get_json()["stopped"]["stopped"] is True
        assert _post(web_app, admin_client, action="queue-resume-all").status_code == 303
        assert not stopall.state()["stopped"] and not pause.is_paused(conn)
        assert "Everything is stopped" not in admin_client.get("/admin/queue").get_data(as_text=True)
    finally:
        stopall.resume_all(conn)


def test_stop_everything_works_when_redis_is_not_answering(web_app, admin_client, monkeypatch):
    from planetgen.api import workqueue
    from planetgen.queue import stopall

    monkeypatch.setattr(workqueue, "_connection", lambda: None)
    monkeypatch.setattr(stopall, "kill_processes", lambda *a, **k: {"killed": 0, "survivors": []})
    try:
        assert _post(web_app, admin_client, action="queue-stop-all").status_code == 303
        page = admin_client.get("/admin/queue").get_data(as_text=True)
        assert "Everything is stopped" in page and "Redis isn't answering" in page
    finally:
        stopall.resume_all()


def test_no_queue_tile_is_a_bare_number(admin_client, conn):
    """UX.212: a tile's value names what it counts."""
    import re

    _run(conn, queued=False)
    html = admin_client.get("/admin/queue").get_data(as_text=True)
    values = re.findall(r'<span class="stat-value">([^<]*)</span>', html)
    assert len(values) == 3
    assert [v for v in values if re.fullmatch(r"[\d.,\s/]*", v.strip()) and "/" not in v] == []
    assert any(v.endswith(("worker", "workers")) for v in values) and any(v.endswith(("run", "runs")) for v in values)
