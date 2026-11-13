# tests/test_queue_page.py

"""
The admin Queue page and its API read and control the runs in Redis (ADM.10, PERF.99): the list, one run with its stages
and log, cancel, the stale-run Clear, delete, retry, the log download, the whole-queue pause and the worker count. The
controls themselves (`runcontrol`) are tested in `tests/test_run_control.py`.

Tests that take the `mysql_config` fixture (see `conftest.py`) are skipped, not failed, when no MySQL test server is
configured/reachable; those that take `redis_server` when there is no Redis.
"""

import ast
import inspect

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
    assert runrecord.read_run(conn, run_id)["state"] == "canceled"
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


@pytest.fixture
def held_pair(conn):
    """A queued job and a job held for it, on the check queue where no worker runs them (PERF.248)."""
    from planetgen.queue import requires
    queue = redisqueue.queue(redisqueue.CHECK_QUEUE, conn)
    queue.empty()
    requires.drop_held(conn, redisqueue.CHECK_QUEUE)
    first = requires.enqueue(queue, "tests.requires_jobs.mark", "first", job_id="held-test-first", job_timeout=60)
    second = requires.enqueue(queue, "tests.requires_jobs.mark", "second", requires=[first], job_id="held-test-second", job_timeout=60)
    yield first, second
    requires.drop_held(conn, redisqueue.CHECK_QUEUE)
    queue.empty()
    for job in (first, second):
        conn.delete(f"rq:job:{job.id}", f"rq:job:{job.id}:dependents")


def test_the_api_lists_the_jobs_held_for_a_job_they_require_and_not_the_queued_one(admin_client, held_pair):
    first, second = held_pair
    body = admin_client.get("/api/admin/work/held").get_json()
    assert body["total"] == 1
    item = body["items"][0]
    assert item["id"] == second.id and item["waits_for"] == [first.id] and item["queue"] == redisqueue.CHECK_QUEUE
    assert item["label"].endswith("mark") and item["run"] is None and item["stage"] is None


def test_the_held_jobs_need_an_admin(web_app):
    assert web_app.test_client().get("/api/admin/work/held").status_code == 401


def test_the_queue_page_has_a_held_jobs_panel_only_when_a_job_is_held(admin_client, conn):
    from planetgen.queue import requires
    requires.drop_held(conn, redisqueue.CHECK_QUEUE)
    assert "Held jobs" not in admin_client.get("/admin/queue").get_data(as_text=True)


def test_the_queue_page_names_the_held_job_and_what_it_waits_for(admin_client, held_pair):
    first, second = held_pair
    html = admin_client.get("/admin/queue").get_data(as_text=True)
    assert "Held jobs" in html and "1 job is held" in html and second.id in html and first.id in html and "check" in html


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
    assert status["text"] == "Running" and "Stage 1 of 2: Sectors" in progress["text"] and "1 of 5 units" in progress["text"]
    labels = {part["text"]: part["href"] for part in actions["parts"] if isinstance(part, dict)}
    assert labels["Cancel"] == f"/admin/queue/confirm/cancel?run={run_id}"
    assert "Clear" not in labels, "a run that answers is canceled, not cleared"


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
    assert queue_page._progress_text(run) == "Stage 2 of 3: Settle, 4 of 10 units"
    assert queue_page._progress_text({**run, "live": False}) == "3 stages, 4 of 10 units"


def test_the_run_page_shows_its_stages_and_log(admin_client, conn):
    run_id = _run(conn, queued=False)
    html = admin_client.get(f"/admin/queue/{run_id}").get_data(as_text=True)
    assert 'id="stage-1"' in html and 'id="stage-2"' in html and "1 of 4" in html
    assert "sector 1,0,0 saved" in html and "Download the log" in html
    assert f"/admin/queue/confirm/cancel?run={run_id}" in html
    assert admin_client.get("/admin/queue/20260101-000000-00000000").status_code == 404
    assert admin_client.get("/admin/queue/not-a-run").status_code == 404


def _publish_unit(conn, run_id, unit, **fields):
    """A verdict on a running unit as the follower publishes it (`unitwatch.view`)."""
    import json
    shown = {"state": "working", "severity": "", "reason": "", "actionable": False, "agreed": False, "counter": True, "rate": None, "projected_finish_seconds": None,
             "still_seconds": 0.0, "over_limit": False, "phase": "", "waiting_on": "", "done": 0.0, "estimate": 0.0, "restarts": 0, "resident_mb": None,
             "beat_age_seconds": 1.0, "running_seconds": 10.0, **fields}
    conn.hset(runrecord.unitstate_key(run_id, 0), unit, json.dumps(shown))


def test_the_run_page_lists_the_running_units_with_their_beat_counter_phase_and_state_the_worst_first(admin_client, conn):
    """OPS.139: the Queue run page draws the verdicts `GET /api/admin/work/<id>/units` serves, the units that need a look on top."""
    run_id = _run(conn, queued=False)
    _publish_unit(conn, run_id, 0, state="working", phase="drawing tier 3", done=1200.0, estimate=4000.0, rate=12.5, resident_mb=812.0, beat_age_seconds=2.0)
    _publish_unit(conn, run_id, 1, state="quiet", actionable=True, agreed=True, reason="no beat for 400 s", phase="waiting", waiting_on="gate", counter=False, beat_age_seconds=400.0)
    _publish_unit(conn, run_id, 2, state="slow", severity="flag", done=10.0, estimate=1000.0, rate=0.02, projected_finish_seconds=3600.0)
    html = admin_client.get(f"/admin/queue/{run_id}").get_data(as_text=True)
    assert "Running units" in html and "3 units: 1 working, 1 slow, 1 quiet" in html
    assert html.index("unit-quiet") < html.index("unit-slow") < html.index("unit-working")
    assert "1,200 of 4,000" in html and "drawing tier 3" in html and "12.5 a second" in html and "812 MB" in html
    assert "waiting (waiting on gate)" in html and "no counter" in html and "6 m 40 s ago" in html and "to be stopped" in html
    assert "50.0 s each" in html and "1 h 00 m" in html


def test_the_run_page_has_no_unit_table_without_running_units_or_once_the_run_is_over(admin_client, conn):
    run_id = _run(conn, queued=False)
    assert "Running units" not in admin_client.get(f"/admin/queue/{run_id}").get_data(as_text=True)
    _publish_unit(conn, run_id, 0)
    assert "Running units" in admin_client.get(f"/admin/queue/{run_id}").get_data(as_text=True)
    runcontrol.cancel_run(conn, run_id, who="boss")
    assert "Running units" not in admin_client.get(f"/admin/queue/{run_id}").get_data(as_text=True)


def test_the_unit_rows_are_capped_and_the_rest_counted(web_app):
    from planetgen.web import queue_page
    run = {"stages": [{"label": "Sectors"}]}
    shown = [{"stage": 0, "unit": number, "state": "working", "counter": True, "beat_age_seconds": 1.0} for number in range(queue_page.UNIT_ROWS_SHOWN + 3)]
    rows, summary, hidden = queue_page._unit_rows(run, shown)
    assert len(rows) == queue_page.UNIT_ROWS_SHOWN and hidden == 3 and summary.startswith(f"{len(shown)} units: ")
    assert queue_page._unit_rows(run, []) == ([], "", 0)


def test_the_run_page_log_lines_are_the_display_lines_of_the_records_the_download_holds(admin_client, conn):
    """UX.123: the page draws each event with the logging layer's display line, the download writes the same record as a file line."""
    from planetgen.util import logrecord
    run_id = _run(conn, queued=False)
    records = runrecord.log_records(conn, run_id)
    html = admin_client.get(f"/admin/queue/{run_id}").get_data(as_text=True)
    assert records
    for record in records:
        assert html_escape(logrecord.display_line(record)) in html
    assert " INFO     queue: sector 1,0,0 saved" in html
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
    assert "sector 1,0,0 saved" in text and "canceled by boss" in text
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
    assert runrecord.read_run(conn, run_id)["state"] == "canceled"
    assert runrecord.active_run(conn) is None
    assert "Canceled" in admin_client.get(response.headers["Location"]).get_data(as_text=True)
    assert admin_client.get("/admin/queue/confirm/nonsense").status_code == 404
    assert _post(web_app, admin_client, action="queue-pause").status_code == 303
    assert "The queue is paused" in admin_client.get("/admin/queue").get_data(as_text=True)
    _post(web_app, admin_client, action="queue-resume")
    assert "The queue is paused" not in admin_client.get("/admin/queue").get_data(as_text=True)


def test_a_solo_hold_is_named_on_the_page_with_who_holds_it_and_can_be_lifted_there(web_app, admin_client, control_config, conn, monkeypatch):
    """OPS.145: the pause that holds the queue for a restarted unit says why and by whom, and the button lifts it."""
    from planetgen.queue import pause
    monkeypatch.setattr(pause, "_solo_is_over", lambda connection, solo: False)
    assert pause.hold_solo("20261109-120000-abcd", 3, 7, "job-7", conn)
    page = admin_client.get("/admin/queue").get_data(as_text=True)
    assert "The queue is held by run 20261109-120000-abcd" in page and "solo run of unit 7 (stage 3) after a restart" in page and "Lift the hold" in page
    assert "The queue is paused" not in page
    assert _post(web_app, admin_client, action="queue-resume").status_code == 303
    assert not pause.is_paused(conn) and "The queue is held" not in admin_client.get("/admin/queue").get_data(as_text=True)


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
        "a live run is canceled, not retried"
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
        assert runrecord.read_run(conn, run_id)["state"] == "canceled"
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


def test_the_units_route_lists_what_the_follower_last_judged(admin_client, control_config, conn):
    """OPS.138: `GET /api/admin/work/<id>/units`."""
    import json
    from planetgen.queue import unitbeat, unitwatch
    run_id = _run(conn, queued=False)
    now = unitbeat.server_time(conn)
    for unit, still in ((2, 2000.0), (0, 0.0)):
        conn.hset(runrecord.progress_key(run_id, 0), unit, json.dumps({"done": 5.0, "estimate": 100.0, "beat": now - 1.0, "started": now - 2500.0, "moved": now - still,
                                                                       "phase": "p", "waiting_on": "", "restarts": 0}))
    unitwatch.UnitWatch(conn, run_id, 0, "test", record_error=lambda *a, **k: None, stuck=False).check()
    body = admin_client.get(f"/api/admin/work/{run_id}/units").get_json()
    assert [(item["stage"], item["unit"], item["state"]) for item in body["units"]] == [(0, 0, "working"), (0, 2, "quiet")]
    assert admin_client.get("/api/admin/work/20260101-000000-00000000/units").status_code == 404


# --- OPS.144: the killed flag and Resume run ----------------------------------------------------------------

def test_a_run_ended_from_outside_its_work_is_flagged_killed_and_one_that_failed_by_its_own_code_is_not(conn):
    """A canceled run, a cleared one, a unit RQ ended at its time limit and a unit whose worker was lost are killed; a unit that raised an error of its own is a failure, not a kill."""
    class JobTimeoutException(Exception):
        pass

    class Job:
        func_name = "x"

        def __init__(self, run_id):
            self.args = (run_id, 0, 1)

    canceled = _run(conn, title="canceled", queued=False)
    runcontrol.cancel_run(conn, canceled, who="boss")
    cleared = _run(conn, title="cleared", queued=False)
    runrecord.finish_run(conn, cleared, "failed", error="cleared by boss: it showed no sign of life", killed=True)
    timed_out = _run(conn, title="timed out", queued=False)
    units.unit_failed(Job(timed_out), conn, JobTimeoutException, JobTimeoutException("Task exceeded maximum timeout value (1885 seconds)"), None)
    broke = _run(conn, title="broke", queued=False)
    units.unit_failed(Job(broke), conn, ValueError, ValueError("bad value"), None)
    flags = {title: runrecord.read_run(conn, run_id) for title, run_id in
             (("canceled", canceled), ("cleared", cleared), ("timed out", timed_out), ("broke", broke))}
    assert [flags[title]["killed"] for title in ("canceled", "cleared", "timed out", "broke")] == [True, True, True, False]
    assert flags["canceled"]["killed_reason"] == "canceled by boss" and flags["broke"]["killed_reason"] is None
    assert "1885" in flags["timed out"]["killed_reason"]


def test_mark_killed_keeps_the_first_reason_and_works_on_a_live_run(conn):
    run_id = _run(conn, queued=False)
    assert runrecord.read_run(conn, run_id)["killed"] is False
    runrecord.mark_killed(conn, run_id, "stuck: counter frozen")
    runrecord.mark_killed(conn, run_id, "something later")
    run = runrecord.read_run(conn, run_id)
    assert run["killed"] and run["killed_reason"] == "stuck: counter frozen" and run["state"] == "running"


def test_resume_run_is_offered_on_a_killed_run_that_can_be_carried_on_and_starts_finish_run(web_app, admin_client, control_config, conn, monkeypatch):
    started = []
    monkeypatch.setattr(jobs, "start_job", lambda kind, title, steps, **kw: started.append((kind, title, steps, kw)) or "20261001-000000-abcd")
    from planetgen.web import generate_page
    asked = []

    def finish_status(database, run_id=None):
        asked.append(run_id)
        return {"error": None, "run": run_id, "title": "Generate sectors", "stopped_in": "Scatter", "carried": ["A", "B"], "again": ["Scatter", "Index"]}

    monkeypatch.setattr(generate_page, "finish_status", finish_status)
    killed = _run(conn, title="Generate sectors", database=control_config.database, queued=False)
    runrecord.finish_run(conn, killed, "failed", error="unit 1 of stage Scatter failed: JobTimeoutException", killed=True)
    plain = _run(conn, title="Plain failure", database=control_config.database, queued=False)
    runrecord.finish_run(conn, plain, "failed", error="unit 1 failed: ValueError")
    page = admin_client.get(f"/admin/queue/{killed}").get_data(as_text=True)
    assert f"/admin/queue/confirm/resume?run={killed}" in page and "This run was killed" in page
    assert "/admin/queue/confirm/resume" not in admin_client.get(f"/admin/queue/{plain}").get_data(as_text=True)
    confirm = admin_client.get(f"/admin/queue/confirm/resume?run={killed}").get_data(as_text=True)
    assert "keeps 2 finished stages" in confirm and "Scatter" in confirm
    response = _post(web_app, admin_client, action="resume", run=killed)
    assert response.status_code == 303 and "/generate" in response.headers["Location"]
    kind, title, steps, kw = started[-1]
    assert kind == "finish_run" and steps[0]["argv"][-3:] == ["finish-run", "--run", killed] and kw["database"] == control_config.database
    assert asked and set(asked) == {killed}


def test_resume_run_is_refused_when_the_run_cannot_be_carried_on(web_app, admin_client, control_config, conn, monkeypatch):
    monkeypatch.setattr(jobs, "start_job", lambda *a, **k: pytest.fail("should not start"))
    from planetgen.web import generate_page
    monkeypatch.setattr(generate_page, "finish_status", lambda database, run_id=None: {"error": "The galaxy was planned again."})
    killed = _run(conn, title="Generate sectors", database=control_config.database, queued=False)
    runrecord.finish_run(conn, killed, "canceled", error="canceled by boss")
    assert "/admin/queue/confirm/resume" not in admin_client.get(f"/admin/queue/{killed}").get_data(as_text=True)
    response = admin_client.get(f"/admin/queue/confirm/resume?run={killed}")
    assert response.status_code == 303
    response = _post(web_app, admin_client, action="resume", run=killed)
    assert "planned again" in admin_client.get(response.headers["Location"]).get_data(as_text=True)


def test_the_jobs_route_and_the_page_list_the_jobs_that_are_running_by_their_beats(admin_client, control_config, conn):
    """OPS.142: `GET /api/admin/work/jobs` and the Jobs running now panel."""
    from planetgen.queue import jobbeat
    assert admin_client.get("/api/admin/work/jobs").get_json()["jobs"] == []
    with jobbeat.beating("web", "job-7", label="Fill the map", conn=conn) as beat:
        beat.phase("step 1 of 2: Fill")
        beat.write()
        body = admin_client.get("/api/admin/work/jobs").get_json()
        assert [(job["kind"], job["job"], job["label"], job["state"], job["phase"]) for job in body["jobs"]] == [("web", "job-7", "Fill the map", "working", "step 1 of 2: Fill")]
        html = admin_client.get("/admin/queue").get_data(as_text=True)
        assert "Jobs running now" in html and "Generate page job" in html and "Fill the map" in html and "step 1 of 2: Fill" in html
    assert "Jobs running now" not in admin_client.get("/admin/queue").get_data(as_text=True)


def test_the_two_duration_formats_stay_apart():
    """The run, stage, ETA and heartbeat figures use `_duration` ("3 m 05 s", "12.4 s"); only the other-job rows use `_beat_gap` ("4 min"). OPS.142 once defined a second `_duration` that replaced the first at import."""
    assert queue_page._duration(185) == "3 m 05 s" and queue_page._duration(12.4) == "12.4 s" and queue_page._duration(None) == ""
    assert queue_page._beat_gap(5) == "5 s" and queue_page._beat_gap(240) == "4 min" and queue_page._beat_gap(7800) == "2 h 10 min"
    names = [node.name for node in ast.parse(inspect.getsource(queue_page)).body if isinstance(node, ast.FunctionDef)]
    assert len(names) == len(set(names)), "a function is defined twice in web/queue_page.py"


# --- OPS.158: what the queue has learned per job kind -----------------------------------------------------------------

@pytest.fixture
def gate_keys(key_prefix, monkeypatch):
    from planetgen.queue import memorygate, profileclaim
    monkeypatch.setattr(memorygate, "PREFIX", f"{key_prefix}:gate")
    monkeypatch.setenv(profileclaim.ENV_VAR, "1")


def test_the_tuning_api_needs_an_admin(web_app):
    client = web_app.test_client()
    assert client.get("/api/admin/work/tuning").status_code == 401
    assert client.post("/api/admin/work/tuning", json={"action": "clear", "kind": "unit:x"}).status_code == 401


def test_the_tuning_api_pins_unpins_and_clears_a_kind_and_logs_each(admin_client, control_config, conn, gate_keys):
    from planetgen.queue import memorygate, profileclaim
    memorygate.learn(conn, "unit:link-neighbors", 900.0)
    body = admin_client.post("/api/admin/work/tuning", json={"action": "pin", "kind": "unit:link-neighbors", "mb": 1500}).get_json()
    (row,) = [r for r in body["kinds"] if r["kind"] == "unit:link-neighbors"]
    assert row["pinned_mb"] == 1500.0 and row["typed_mb"] > 0 and body["self_tuning"] is True
    assert body["changes"][0]["what"] == "pin" and body["changes"][0]["new"] == 1500.0
    assert profileclaim.claim(conn, "unit:link-neighbors", 700.0, token="r1", read_rows=lambda kind: []) == 1500.0
    assert memorygate.claim_mb(conn, "unit:link-neighbors", 700.0, known={}) == 1500.0
    assert admin_client.post("/api/admin/work/tuning", json={"action": "unpin", "kind": "unit:link-neighbors"}).status_code == 200
    assert admin_client.post("/api/admin/work/tuning", json={"action": "unpin", "kind": "unit:link-neighbors"}).status_code == 404
    assert admin_client.post("/api/admin/work/tuning", json={"action": "clear", "kind": "unit:link-neighbors"}).status_code == 200
    assert memorygate.learned(conn).get("unit:link-neighbors") is None
    assert [c["what"] for c in admin_client.get("/api/admin/work/tuning").get_json()["changes"]][:3] == ["clear", "unpin", "pin"]
    assert {"work.tuning.pin", "work.tuning.unpin", "work.tuning.clear"} <= set(_audit(control_config))


@pytest.mark.parametrize("body", [{"action": "pin", "kind": "unit:x"}, {"action": "pin", "kind": "unit:x", "mb": 10}, {"action": "stir", "kind": "unit:x"},
                                  {"action": "clear", "kind": "no spaces or colons"}, {"action": "clear"}])
def test_the_tuning_api_refuses_a_bad_ask(admin_client, gate_keys, body):
    assert admin_client.post("/api/admin/work/tuning", json=body).status_code == 400


def test_the_queue_page_shows_the_learned_figures_and_the_latest_changes(admin_client, conn, gate_keys):
    from planetgen.queue import profileclaim
    assert "What the queue has learned" not in admin_client.get("/admin/queue").get_data(as_text=True)
    profileclaim.pin(conn, "unit:link-neighbors", 1500.0, "boss")
    html = admin_client.get("/admin/queue").get_data(as_text=True)
    assert "What the queue has learned" in html and "unit:link-neighbors" in html and "1,500 MB" in html and "fixed by boss" in html


def test_the_queue_page_controls_pin_release_and_forget(web_app, admin_client, conn, gate_keys):
    from planetgen.queue import profileclaim
    assert _post(web_app, admin_client, action="tuning-pin", kind="unit:link-neighbors", mb="1800").status_code == 303
    assert profileclaim.pinned_mb(conn, "unit:link-neighbors") == 1800.0
    assert _post(web_app, admin_client, action="tuning-pin", kind="unit:link-neighbors", mb="").status_code == 303
    assert profileclaim.pinned_mb(conn, "unit:link-neighbors") == 1800.0
    _post(web_app, admin_client, action="tuning-unpin", kind="unit:link-neighbors")
    assert profileclaim.pinned_mb(conn, "unit:link-neighbors") is None
    _post(web_app, admin_client, action="tuning-clear", kind="unit:link-neighbors")
    assert profileclaim.changes(conn)[0]["what"] == "clear"


def test_a_running_units_own_bar_and_time_left_show_under_its_counter_only_when_it_has_one(admin_client, conn):
    """PERF.147 in the OPS.139 rows: a unit with a shown `progress` draws the bar macro under its counter, one without draws none."""
    run_id = _run(conn, queued=False)
    _publish_unit(conn, run_id, 0, done=300.0, estimate=1000.0,
                  progress={"show": True, "percent": 30.0, "eta_seconds": 120.0, "source": "running"})
    _publish_unit(conn, run_id, 1, done=5.0, estimate=10.0, progress={"show": False, "percent": 50.0, "eta_seconds": None, "source": "running"})
    _publish_unit(conn, run_id, 2, done=1.0, estimate=100.0, progress={"show": True, "percent": 1.0, "eta_seconds": 600.0, "source": "stored"})
    html = admin_client.get(f"/admin/queue/{run_id}").get_data(as_text=True)
    assert html.count('class="task-bar"') == 2
    assert "30%, about 120 s left" in html and "1%, about 600 s left (from earlier runs)" in html
    assert "50%" not in html
