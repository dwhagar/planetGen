# tests/test_work_queue_admin.py

"""
Tests for managing the work queue from the admin page (ADM.10): pause,
resume and cancel on any node of a job tree as a run sees them
(`WorkQueue._obey`), the whole-queue pause, the stale lease, the load
line (`planetgen/queue/load.py`), the API under `/api/admin/work`
and the pages under `/admin/queue`.

Tests that take the `mysql_config` fixture (see `conftest.py`) are
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import math
import os
import threading
import time

import pytest

from planetgen.web.app import create_app
from planetgen.api.config import Config
from planetgen.queue import load as systemLoad, work as workQueue
from planetgen.db import store
from planetgen.admin import auth as adminAuth

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.web import csrf, jobs, queue_page  # noqa: E402


def _square(payload):
    return payload * payload


def _slow(payload):
    time.sleep(0.2)
    return payload


@pytest.fixture
def control_config(mysql_config, monkeypatch):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    monkeypatch.delenv(workQueue.PARENT_ENV_VAR, raising=False)
    monkeypatch.setattr(workQueue, "WAIT_POLL_SECONDS", 0.05)
    monkeypatch.setattr(workQueue, "CONTROL_POLL_SECONDS", 0.0)
    return mysql_config


def _conn(config):
    return store.get_control_connection(config)


def _call(config, fn, *args):
    conn = _conn(config)
    try:
        return fn(conn, *args)
    finally:
        conn.close()


def _state(config, node_id):
    return _call(config, workQueue.get_node, node_id)["state"]


def _wait_for(check, timeout=20):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if check():
            return
        time.sleep(0.02)
    raise AssertionError("timed out")


# ---------------------------------------------------------------------------
# What a run does when asked
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("workers", [1, 2])
def test_a_paused_queue_node_stands_by_until_resumed(control_config, workers):
    done = []
    with workQueue.WorkQueue("Paused", workers=workers, control_config=control_config) as queue:
        node_id = queue.job_id
        assert _call(control_config, workQueue.request_control, node_id, "pause")

        def resume_later():
            _wait_for(lambda: _state(control_config, node_id) == "paused")
            if workers > 1:
                # Standby gives the lease up, so other runs can go.
                lease = _call(control_config, workQueue.queue_status)
                assert lease["holder"] is None
            _call(control_config, workQueue.request_control, node_id, "resume")
            done.append("resumed")

        thread = threading.Thread(target=resume_later)
        thread.start()
        for n in range(3):
            queue.submit("square", str(n), _square, n)
        thread.join()
    assert done == ["resumed"] and queue.finished == 3
    assert _state(control_config, node_id) == "done"


def test_cancelling_a_queue_node_skips_its_tasks_and_the_run_goes_on(control_config):
    with workQueue.job_node("galaxy", "Run", control_config) as root:
        with workQueue.WorkQueue("First batch", workers=1, control_config=control_config) as first:
            first.submit("sector", "0", _square, 0)
            assert _call(control_config, workQueue.request_control, first.job_id, "cancel")
            for n in range(1, 5):
                first.submit("sector", str(n), _square, n)
        with workQueue.WorkQueue("Second batch", workers=1, control_config=control_config) as second:
            second.submit("sector", "9", _square, 9)
    assert first.finished == 1 and second.finished == 1
    tree = _call(control_config, workQueue.load_tree, root.id)
    assert [child["state"] for child in tree["children"]] == ["cancelled", "done"]
    assert tree["state"] == "done"


@pytest.mark.parametrize("workers", [1, 2])
def test_cancelling_a_run_ends_it_from_its_queue(control_config, workers):
    with pytest.raises(workQueue.Cancelled) as raised:
        with workQueue.job_node("galaxy", "Run", control_config) as root:
            with workQueue.WorkQueue("Batch", workers=workers, control_config=control_config) as queue:
                queue.submit("sector", "0", _slow, 0)
                assert _call(control_config, workQueue.request_control, root.id, "cancel")
                for n in range(1, 6):
                    queue.submit("sector", str(n), _slow, n)
    assert raised.value.code == workQueue.CANCELLED_EXIT_CODE
    assert queue.finished < 6
    tree = _call(control_config, workQueue.load_tree, root.id)
    assert tree["state"] == "cancelled" and tree["children"][0]["state"] == "cancelled"
    assert _call(control_config, workQueue.queue_status)["holder"] is None


def test_a_cancelled_run_starts_no_new_phase(control_config):
    with pytest.raises(workQueue.Cancelled):
        with workQueue.job_node("plan", "Plan", control_config) as root:
            _call(control_config, workQueue.request_control, root.id, "cancel")
            with workQueue.job_node("bright-stars", "Bright stars"):
                raise AssertionError("should not start")


def test_pausing_the_whole_queue_holds_every_run(control_config, redis_server):
    assert _call(control_config, workQueue.pause_queue, "boss")
    assert not _call(control_config, workQueue.pause_queue, "boss")
    waits = []

    def resume_later():
        _wait_for(lambda: waits)
        time.sleep(0.2)
        _call(control_config, workQueue.resume_queue)

    thread = threading.Thread(target=resume_later)
    thread.start()
    started = time.monotonic()
    with workQueue.WorkQueue("Waits", workers=2, control_config=control_config, on_wait=waits.append) as queue:
        queue.submit("square", "1", _square, 1)
    thread.join()
    assert waits[0]["paused"] and waits[0]["paused_by"] == "boss"
    assert time.monotonic() - started >= 0.2
    assert not _call(control_config, workQueue.queue_status)["paused"]


def test_a_one_worker_run_also_waits_for_a_paused_queue(control_config):
    _call(control_config, workQueue.pause_queue, "boss")
    timer = threading.Timer(0.3, lambda: _call(control_config, workQueue.resume_queue))
    timer.start()
    started = time.monotonic()
    with workQueue.WorkQueue("Serial", workers=1, control_config=control_config) as queue:
        queue.submit("square", "1", _square, 1)
    timer.join()
    assert time.monotonic() - started >= 0.25


def test_controls_only_reach_live_nodes(control_config):
    with workQueue.job_node("galaxy", "Done soon", control_config) as root:
        pass
    assert not _call(control_config, workQueue.request_control, root.id, "pause")
    assert not _call(control_config, workQueue.request_control, root.id, "cancel")
    with pytest.raises(ValueError):
        _call(control_config, workQueue.request_control, root.id, "explode")


def test_resuming_a_node_resumes_its_paused_subtree(control_config):
    with workQueue.job_node("galaxy", "Run", control_config) as root:
        with workQueue.job_node("population", "Phase") as phase:
            _call(control_config, workQueue.request_control, root.id, "pause")
            _call(control_config, workQueue.request_control, phase.id, "pause")
            assert _call(control_config, workQueue.request_control, root.id, "resume")
            assert _call(control_config, workQueue.get_node, phase.id)["control"] is None


def test_stale_lease_is_cleared_and_live_one_kept(control_config):
    conn = _conn(control_config)
    try:
        with conn:
            conn.execute("INSERT INTO work_jobs (id, root_id, title, holder, state, workers, created_at,"
                         " heartbeat_at) VALUES ('20261001-000000-0000dead', '20261001-000000-0000dead', 'Died',"
                         " 'gone:1:x', 'running', 2, NOW(6), NOW(6) - INTERVAL 5 MINUTE)")
            conn.execute("INSERT INTO work_lease (id, holder, job_id, heartbeat_at) VALUES"
                         " (1, 'gone:1:x', '20261001-000000-0000dead', NOW(6) - INTERVAL 5 MINUTE)")
    finally:
        conn.close()
    status = _call(control_config, workQueue.queue_status)
    assert status["stale"] and status["holder"] == "gone:1:x" and status["heartbeat_age_s"] > 200
    assert _call(control_config, workQueue.clear_stale_lease) == "gone:1:x"
    assert _call(control_config, workQueue.queue_status)["holder"] is None
    assert _state(control_config, "20261001-000000-0000dead") == "cancelled"
    assert _call(control_config, workQueue.clear_stale_lease) is None


def test_workers_active_counts_running_queues(control_config, redis_server):
    with workQueue.job_node("galaxy", "Run", control_config):
        with workQueue.WorkQueue("Pool", workers=3, control_config=control_config):
            status = _call(control_config, workQueue.queue_status)
            assert (status["workers_active"], status["runs_active"]) == (3, 1)
    assert _call(control_config, workQueue.queue_status)["workers_active"] == 0


def test_only_finished_roots_are_deleted(control_config):
    with workQueue.job_node("galaxy", "Run", control_config) as root:
        with workQueue.job_node("population", "Phase") as phase:
            assert not _call(control_config, workQueue.delete_tree, root.id)
    assert not _call(control_config, workQueue.delete_tree, phase.id)
    assert _call(control_config, workQueue.delete_tree, root.id)
    assert _call(control_config, workQueue.get_node, phase.id) is None


# ---------------------------------------------------------------------------
# The load line
# ---------------------------------------------------------------------------

def test_the_load_line_on_this_system():
    load = systemLoad.load_average()
    if hasattr(os, "getloadavg"):
        assert load["kind"] == "load" and len(load["values"]) == 3
        assert load["text"].count(" / ") == 2


def test_cpu_percent_averages_decay_like_the_load_average():
    readings = iter([(0, 0), (50, 100), (50, 200), (50, 300)])
    clock = iter([0.0, 5.0, 10.0, 15.0])
    sampler = systemLoad.CpuSampler(lambda: next(readings), clock=lambda: next(clock))
    sampler._thread = object()  # no background thread in the test
    assert sampler.averages() is None
    sampler.sample()
    sampler.sample()
    assert sampler.averages() == [50.0, 50.0, 50.0]  # first interval: half idle
    sampler.sample()
    one, five, fifteen = sampler.averages()
    assert one == pytest.approx(50 + 50 * (1 - math.exp(-5 / 60)))
    assert 50 < fifteen < five < one < 100


def test_windows_shows_cpu_percent(monkeypatch):
    monkeypatch.delattr(os, "getloadavg", raising=False)
    monkeypatch.setattr(os, "name", "nt")

    class Fake:
        def averages(self):
            return [12.4, 9.0, 8.2]

    monkeypatch.setattr(systemLoad, "_sampler", lambda: Fake())
    assert systemLoad.load_average() == {"values": [12.4, 9.0, 8.2], "kind": "cpu", "text": "12% / 9% / 8%"}


# ---------------------------------------------------------------------------
# The API and the pages
# ---------------------------------------------------------------------------

@pytest.fixture
def web_app(control_config, monkeypatch, tmp_path):
    class RealConfig(Config):
        MYSQL_CONFIG = control_config
        WRITE_MYSQL_CONFIG = control_config
        CONTROL_MYSQL_CONFIG = control_config
        WEB_DATABASE = control_config.database
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))

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


def _post(app, client, **form):
    nonce = "n" * 43
    client.set_cookie(csrf.COOKIE_NAME, nonce)
    session = client.get_cookie("pg_admin_session")
    with app.app_context():
        form[csrf.FIELD_NAME] = csrf._sign(nonce, session.value if session else "")
    return client.post("/admin/queue/action", data=form)


def _audit(config):
    conn = _conn(config)
    try:
        return [row["action"] for row in conn.execute("SELECT action FROM admin_audit_log ORDER BY id").fetchall()]
    finally:
        conn.close()


def test_the_api_needs_an_admin(web_app):
    client = web_app.test_client()
    assert client.get("/api/admin/work").status_code == 401
    assert client.post("/api/admin/work/queue", json={"action": "pause"}).status_code == 401


def test_the_api_lists_trees_and_controls_them(admin_client, control_config):
    with workQueue.job_node("galaxy", "Live run", control_config, argv=["galaxy", "--ring", "1"]) as root:
        with workQueue.WorkQueue("Batch", workers=1, control_config=control_config) as queue:
            queue.expect(3)
            queue.submit("sector", "1,0,0", _square, 1)
            body = admin_client.get("/api/admin/work").get_json()
            assert body["available"] and body["total"] == 1
            item = body["items"][0]
            assert item["id"] == root.id and item["live"] and item["totals"]["tasks"] == 3
            assert "children" not in item
            assert body["status"]["runs_active"] == 1 and body["status"]["load"]["text"]
            tree = admin_client.get(f"/api/admin/work/{queue.job_id}").get_json()["tree"]
            assert tree["id"] == root.id and tree["children"][0]["tasks"][0]["task_key"] == "1,0,0"
            response = admin_client.post(f"/api/admin/work/{queue.job_id}/control", json={"action": "cancel"})
            assert response.get_json() == {"ok": True}
            queue.submit("sector", "1,0,1", _square, 2)
    assert queue.finished == 1
    assert admin_client.post(f"/api/admin/work/{root.id}/control", json={"action": "bogus"}).status_code == 400
    assert admin_client.get("/api/admin/work/not-an-id").status_code == 404
    assert admin_client.get("/api/admin/work/20260101-000000-00000000").status_code == 404
    assert admin_client.post("/api/admin/work/queue", json={"action": "pause"}).get_json() == {"ok": True}
    assert admin_client.post("/api/admin/work/queue", json={"action": "resume"}).get_json() == {"ok": True}
    assert admin_client.post(f"/api/admin/work/{root.id}/delete").get_json() == {"ok": True}
    assert _audit(control_config)[-4:] == ["work.cancel", "work.queue.pause", "work.queue.resume", "work.delete"]


def test_the_queue_page_shows_workers_load_and_jobs(admin_client, control_config):
    with workQueue.job_node("galaxy", "Fill ring 1", control_config):
        with workQueue.WorkQueue("Sectors (ring 1)", workers=2, control_config=control_config):
            html = admin_client.get("/admin/queue").get_data(as_text=True)
    assert "Workers active" in html and "Load average" in html
    assert "Fill ring 1" in html and "Pause the queue" in html
    assert 'href="/admin/queue/confirm/pause?node=' in html


def test_the_jobs_table_route_lists_jobs_with_their_controls(web_app, admin_client, control_config):
    assert web_app.test_client().get("/table/queue-jobs").status_code == 403
    with workQueue.job_node("galaxy", "Fill ring 1", control_config) as root:
        with workQueue.WorkQueue("Sectors (ring 1)", workers=1, control_config=control_config):
            data = admin_client.get("/table/queue-jobs").get_json()
    assert data["total"] == 1
    title, kind, status, started, _duration, _progress, actions = data["rows"][0]
    assert title["text"] == "Fill ring 1" and title["href"].startswith(f"/admin/queue/{root.id}")
    labels = {part["text"]: part["href"] for part in actions["parts"] if isinstance(part, dict)}
    assert labels["Cancel"].startswith("/admin/queue/confirm/cancel?node=")


def test_the_tree_page_expands_down_to_tasks(admin_client, control_config):
    with workQueue.job_node("galaxy", "Fill ring 1", control_config) as root:
        with workQueue.WorkQueue("Sectors (ring 1)", workers=1, control_config=control_config) as queue:
            queue.submit("sector", "1,0,0", _square, 1)
    html = admin_client.get(f"/admin/queue/{queue.job_id}").get_data(as_text=True)
    assert f'id="node-{root.id}"' in html and f'id="node-{queue.job_id}"' in html
    assert "<code>1,0,0</code>" in html and "Delete" in html
    assert admin_client.get("/admin/queue/20260101-000000-00000000").status_code == 404


def test_every_control_is_confirmed_first(web_app, admin_client, control_config):
    with workQueue.job_node("galaxy", "Fill ring 1", control_config) as root:
        page = admin_client.get(f"/admin/queue/confirm/pause?node={root.id}").get_data(as_text=True)
        assert "Pause “Fill ring 1” and everything under it?" in page
        assert 'name="action" value="pause"' in page
        assert _call(control_config, workQueue.get_node, root.id)["control"] is None
        response = _post(web_app, admin_client, action="pause", node=root.id)
        assert response.status_code == 303
        assert _call(control_config, workQueue.get_node, root.id)["control"] == "pause"
        assert "Asked" in admin_client.get(response.headers["Location"]).get_data(as_text=True)
        _post(web_app, admin_client, action="resume", node=root.id)
        assert _call(control_config, workQueue.get_node, root.id)["control"] is None
    assert _post(web_app, admin_client, action="queue-pause").status_code == 303
    assert _call(control_config, workQueue.queue_status)["paused"]
    assert "The queue is paused" in admin_client.get("/admin/queue").get_data(as_text=True)
    _post(web_app, admin_client, action="queue-resume")
    assert not _call(control_config, workQueue.queue_status)["paused"]


def test_the_action_needs_csrf(admin_client):
    assert admin_client.post("/admin/queue/action", data={"action": "queue-pause"}).status_code == 400


def test_retry_reruns_a_cli_run_and_a_failed_sector(web_app, admin_client, control_config, monkeypatch):
    started = []
    monkeypatch.setattr(jobs, "start_job", lambda kind, title, steps, **kw: started.append((kind, title, steps, kw))
                        or "20261001-000000-abcd")
    with pytest.raises(ValueError):
        with workQueue.job_node("galaxy", "planetgen galaxy --ring 2", control_config,
                                argv=["galaxy", "--ring", "2"], database=control_config.database) as root:
            with workQueue.WorkQueue("Sectors", workers=1, control_config=control_config) as queue:
                queue.submit("sector", "2,0,5", _broken, 0)
    page = admin_client.get(f"/admin/queue/{root.id}").get_data(as_text=True)
    assert "/admin/queue/confirm/retry?node=" in page
    _post(web_app, admin_client, action="retry", node=root.id)
    kind, title, steps, kw = started[-1]
    assert kind == "galaxy" and steps[0]["argv"][-3:] == ["galaxy", "--ring", "2"]
    assert kw["database"] == control_config.database
    task_id = _call(control_config, workQueue.load_tree, root.id)["children"][0]["tasks"][0]["id"]
    confirm = admin_client.get(f"/admin/queue/confirm/retry?node={queue.job_id}&task={task_id}")
    assert "Retry sector 2,0,5?" in confirm.get_data(as_text=True)
    _post(web_app, admin_client, action="retry", node=queue.job_id, task=str(task_id))
    assert started[-1][2][0]["argv"][-7:] == ["galaxy", "--ring", "2", "--layer", "0", "--slot", "5"]


def test_retry_refuses_another_database(web_app, admin_client, control_config, monkeypatch):
    monkeypatch.setattr(jobs, "start_job", lambda *a, **k: pytest.fail("should not start"))
    with pytest.raises(ValueError):
        with workQueue.job_node("galaxy", "Elsewhere", control_config, argv=["galaxy", "--ring", "2"],
                                database="some_other_galaxy") as root:
            raise ValueError("boom")
    response = _post(web_app, admin_client, action="retry", node=root.id)
    assert "some_other_galaxy" in admin_client.get(response.headers["Location"]).get_data(as_text=True)


def _broken(payload):
    raise ValueError("sector broke")


def test_retry_plan_for_a_web_job_starts_at_the_step_that_failed(tmp_path, monkeypatch):
    root_dir = tmp_path / "jobs"
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(root_dir))
    job_id = jobs.start_job("new_galaxy", "New galaxy", [
        {"label": "Reset", "argv": ["x"]}, {"label": "Plan", "argv": ["y"]}, {"label": "Fill", "argv": ["z"]},
    ], spawn=False)
    path = root_dir / job_id
    (path / "state.json").write_text('{"status": "failed", "step": 2, "started_at": 1, "finished_at": 2}')
    (root_dir / jobs.LOCK_NAME).unlink()
    job, steps = jobs.remaining_steps(job_id)
    assert job["status"] == "failed" and [step["label"] for step in steps] == ["Plan", "Fill"]
    root = {"kind": "web-job", "web_job_id": job_id, "title": "New galaxy", "argv": None, "database_name": None}
    plan = queue_page._retry_plan(root)
    assert plan["describe"] == "Plan, then Fill" and plan["title"] == "Retry: New galaxy"
    (path / "state.json").write_text('{"status": "succeeded", "step": 3, "started_at": 1, "finished_at": 2}')
    assert jobs.remaining_steps(job_id) is None
