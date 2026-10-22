# tests/test_web_job_queue.py

"""
The one job that waits behind the running one (ADM.123, `planetgen.web.jobs.queue_next` and `promote_queued`): it
starts by itself when the job before it succeeded or was cancelled, is held when that job failed, can be started
or discarded by an admin from the job strip, and a second one is refused.

Run with: pytest tests/test_web_job_queue.py
"""

import json
import os
import re
import time

import pytest

from planetgen.web import jobs
from tests.test_web_generate import DB, _post, client, jobs_root, no_spawn, site  # noqa: F401 -- fixtures

STEP = [{"label": "x", "argv": ["true"]}]
REAL_START = jobs.start_job


@pytest.fixture
def started(monkeypatch):
    """Jobs the queue starts, written without running them."""
    seen = []

    def fake_start(kind, title, steps, **kwargs):
        seen.append(title)
        return REAL_START(kind, title, steps, **{**kwargs, "spawn": False})

    monkeypatch.setattr(jobs, "start_job", fake_start)
    return seen


def _end(job_id, status):
    """Marks a job over with `status`, as its runner does, and frees the lock."""
    path = os.path.join(jobs.jobs_dir(), job_id)
    with open(os.path.join(path, "state.json"), "w") as f:
        json.dump({"status": status, "pid": None, "step": 1, "started_at": time.time() - 5,
                   "finished_at": time.time(), "exit_code": 0 if status == "succeeded" else 1}, f)
    jobs._release_lock(jobs.jobs_dir(), job_id)


def _running(started):
    first = REAL_START("plan", "Plan", STEP, spawn=False)
    jobs.queue_next("reset", "Reset", STEP, env={"PLANETGEN_MYSQL_PASSWORD": "pw"}, admin="boss", database="d")
    return first


def test_a_job_waits_behind_the_running_one_and_only_one_can(jobs_root, started):
    first = _running(started)
    queued = jobs.queued_job()
    assert queued["title"] == "Reset" and queued["after_job"] == first and queued["held"] is None
    assert "env" not in queued and "steps" not in queued
    assert oct(os.stat(os.path.join(jobs_root, jobs.QUEUED_NAME)).st_mode & 0o777) == "0o600"
    with pytest.raises(jobs.QueueFull):
        jobs.queue_next("reset", "Again", STEP)
    assert jobs.promote_queued() is None and started == [], "it waits while the first job runs"


@pytest.mark.parametrize("status", ["succeeded", "cancelled"])
def test_it_starts_when_the_job_before_it_succeeded_or_was_cancelled(jobs_root, started, status):
    first = _running(started)
    _end(first, status)
    new = jobs.promote_queued()
    assert new and jobs.active_job()["id"] == new and started == ["Reset"]
    assert jobs.queued_job() is None
    with open(os.path.join(jobs_root, new, "job.json")) as f:
        assert json.load(f)["env"] == {"PLANETGEN_MYSQL_PASSWORD": "pw"}, "the database environment is kept"


@pytest.mark.parametrize("status", ["failed", "interrupted"])
def test_it_is_held_when_the_job_before_it_failed_until_an_admin_starts_it(jobs_root, started, status):
    first = _running(started)
    _end(first, status)
    assert jobs.promote_queued() is None and started == []
    assert f"did not succeed ({status})" in jobs.queued_job()["held"]
    assert jobs.promote_queued(force=True) and started == ["Reset"] and jobs.queued_job() is None


def test_discard_and_a_second_start_are_safe(jobs_root, started):
    _running(started)
    assert jobs.discard_queued() is True and jobs.discard_queued() is False
    assert jobs.promote_queued() is None


def test_a_busy_lock_puts_the_job_back(jobs_root, started, monkeypatch):
    first = _running(started)
    _end(first, "succeeded")

    def busy(*args, **kwargs):
        raise jobs.JobBusy(None)

    monkeypatch.setattr(jobs, "start_job", busy)
    assert jobs.promote_queued() is None and jobs.queued_job()["title"] == "Reset"
    assert not os.path.exists(os.path.join(jobs_root, jobs.QUEUED_NAME + ".starting"))


def test_the_page_queues_a_second_job_and_the_strip_shows_and_starts_it(site, client, no_spawn, jobs_root):
    assert _post(client, action="plan").status_code == 303
    assert _post(client, action="reset", confirm=DB).status_code == 303
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "data-strip-queued" in html and not re.search(r"data-strip-queued[^>]* hidden", html)
    assert "Starts when the running job ends." in html
    (first,) = [job["id"] for job in jobs.list_jobs()]
    _end(first, "failed")
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "did not succeed" in html and not re.search(r"data-strip-start[^>]* hidden", html)
    resp = _post(client, action="start_queued")
    assert resp.status_code == 303 and jobs.queued_job() is None and jobs.active_job()["title"] == "Reset the galaxy"


def test_discard_from_the_strip(site, client, no_spawn, jobs_root):
    _post(client, action="plan")
    _post(client, action="reset", confirm=DB)
    resp = _post(client, headers={"Referer": "http://localhost/admin/stats?x=1"}, action="discard_queued")
    assert resp.status_code == 303 and resp.headers["Location"] == "/admin/stats?x=1"
    assert jobs.queued_job() is None
    resp = _post(client, headers={"Referer": "http://evil.example/admin"}, action="discard_queued")
    assert resp.headers["Location"] == "/admin/generate"
