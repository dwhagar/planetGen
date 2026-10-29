# tests/test_web_job_queue.py

"""
The jobs that wait behind the running one in the one waiting list (ADM.123, PERF.154; `planetgen.web.jobs.queue_next`
and `promote_queued`): they start one at a time, in the order they arrived, by themselves when the job before them
succeeded or was cancelled, are held when that job failed, can be started or discarded by an admin from the job strip,
and one more than `queue.max_waiting` is refused naming the list.

Run with: pytest tests/test_web_job_queue.py
"""

import json
import os
import re
import time

import pytest

from planetgen.queue import redisqueue, runline, runrecord
from planetgen.util.settings import get_settings
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
    jobs.release_lock(job_id)


def _running(started):
    first = REAL_START("plan", "Plan", STEP, spawn=False)
    jobs.queue_next("reset", "Reset", STEP, env={"PLANETGEN_MYSQL_PASSWORD": "pw"}, admin="boss", database="d")
    return first


def test_a_job_waits_behind_the_running_one_and_its_environment_stays_out_of_redis(jobs_root, started):
    first = _running(started)
    (queued,) = jobs.queued_jobs()
    assert queued["title"] == "Reset" and queued["after_job"] == first and queued["held"] is None
    assert queued["position"] == 1 and "env" not in queued and "steps" not in queued
    path = os.path.join(jobs_root, f"{jobs.QUEUED_PREFIX}{queued['ticket']}.json")
    assert oct(os.stat(path).st_mode & 0o777) == "0o600"
    assert jobs.promote_queued() is None and started == [], "it waits while the first job runs"
    connection = redisqueue.connect()
    stored = [value for key in connection.scan_iter(match=f"{runrecord.PREFIX}:line:*")
              if connection.type(key) == b"hash" for value in connection.hgetall(key).values()]
    assert stored and not any(b"pw" in value for value in stored), "the database password is only in the private file"


def test_several_wait_in_the_order_they_came_and_start_one_at_a_time(jobs_root, started):
    first = _running(started)
    jobs.queue_next("galaxy", "Second", STEP, admin="boss")
    jobs.queue_next("galaxy", "Third", STEP, admin="boss")
    assert [(job["title"], job["position"]) for job in jobs.queued_jobs()] == [("Reset", 1), ("Second", 2), ("Third", 3)]
    _end(first, "succeeded")
    second = jobs.promote_queued()
    assert started == ["Reset"] and jobs.active_job()["id"] == second
    assert [job["title"] for job in jobs.queued_jobs()] == ["Second", "Third"]
    assert jobs.queued_jobs()[0]["after_job"] == second, "the next one waits for the job that has just started"
    assert jobs.promote_queued() is None, "the next waits for that one"
    _end(second, "succeeded")
    jobs.promote_queued()
    assert started == ["Reset", "Second"]


def test_a_ninth_entry_is_refused_naming_the_list(jobs_root, started, monkeypatch):
    monkeypatch.setattr(get_settings().queue, "max_waiting", 3)
    REAL_START("plan", "Plan", STEP, spawn=False)
    for title in ("One", "Two", "Three"):
        jobs.queue_next("galaxy", title, STEP, admin="boss")
    with pytest.raises(runline.LineFull) as refused:
        jobs.queue_next("galaxy", "Four", STEP, admin="boss")
    assert '"One" (boss)' in str(refused.value) and '"Three" (boss)' in str(refused.value) and "holds 3" in str(refused.value)
    assert [job["title"] for job in jobs.queued_jobs()] == ["One", "Two", "Three"]
    assert not [name for name in os.listdir(jobs_root) if name.startswith(jobs.QUEUED_PREFIX)][3:], "no recipe left for the refused job"


def test_the_default_list_holds_eight(jobs_root, started):
    assert get_settings().queue.max_waiting == 8
    REAL_START("plan", "Plan", STEP, spawn=False)
    for number in range(8):
        jobs.queue_next("galaxy", f"Job {number}", STEP)
    with pytest.raises(runline.LineFull, match="holds 8"):
        jobs.queue_next("galaxy", "Job 8", STEP)


def test_submit_starts_a_job_or_takes_a_place(jobs_root, started):
    first = jobs.submit("plan", "Plan", STEP)
    assert first.job_id and first.ticket is None
    second = jobs.submit("reset", "Reset", STEP)
    assert second.job_id is None and second.position == 1 and second.running == "Plan"
    third = jobs.submit("galaxy", "Galaxy", STEP)
    assert third.position == 2


@pytest.mark.parametrize("status", ["succeeded", "cancelled"])
def test_it_starts_when_the_job_before_it_succeeded_or_was_cancelled(jobs_root, started, status):
    first = _running(started)
    _end(first, status)
    new = jobs.promote_queued()
    assert new and jobs.active_job()["id"] == new and started == ["Reset"]
    assert jobs.queued_jobs() == []
    assert not [name for name in os.listdir(jobs_root) if name.startswith(jobs.QUEUED_PREFIX)], "its recipe is gone"
    with open(os.path.join(jobs_root, new, "job.json")) as f:
        assert json.load(f)["env"] == {"PLANETGEN_MYSQL_PASSWORD": "pw"}, "the database environment is kept"


@pytest.mark.parametrize("status", ["failed", "interrupted"])
def test_it_is_held_when_the_job_before_it_failed_until_an_admin_starts_it(jobs_root, started, status):
    first = _running(started)
    _end(first, status)
    assert jobs.promote_queued() is None and started == []
    assert f"did not succeed ({status})" in jobs.queued_jobs()[0]["held"]
    assert jobs.promote_queued(force=True) and started == ["Reset"] and jobs.queued_jobs() == []


def test_discard_by_ticket_and_a_second_discard_are_safe(jobs_root, started):
    _running(started)
    other = jobs.queue_next("galaxy", "Other", STEP)
    assert jobs.discard_queued(ticket=other) is True and jobs.discard_queued(ticket=other) is False
    assert [job["title"] for job in jobs.queued_jobs()] == ["Reset"]
    assert jobs.discard_queued() is True and jobs.discard_queued() is False
    assert jobs.promote_queued() is None
    assert not [name for name in os.listdir(jobs_root) if name.startswith(jobs.QUEUED_PREFIX)]


def test_discard_all_empties_the_list(jobs_root, started):
    _running(started)
    jobs.queue_next("galaxy", "Other", STEP)
    assert jobs.discard_all_queued() == 2 and jobs.queued_jobs() == []


def test_a_busy_lock_puts_the_job_back_in_its_place(jobs_root, started, monkeypatch):
    first = _running(started)
    jobs.queue_next("galaxy", "Behind", STEP)
    _end(first, "succeeded")

    def busy(*args, **kwargs):
        raise jobs.JobBusy(None)

    monkeypatch.setattr(jobs, "start_job", busy)
    assert jobs.promote_queued() is None
    assert [job["title"] for job in jobs.queued_jobs()] == ["Reset", "Behind"]


def test_a_job_whose_recipe_is_gone_is_held_not_lost(jobs_root, started):
    first = _running(started)
    _end(first, "succeeded")
    (queued,) = jobs.queued_jobs()
    os.remove(os.path.join(jobs_root, f"{jobs.QUEUED_PREFIX}{queued['ticket']}.json"))
    assert jobs.promote_queued() is None
    assert "recipe is gone" in jobs.queued_jobs()[0]["held"]


def test_the_page_queues_a_second_job_and_the_strip_shows_and_starts_it(site, client, no_spawn, jobs_root):
    assert _post(client, action="plan").status_code == 303
    assert _post(client, action="reset", confirm=DB).status_code == 303
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "data-strip-queued" in html and not re.search(r"data-strip-queued[^>]* hidden", html)
    assert "Starts when the running job ends." in html
    (first,) = [job["id"] for job in jobs.list_jobs()]
    _end(first, "failed")
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "did not succeed" in html and 'value="start_queued"' in html
    ticket = jobs.queued_jobs()[0]["ticket"]
    resp = _post(client, action="start_queued", ticket=ticket)
    assert resp.status_code == 303 and jobs.queued_jobs() == [] and jobs.active_job()["title"] == "Reset the galaxy"


def test_the_strip_lists_the_waiting_in_order_with_their_places(site, client, no_spawn, jobs_root):
    _post(client, action="plan")
    _post(client, action="reset", confirm=DB)
    _post(client, action="plan")
    html = client.get("/admin/generate").get_data(as_text=True)
    assert html.index("Reset the galaxy") < html.rindex("Starts after the 1 job ahead of it.")
    assert "Starts after the 1 job ahead of it." in html and "Starts when the running job ends." in html


def test_a_full_list_refuses_the_page_and_a_script_with_the_list_named(site, client, no_spawn, jobs_root, monkeypatch):
    monkeypatch.setattr(get_settings().queue, "max_waiting", 1)
    _post(client, action="plan")
    assert _post(client, action="reset", confirm=DB).status_code == 303
    refused = _post(client, action="plan")
    assert refused.status_code == 409 and "The line is full" in refused.get_data(as_text=True)
    scripted = _post(client, headers={"Accept": "application/json"}, action="plan")
    assert scripted.status_code == 409 and "The line is full" in scripted.get_json()["error"]


def test_a_script_that_finds_a_job_running_is_queued_and_told_its_place(site, client, no_spawn, jobs_root):
    _post(client, action="plan")
    queued = _post(client, headers={"Accept": "application/json"}, action="reset", confirm=DB)
    assert queued.status_code == 202
    body = queued.get_json()
    assert body["queued"] is True and body["position"] == 1 and body["ticket"] and "job" not in body


def test_discard_from_the_strip(site, client, no_spawn, jobs_root):
    _post(client, action="plan")
    _post(client, action="reset", confirm=DB)
    ticket = jobs.queued_jobs()[0]["ticket"]
    resp = _post(client, headers={"Referer": "http://localhost/admin/stats?x=1"}, action="discard_queued", ticket=ticket)
    assert resp.status_code == 303 and resp.headers["Location"] == "/admin/stats?x=1"
    assert jobs.queued_jobs() == []
    resp = _post(client, headers={"Referer": "http://evil.example/admin"}, action="discard_queued")
    assert resp.headers["Location"] == "/admin/generate"
