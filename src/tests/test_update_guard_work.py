# tests/test_update_guard_work.py

"""
OPS.53: before an update changes anything it lists the work that runs and stops unless the operator agrees to
cancel it, and from then on nothing new starts. Fake work for the decisions (`tests/test_update_guard.py` has
`FakeWork`), real processes, a real Generate job directory, the work queue's rows and a Redis queue for the
sources.
"""

import math
import os
import subprocess
import sys
import time

import pytest
import rq

from planetgen.admin import maintenance, updateguard
from planetgen.db import store
from planetgen.queue import api_jobs, redisqueue, work, workers
from tests.test_update_guard import FakeWork, begin, lines, say  # noqa: F401 -- fixtures


def sure_nothing_changed(fake=None):
    assert not maintenance.is_on() and updateguard.read_lock() is None
    if fake is not None:
        assert fake.cancelled == []


# --- the decision ---------------------------------------------------------------------------------------------

def test_with_nothing_running_the_update_goes_on_without_asking(say):
    asked = []
    assert begin(sources=[FakeWork()], interactive=True, answer=lambda text: asked.append(text) or "n")
    assert asked == []


def test_running_work_is_listed_and_without_a_terminal_the_update_refuses(say, lines):
    fake = FakeWork("Generate sectors 12-14", "Orbital update", queued=("Orbital update",))
    assert not begin(sources=[fake], write=say)
    text = "\n".join(lines)
    assert "Generate sectors 12-14 (generate job-0, running, started 2026-10-16T03:00:00Z; holds the work lease)" in text
    assert "Orbital update (generate job-1, queued" in text
    assert "--cancel-running" in text and "nothing was changed" in text
    sure_nothing_changed(fake)


def test_the_flag_is_the_operators_word_to_cancel(say, lines):
    fake = FakeWork("Generate sectors", "Warm-up")
    assert begin(sources=[fake], cancel_running=True, write=say)
    assert sorted(fake.cancelled) == ["job-0", "job-1"] and fake.items() == []
    assert maintenance.read()["cause"] == updateguard.UPDATE_CAUSE and updateguard.read_lock() is not None
    assert "Everything stopped" in "\n".join(lines)


def test_on_a_terminal_a_no_changes_nothing_and_a_yes_cancels(say):
    fake = FakeWork("Generate sectors")
    questions = []
    assert not begin(sources=[fake], interactive=True, answer=lambda text: questions.append(text) or "n", write=say)
    assert len(questions) == 1 and "Cancel these" in questions[0]
    sure_nothing_changed(fake)
    assert begin(sources=[fake], interactive=True, answer=lambda text: "Yes", write=say)
    assert fake.cancelled == ["job-0"]


def test_a_queued_job_is_listed_and_cancelled_like_a_running_one(say, lines):
    fake = FakeWork("Backfill", queued=("Backfill",))
    assert not begin(sources=[fake], write=say)
    assert "queued" in "\n".join(lines)
    assert begin(sources=[fake], cancel_running=True, write=say)
    assert fake.cancelled == ["job-0"]


def test_a_job_that_does_not_stop_in_time_is_reported_and_the_update_stops_before_any_change(say, lines):
    fake = FakeWork("Stubborn run", stops=False)
    assert not begin(sources=[fake], cancel_running=True, write=say, timeout=0.1)
    text = "\n".join(lines)
    assert "did not stop within" in text and "Stubborn run" in text and "nothing was changed" in text
    assert fake.cancelled == ["job-0"]       # asked once, not every poll
    sure_nothing_changed()


def test_work_that_starts_between_the_first_look_and_the_lock_is_not_cancelled_without_a_word(say, lines):
    class Late(FakeWork):
        looks = 0

        def items(self):
            Late.looks += 1
            return [] if Late.looks == 1 else super().items()

    late = Late("Started a moment ago")
    assert not begin(sources=[late], write=say)
    assert "started while the update began" in "\n".join(lines)
    sure_nothing_changed(late)


def test_work_that_starts_after_the_word_is_cancelled_too(say):
    class Late(FakeWork):
        looks = 0

        def items(self):
            Late.looks += 1
            return super().items() if Late.looks > 1 else []

    late = Late("Started a moment ago")
    assert begin(sources=[late], cancel_running=True, write=say)
    assert late.cancelled == ["job-0"]


# --- nothing new starts during the update --------------------------------------------------------------------

@pytest.fixture
def updating():
    assert begin(sources=[])
    yield
    updateguard.end(os.getpid(), ok=True)


def test_a_generate_job_offered_during_the_update_is_not_started(updating, tmp_path):
    from planetgen.web import jobs
    with pytest.raises(updateguard.UpdateRunning, match="A Generate job can't start"):
        jobs.start_job("galaxy", "Generate", [{"label": "x", "argv": ["true"]}], root=str(tmp_path), spawn=False)
    assert os.listdir(tmp_path) == []


def test_a_queue_job_offered_during_the_update_is_not_queued(updating):
    with pytest.raises(api_jobs.NoQueue, match="the site is being updated"):
        api_jobs.submit(math.pow, 2, 2)


def test_a_run_offered_during_the_update_does_not_start(updating):
    with pytest.raises(updateguard.UpdateRunning, match="A run can't start"):
        work.open_node("plan", "A run")
    assert work.current_node() is None


def test_the_workers_take_no_job_during_the_update(updating, capsys):
    from planetgen.cli import worker
    assert workers.ensure_workers("planetgen-bulk", connection=object(), url="redis://never") == 0
    assert worker.main(["planetgen-bulk", "--burst"]) == updateguard.EX_TEMPFAIL
    assert "A worker can't start" in capsys.readouterr().err


def test_a_site_a_failed_update_left_closed_refuses_the_same_work(tmp_path):
    from planetgen.web import jobs
    assert begin(sources=[])
    updateguard.end(os.getpid(), ok=False, step="step 5/9", write=lambda text: None)
    with pytest.raises(updateguard.UpdateRunning):
        jobs.start_job("galaxy", "Generate", [{"label": "x", "argv": ["true"]}], root=str(tmp_path), spawn=False)


# --- the real sources -------------------------------------------------------------------------------------------

REGISTER = ("import sys, time; from planetgen.admin import updateguard\n"
            "with updateguard.running('orbits', 'Orbital update'):\n"
            "    print('registered', flush=True); time.sleep(120)\n")


def registered_process():
    child = subprocess.Popen([sys.executable, "-c", REGISTER], stdout=subprocess.PIPE, text=True,
                             env={**os.environ})
    assert child.stdout.readline().strip() == "registered"
    return child


def test_a_registered_run_is_listed_stopped_and_forgotten(say):
    child = registered_process()
    try:
        source = updateguard.MarkerSource()
        found = source.items()
        assert [(item["kind"], item["pid"]) for item in found] == [("orbits", child.pid)]
        assert not begin(sources=[source], write=say)
        assert child.poll() is None
        assert begin(sources=[source], cancel_running=True, write=say, timeout=30, poll=0.1)
        assert child.wait(timeout=30) != 0 and source.items() == []
    finally:
        child.kill()
        child.wait()


def test_a_marker_of_a_process_that_died_is_ignored_and_removed():
    child = registered_process()
    child.kill()
    child.wait()
    source = updateguard.MarkerSource()
    assert source.items() == []
    assert os.listdir(updateguard.running_dir()) == []


def test_a_generate_job_is_listed_with_its_title_and_cancelled_through_its_own_cancel_file(tmp_path, monkeypatch):
    from planetgen.web import jobs
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path))
    runner = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)", "fake-runner"])
    try:
        job_id = jobs.start_job("galaxy", "Generate sectors 1-9", [{"label": "step", "argv": ["true"]}],
                                root=str(tmp_path), spawn=False)
        with open(os.path.join(tmp_path, job_id, "state.json"), "w", encoding="utf-8") as handle:
            handle.write('{"status": "running", "pid": %d, "started_at": %f, "step": 1}' % (runner.pid, time.time()))
        # the job's runner is recognised by its command line naming the job
        monkeypatch.setattr(jobs, "_runner_alive", lambda pid, name: pid == runner.pid)
        source = updateguard.GenerateJobSource()
        found = source.items()
        assert [(item["id"], item["title"], item["queued"]) for item in found] == [(job_id, "Generate sectors 1-9", False)]
        assert source.cancel(found[0]) is True
        assert os.path.exists(os.path.join(tmp_path, job_id, jobs.CANCEL_NAME))
    finally:
        runner.kill()
        runner.wait()


def test_no_jobs_directory_is_listed_as_nothing_and_not_created(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "not-there"))
    assert updateguard.GenerateJobSource().items() == []
    assert not (tmp_path / "not-there").exists()


@pytest.fixture
def control_config(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.close()
    return mysql_config


def test_a_work_queue_run_is_listed_and_cancelled_through_its_control_request(control_config):
    source = updateguard.WorkRunSource(control_config)
    assert source.items() == []
    node = work.open_node("plan", "Generate sectors by command", control_config=control_config)
    try:
        found = source.items()
        assert [(item["kind"], item["title"], item["queued"]) for item in found] == [
            ("plan", "Generate sectors by command", False)]
        assert source.cancel(found[0]) is True
        conn = store.get_control_connection(control_config)
        try:
            assert conn.execute("SELECT control FROM work_jobs WHERE id = ?", (node.id,)).fetchone()["control"] == "cancel"
        finally:
            conn.close()
    finally:
        work.close_node(node, "cancelled")
    assert source.items() == []


def test_a_run_of_a_generate_job_is_left_to_the_job(control_config):
    node = work.open_node("web-job", "Generate", control_config=control_config, web_job_id="20261016-000000-abcd")
    try:
        assert updateguard.WorkRunSource(control_config).items() == []
    finally:
        work.close_node(node, "done")


def test_without_a_control_database_no_run_is_listed():
    assert updateguard.WorkRunSource(None).items() == []


def test_queued_and_running_jobs_on_the_shared_queues_are_listed_and_cancelled(redis_server, monkeypatch):
    unique = f"planetgen-update-guard-test-{os.getpid()}"
    monkeypatch.setattr(redisqueue, "QUEUES", (unique,))
    connection = redisqueue.connect(redis_server)
    queue = redisqueue.queue(unique, connection)
    waiting = queue.enqueue(math.pow, 2, 2)
    running = queue.enqueue(math.pow, 2, 3)
    queue.remove(running)
    running.set_status(rq.job.JobStatus.STARTED)
    connection.zadd(queue.started_job_registry.key, {running.id: int(time.time()) + 60})
    stopped = []
    monkeypatch.setattr(workers.command, "send_stop_job_command", lambda conn, job_id: stopped.append(job_id))
    try:
        source = updateguard.QueuedJobSource(redis_server)
        found = {item["id"]: item for item in source.items()}
        assert set(found) == {waiting.id, running.id}
        assert found[waiting.id]["queued"] is True and found[running.id]["queued"] is False
        for item in found.values():
            source.cancel(item)
        assert waiting.get_status(refresh=True) == "canceled" and stopped == [running.id]
    finally:
        queue.empty()
        connection.zrem(queue.started_job_registry.key, running.id)
