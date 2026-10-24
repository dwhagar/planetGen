# tests/test_reset_jobs.py

"""
A job that a killed worker left "running" in Redis no longer holds the update for ever (`workers.reap_dead_jobs`,
called by the update guard and at every worker's start), and `planetgen.cli.reset_jobs` (with `update.sh
--reset-jobs`) resets every job and the queue.
"""

import math
import os
import subprocess
import sys
import uuid

import pytest
import rq

from planetgen.admin import updateguard
from planetgen.cli import reset_jobs
from planetgen.db import store
from planetgen.queue import pause, redisqueue, runrecord, work, workers
from tests.test_update_guard import FakeWork, REPO, run_shell

pytestmark = pytest.mark.real_queue

ALL_FOREVER = "+inf"
"""str: The score of a job queued with no time limit in the started registry."""


@pytest.fixture
def queue_of_its_own(redis_server, monkeypatch):
    """A queue name of this test's own on the test server, emptied after; `(queue, connection)`."""
    name = f"planetgen-api-test-{uuid.uuid4().hex[:10]}"
    connection = redisqueue.connect(redis_server)
    queue = redisqueue.queue(name, connection)
    queue.registered_test_workers = []
    connection.sadd("rq:queues", queue.key)
    yield queue, connection
    for worker in queue.registered_test_workers:
        worker.register_death()
    queue.empty()
    connection.delete(queue.started_job_registry.key, queue.failed_job_registry.key)
    connection.srem("rq:queues", f"rq:queue:{name}")


def dead_pid():
    child = subprocess.Popen([sys.executable, "-c", "pass"])
    child.wait()
    return child.pid


def register_worker(queue, connection, name, pid):
    worker = rq.Worker([queue], connection=connection, name=name, serializer=queue.serializer)
    worker.pid = pid
    worker.register_birth()
    getattr(queue, "registered_test_workers", []).append(worker)
    return worker


def running_job(queue, connection, owner, retries=None):
    """A job the worker `owner` started and never ended, in the started registry with no expiry."""
    job = queue.enqueue(math.pow, 2, 2, job_timeout=600, retry=retries)
    queue.remove(job)
    job.set_status(rq.job.JobStatus.STARTED)
    job.worker_name = owner
    job.save()
    connection.zadd(queue.started_job_registry.key, {f"{job.id}:{uuid.uuid4().hex[:8]}": ALL_FOREVER})
    return job


def test_a_job_whose_worker_process_is_gone_is_failed_and_leaves_the_started_registry(queue_of_its_own):
    queue, connection = queue_of_its_own
    gone = f"gone-{uuid.uuid4().hex[:8]}"
    register_worker(queue, connection, gone, dead_pid())
    job = running_job(queue, connection, gone)
    assert workers.reap_dead_jobs(connection) == [(queue.name, job.id)]
    assert job.get_status(refresh=True) == "failed"
    assert queue.started_job_registry.get_job_ids(cleanup=False) == []
    assert job.id in queue.failed_job_registry.get_job_ids()
    assert workers.reap_dead_jobs(connection) == []


def test_a_job_whose_worker_key_has_expired_is_cleared_too(queue_of_its_own):
    queue, connection = queue_of_its_own
    job = running_job(queue, connection, "never-registered")
    assert workers.reap_dead_jobs(connection) == [(queue.name, job.id)]
    assert job.get_status(refresh=True) == "failed"


def test_a_job_with_a_retry_left_is_queued_again(queue_of_its_own):
    queue, connection = queue_of_its_own
    job = running_job(queue, connection, "never-registered", retries=rq.Retry(max=1))
    assert workers.reap_dead_jobs(connection) == [(queue.name, job.id)]
    assert job.get_status(refresh=True) == "queued" and job.id in queue.get_job_ids()
    assert queue.started_job_registry.get_job_ids(cleanup=False) == []


def test_a_job_whose_worker_is_alive_is_left_running(queue_of_its_own):
    queue, connection = queue_of_its_own
    live = f"live-{uuid.uuid4().hex[:8]}"
    register_worker(queue, connection, live, os.getpid())
    job = running_job(queue, connection, live)
    assert workers.reap_dead_jobs(connection) == []
    assert job.get_status(refresh=True) == "started"
    assert queue.started_job_registry.get_job_ids(cleanup=False) == [job.id]


def test_a_registry_entry_whose_job_has_expired_is_dropped(queue_of_its_own):
    queue, connection = queue_of_its_own
    connection.zadd(queue.started_job_registry.key, {"0123456789abcdef:abcd": ALL_FOREVER})
    assert workers.reap_dead_jobs(connection) == [(queue.name, "0123456789abcdef")]
    assert queue.started_job_registry.get_job_ids(cleanup=False) == []


def test_the_update_guard_does_not_list_a_job_whose_worker_is_gone(redis_server, monkeypatch):
    unique = f"planetgen-api-test-{uuid.uuid4().hex[:10]}"
    monkeypatch.setattr(redisqueue, "QUEUES", (unique,))
    connection = redisqueue.connect(redis_server)
    queue = redisqueue.queue(unique, connection)
    try:
        job = running_job(queue, connection, "never-registered")
        assert job.id and updateguard.QueuedJobSource(redis_server).items() == []
        assert job.get_status(refresh=True) == "failed"
    finally:
        queue.empty()
        connection.delete(queue.started_job_registry.key, queue.failed_job_registry.key)
        connection.srem("rq:queues", f"rq:queue:{unique}")


def test_a_worker_clears_the_jobs_of_a_killed_worker_when_it_starts(queue_of_its_own, redis_server):
    queue, connection = queue_of_its_own
    job = running_job(queue, connection, "never-registered")
    done = subprocess.run([sys.executable, "-m", "planetgen.cli.worker", "--burst", "--url", redis_server, queue.name],
                          capture_output=True, text=True, timeout=300, cwd=REPO,
                          env={**os.environ, "PYTHONPATH": os.path.join(REPO, "src")})
    assert done.returncode == 0, done.stderr
    assert job.get_status(refresh=True) == "failed"


# --- planetgen.cli.reset_jobs ---------------------------------------------------------------------------------

def database_of_its_own(redis_server):
    """The test server's Redis URL with the database number 14, which no other test uses."""
    return redis_server.rsplit("/", 1)[0] + "/14" if redis_server.count("/") > 2 else redis_server + "/14"


def test_reset_jobs_clears_the_queue_the_jobs_and_the_lease_and_keeps_the_rest(redis_server, mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.close()
    url = database_of_its_own(redis_server)
    connection = redisqueue.connect(url)
    connection.flushdb()
    queue = redisqueue.queue(redisqueue.INTERACTIVE_QUEUE, connection)
    waiting = queue.enqueue(math.pow, 2, 2, job_timeout=600)
    running_job(queue, connection, "never-registered")
    connection.set(workers.WANTED_KEY, 3)
    connection.set("planetgen:page:kept", "1")
    connection.set(workers.STARTING_PREFIX + redisqueue.INTERACTIVE_QUEUE, 2)
    pause.pause("boss", connection)
    stale_run = runrecord.start_run(connection, "Stale run", [runrecord.Stage("mass", "Mass", units=3)])
    with pytest.raises(runrecord.Busy):
        runrecord.start_run(connection, "Blocked", [runrecord.Stage("mass", "Mass")])
    node = work.open_node("plan", "Stuck run", control_config=mysql_config)
    control = store.get_control_connection(mysql_config)
    try:
        control.execute("INSERT IGNORE INTO work_lease (id) VALUES (1)")
        control.execute("UPDATE work_lease SET holder = 'host:1:x', job_id = ?, heartbeat_at = NOW(6) "
                        "WHERE id = 1", (node.id,))
        control.commit()
        lines = []
        stuck = FakeWork("Stuck run")

        assert reset_jobs.reset_jobs(mysql_config, redis_url=url, dry_run=True, sources=[stuck], write=lines.append,
                                     control=mysql_config)
        assert stuck.cancelled == [] and connection.exists(waiting.key) == 1
        assert any("Redis: " in line and "key(s)" in line for line in lines)
        assert any("1 job(s)" in line and "held by host:1:x" in line for line in lines)

        questions = []
        assert not reset_jobs.reset_jobs(mysql_config, redis_url=url, sources=[stuck], write=lines.append,
                                         answer=lambda text: questions.append(text) or "n", control=mysql_config)
        assert len(questions) == 1 and stuck.cancelled == [] and connection.exists(waiting.key) == 1

        lines.clear()
        assert reset_jobs.reset_jobs(mysql_config, redis_url=url, assume_yes=True, sources=[stuck], write=lines.append,
                                     control=mysql_config)
        assert stuck.cancelled == ["job-0"]
        assert list(connection.scan_iter(match="rq:*")) == []
        assert list(connection.scan_iter(match="planetgen:run:*")) == [] and runrecord.active_run(connection) is None
        assert runrecord.read_run(connection, stale_run) is None
        assert list(connection.scan_iter(match=workers.STARTING_PREFIX + "*")) == []
        assert connection.get(workers.WANTED_KEY) == b"3" and connection.get("planetgen:page:kept") == b"1"
        assert control.execute("SELECT COUNT(*) AS n FROM work_jobs").fetchone()["n"] == 0
        lease = control.execute("SELECT holder, job_id FROM work_lease WHERE id = 1").fetchone()
        assert (lease["holder"], lease["job_id"]) == (None, None)
        assert not pause.is_paused(connection)          # the pause flag goes with the rest
        text = "\n".join(lines)
        assert "Redis: removed" in text and "Control database: removed 1 job(s)" in text and "Jobs and queue reset" in text

        lines.clear()
        assert reset_jobs.reset_jobs(mysql_config, redis_url=url, assume_yes=True, sources=[], write=lines.append,
                                     control=mysql_config)
        assert lines[-1] == "Nothing to reset."
        assert runrecord.start_run(connection, "Fresh", [runrecord.Stage("mass", "Mass")])      # no Busy
    finally:
        work.close_node(node, "cancelled")      # the run this process opened, whose rows the reset removed
        control.close()
        connection.flushdb()


def test_reset_jobs_kills_a_worker_process_of_this_host(redis_server):
    url = database_of_its_own(redis_server)
    connection = redisqueue.connect(url)
    connection.flushdb()
    sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(600)"], start_new_session=True)
    try:
        queue = redisqueue.queue(redisqueue.INTERACTIVE_QUEUE, connection)
        register_worker(queue, connection, "sleeper-worker", sleeper.pid)
        lines = []
        assert reset_jobs.reset_jobs(None, redis_url=url, assume_yes=True, sources=[], write=lines.append)
        assert sleeper.wait(timeout=30) == -9
        assert any("Killing worker sleeper-worker" in line for line in lines)
    finally:
        sleeper.kill()
        sleeper.wait()
        connection.flushdb()


def test_reset_jobs_without_a_terminal_or_the_word_yes_does_nothing():
    done = subprocess.run([sys.executable, "-m", "planetgen.cli.reset_jobs"], capture_output=True, text=True,
                          timeout=120, stdin=subprocess.DEVNULL, cwd=REPO,
                          env={**os.environ, "PYTHONPATH": os.path.join(REPO, "src")})
    assert done.returncode == 1 and "--yes" in done.stderr


# --- update.sh --reset-jobs -----------------------------------------------------------------------------------

LOGGING_PYTHON = r"""
        printf '#!/bin/sh\necho "python $*" >> "$LOG"\nexit "${PY_STATUS:-0}"\n' > "$VENV/bin/python"
"""


def test_update_sh_resets_the_jobs_once_before_it_takes_the_guard(tmp_path):
    log = tmp_path / "calls"
    done = run_shell(LOGGING_PYTHON + '''
        RESET_JOBS=1
        ensure_update_guard; ensure_update_guard
    ''', status_out=str(log))
    assert done.returncode == 0, done.stderr
    calls = log.read_text().splitlines()
    assert calls[0] == "python -m planetgen.cli.reset_jobs --yes"
    assert calls[1].startswith("update_guard begin") and len(calls) == 2


def test_update_sh_without_the_flag_resets_nothing(tmp_path):
    log = tmp_path / "calls"
    done = run_shell(LOGGING_PYTHON + "ensure_update_guard", status_out=str(log))
    assert done.returncode == 0, done.stderr
    assert [call.split()[0] for call in log.read_text().splitlines()] == ["update_guard"]


def test_a_failed_reset_stops_the_update_before_the_guard(tmp_path):
    log = tmp_path / "calls"
    done = run_shell(LOGGING_PYTHON + '''
        RESET_JOBS=1
        PY_STATUS=1 ensure_update_guard
        echo unreachable
    ''', status_out=str(log))
    assert done.returncode == 1 and "unreachable" not in done.stdout
    assert "could not be reset" in done.stderr
    assert log.read_text().splitlines() == ["python -m planetgen.cli.reset_jobs --yes"]


def test_update_sh_knows_the_flag_and_says_so_in_its_help_and_the_deployment_notes():
    with open(os.path.join(REPO, "update.sh"), encoding="utf-8") as handle:
        text = handle.read()
    assert "--reset-jobs) RESET_JOBS=1" in text and "[--cancel-running] [--reset-jobs]" in text
    with open(os.path.join(REPO, "docs", "deployment", "README.md"), encoding="utf-8") as handle:
        assert "--reset-jobs" in handle.read()
