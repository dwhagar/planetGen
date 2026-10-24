# src/tests/test_stopall.py

"""
Stop everything (PERF.186): the hold keeps anything from starting, `stop_all` ends the processes and clears the queues,
Cancel ends a job whose runner is gone, and a step that ignores the request to stop is killed.
"""

import json
import os
import signal
import subprocess
import sys
import time

import pytest

from planetgen.admin import updateguard
from planetgen.cli import queue as queue_cli, worker as worker_cli
from planetgen.queue import pause, redisqueue, stopall, workers
from planetgen.web import job_runner, jobs
from planetgen.web.lib import tileprefetch

PYTHON = sys.executable


@pytest.fixture(autouse=True)
def no_hold():
    """Every test starts and ends with no hold on."""
    try:
        os.remove(stopall.path())
    except OSError:
        pass
    yield
    try:
        os.remove(stopall.path())
    except OSError:
        pass


def _only_my_children(monkeypatch):
    """`stopall` sees only the processes this test started, so the host's own workers are left alone."""
    real = stopall.process_table

    def table():
        full = real()
        mine = {os.getpid()}
        grew = True
        while grew:
            grew = False
            for pid, (ppid, _args) in full.items():
                if ppid in mine and pid not in mine:
                    mine.add(pid)
                    grew = True
        mine.discard(os.getpid())
        return {pid: row for pid, row in full.items() if pid in mine}
    monkeypatch.setattr(stopall, "process_table", table)


def _spawn_worker(code):
    """A process whose command line names the worker module, as a real worker's does."""
    return subprocess.Popen([PYTHON, "-c", code, "planetgen.cli.worker"], start_new_session=True)


def _gone(pid, seconds=10):
    deadline = time.time() + seconds
    while time.time() < deadline:
        if not stopall._alive(pid):
            return True
        time.sleep(0.1)
    return False


def test_the_hold_is_a_file_that_names_who_and_when():
    assert stopall.state()["stopped"] is False
    stopall.refuse_while_stopped("A worker")        # nothing is raised without the hold
    assert stopall._write_hold("boss") is True
    assert stopall._write_hold("someone else") is False     # the first hold is kept
    held = stopall.state()
    assert held["stopped"] and held["by"] == "boss" and held["at"]
    with pytest.raises(stopall.Stopped, match="stopped by boss"):
        stopall.refuse_while_stopped("A Generate job")


def test_a_damaged_hold_file_still_holds():
    os.makedirs(os.path.dirname(stopall.path()), exist_ok=True)
    with open(stopall.path(), "w", encoding="utf-8") as handle:
        handle.write("{not json")
    assert stopall.is_stopped()


def test_nothing_starts_while_everything_is_stopped(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))
    stopall._write_hold("boss")
    with pytest.raises(updateguard.UpdateRunning, match="everything was stopped"):
        updateguard.refuse_while_updating("A run")      # the guard every start of a run, a job or the orbit update calls
    with pytest.raises(updateguard.UpdateRunning):
        jobs.start_job("galaxy", "Held", [{"label": "x", "argv": [PYTHON, "-c", "pass"]}], spawn=False)
    assert workers.ensure_workers(redisqueue.MAINTENANCE_QUEUE) == 0
    assert worker_cli.main(["--burst"]) == updateguard.EX_TEMPFAIL


def test_the_tile_prefetch_queues_nothing_while_stopped():
    class Queue:
        queued = 0

        def enqueue(self, *args, **kwargs):
            Queue.queued += 1
    tileprefetch._enqueue_chunks(Queue(), "db", ["a"], "prefetch")
    assert Queue.queued == 1
    stopall._write_hold("boss")
    tileprefetch._enqueue_chunks(Queue(), "db", ["b"], "prefetch")
    assert Queue.queued == 1


def test_find_targets_takes_workers_and_everything_they_started_and_spares_the_rest():
    table = {
        10: (1, "/venv/bin/python -m planetgen.cli.worker --url redis://x planetgen-bulk"),
        11: (10, "/venv/bin/python -m planetgen.cli.worker --url redis://x planetgen-bulk"),    # its work horse
        12: (11, "sleep 100"),                                                                    # whatever it started
        20: (1, "/venv/bin/python -m planetgen.cli.generate sectors --radius-pc 12"),
        21: (20, "/venv/bin/python -c from multiprocessing.spawn import spawn_main"),
        30: (1, "/venv/bin/python -c from multiprocessing.spawn import spawn_main"),             # an orphaned pool child
        40: (1, "/venv/bin/python -m planetgen.cli.migrate"),
        41: (1, "apache2 -k start"),
        42: (1, "/usr/sbin/mysqld"),
        50: (1, "/venv/bin/python -m planetgen.cli.worker planetgen-bulk"),                      # protected
    }
    assert sorted(stopall.find_targets(table, protected={50})) == [10, 11, 12, 20, 21, 30]
    order = stopall.find_targets(table, protected={50})
    assert order.index(10) < order.index(11) < order.index(12)      # parents before children


def test_kill_processes_ends_a_worker_its_children_and_one_that_ignores_sigterm(monkeypatch):
    _only_my_children(monkeypatch)
    polite = _spawn_worker("import subprocess, time; subprocess.Popen(['sleep', '60']); time.sleep(60)")
    stubborn = _spawn_worker("import signal, time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(60)")
    bystander = subprocess.Popen([PYTHON, "-c", "import time; time.sleep(60)"], start_new_session=True)
    time.sleep(1.0)
    try:
        done = stopall.kill_processes(term_seconds=1.0)
        assert done["killed"] >= 3          # the worker, its sleep and the stubborn one
        assert done["survivors"] == []
        assert polite.wait(timeout=10) is not None
        assert stubborn.wait(timeout=10) == -signal.SIGKILL
        assert bystander.poll() is None     # not a planetgen process: left alone
    finally:
        for proc in (polite, stubborn, bystander):
            if proc.poll() is None:
                proc.kill()
            proc.wait()


def test_stop_all_holds_cancels_empties_the_queues_and_kills(redis_server, monkeypatch, tmp_path):
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))
    _only_my_children(monkeypatch)
    connection = redisqueue.connect(redis_server)
    queue = redisqueue.queue(redisqueue.MAINTENANCE_QUEUE, connection)
    queue.empty()
    queue.enqueue("time.sleep", 1, job_timeout=60)
    queue.enqueue("time.sleep", 1, job_timeout=60)
    worker = _spawn_worker("import time; time.sleep(60)")
    time.sleep(0.5)
    try:
        done = stopall.stop_all("boss", connection=connection)
        assert done["already"] is False
        assert done["queued"] >= 2
        assert done["killed"] >= 1
        assert queue.count == 0
        assert stopall.state()["by"] == "boss"
        assert pause.is_paused(connection)
        assert _gone(worker.pid)
        # while it is on hold nothing starts, and resuming lifts both the hold and the queue pause
        with pytest.raises(stopall.Stopped):
            stopall.refuse_while_stopped("A worker")
        assert stopall.resume_all(connection) is True
        assert not stopall.is_stopped() and not pause.is_paused(connection)
        assert stopall.resume_all(connection) is False
    finally:
        worker.kill()
        worker.wait()


def test_stop_all_without_redis_still_holds_and_kills(monkeypatch, tmp_path):
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))
    monkeypatch.setattr(redisqueue, "connect", lambda *a, **k: (_ for _ in ()).throw(OSError("no Redis")))
    _only_my_children(monkeypatch)
    worker = _spawn_worker("import time; time.sleep(60)")
    time.sleep(0.5)
    try:
        done = stopall.stop_all("boss")
        assert stopall.is_stopped()
        assert done["killed"] >= 1 and _gone(worker.pid)
    finally:
        worker.kill()
        worker.wait()


def test_stop_all_without_a_hold_is_what_update_sh_kill_all_runs(monkeypatch, tmp_path):
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))
    monkeypatch.setattr(redisqueue, "connect", lambda *a, **k: (_ for _ in ()).throw(OSError("no Redis")))
    _only_my_children(monkeypatch)
    worker = _spawn_worker("import time; time.sleep(60)")
    time.sleep(0.5)
    try:
        done = stopall.stop_all("update.sh", hold=False)
        assert done["killed"] >= 1 and _gone(worker.pid)
        assert not stopall.is_stopped()
    finally:
        worker.kill()
        worker.wait()


def test_the_command_stops_resumes_and_reports(monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))
    monkeypatch.setattr(redisqueue, "connect", lambda *a, **k: (_ for _ in ()).throw(OSError("no Redis")))
    monkeypatch.setattr(stopall, "kill_processes", lambda *a, **k: {"killed": 3, "survivors": []})
    assert queue_cli.main(["stop-all"]) == 0
    assert "ended 3 process(es)" in capsys.readouterr().out
    assert queue_cli.main(["status"]) == 0
    assert "on hold since" in capsys.readouterr().out
    assert queue_cli.main(["resume"]) == 0
    assert "Resumed" in capsys.readouterr().out
    assert queue_cli.main(["stop-all", "--no-hold"]) == 0
    assert not stopall.is_stopped()
    monkeypatch.setattr(stopall, "kill_processes", lambda *a, **k: {"killed": 1, "survivors": [4242]})
    assert queue_cli.main(["stop-all"]) == 1        # something could not be ended: say so in the status
    assert "4242" in capsys.readouterr().out


def _job(tmp_path, monkeypatch, step):
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))
    job_id = jobs.start_job("galaxy", "Stuck", [{"label": "Runs", "argv": [PYTHON, "-c", step]}], spawn=False)
    return job_id, os.path.join(jobs.jobs_dir(), job_id)


def test_cancel_of_a_job_whose_runner_died_ends_its_step_and_frees_the_lock(tmp_path, monkeypatch):
    job_id, path = _job(tmp_path, monkeypatch, "pass")
    step = subprocess.Popen([PYTHON, "-c", "import time; time.sleep(60)"], start_new_session=True)
    try:
        with open(os.path.join(path, "state.json"), "w", encoding="utf-8") as handle:
            json.dump({"status": "running", "pid": 999999, "step": 1, "step_pid": step.pid, "started_at": time.time()},
                      handle)
        assert jobs.get_job(job_id)["status"] == "interrupted"      # what the page shows: the runner is gone
        assert jobs.cancel_job(job_id, who="boss") is True
        assert step.wait(timeout=10) is not None                   # its step was still running; now it is not
        assert jobs.get_job(job_id)["status"] == "cancelled"
        assert jobs.active_job() is None
        assert jobs.cancel_job(job_id) is False                    # nothing left to cancel
    finally:
        if step.poll() is None:
            step.kill()
            step.wait()


def test_a_step_that_ignores_the_request_to_stop_is_killed(tmp_path, monkeypatch):
    stubborn = "import signal, time; signal.signal(signal.SIGTERM, signal.SIG_IGN); print('up', flush=True); time.sleep(120)"
    job_id, path = _job(tmp_path, monkeypatch, stubborn)
    monkeypatch.setattr(job_runner, "CANCEL_KILL_SECONDS", 1)
    with open(os.path.join(path, "cancel"), "w", encoding="utf-8") as handle:
        handle.write("now")
    started = time.time()
    assert job_runner.run(path) == 1
    assert time.time() - started < 60
    assert jobs.get_job(job_id)["status"] == "cancelled"
