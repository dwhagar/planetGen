# src/tests/test_stopall.py

"""
Stop everything (PERF.186): the hold keeps anything from starting, `stop_all` ends the processes and clears the queues,
Cancel ends a job whose runner is gone, and a step that ignores the request to stop is killed.
"""

import os
import signal
import subprocess
import sys
import threading
import time

import pytest

from planetgen.admin import updateguard
from planetgen.cli import queue as queue_cli
from planetgen.queue import pause, redisqueue, stopall
from planetgen.web import job_runner, jobs
from planetgen.web.lib import tileprefetch
from tests import page_jobs

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
    """A process whose command line names a module `stopall` ends (`planetgen.cli.warm`), as a real command's does (OPS.235: there is no worker module any more)."""
    return subprocess.Popen([PYTHON, "-c", code, "planetgen.cli.warm"], start_new_session=True)


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


def test_the_tile_prefetch_queues_nothing_while_stopped(monkeypatch):
    from planetgen.queue import lightjobs
    added = []
    monkeypatch.setattr(lightjobs, "add_job", lambda *args, **kwargs: added.append(args) or "0" * 16)
    assert tileprefetch._add_chunks(None, "db", "gen", ["a"]) == 1
    assert len(added) == 1
    stopall._write_hold("boss")
    assert tileprefetch._add_chunks(None, "db", "gen", ["b"]) == 0
    assert len(added) == 1


def test_find_targets_takes_workers_and_everything_they_started_and_spares_the_rest():
    table = {
        10: (1, "/venv/bin/python -m planetgen.cli.warm --url redis://x planetgen-bulk"),
        11: (10, "/venv/bin/python -m planetgen.cli.warm --url redis://x planetgen-bulk"),    # its work horse
        12: (11, "sleep 100"),                                                                    # whatever it started
        20: (1, "/venv/bin/python -m planetgen.cli.generate sectors --radius-pc 12"),
        21: (20, "/venv/bin/python -c from multiprocessing.spawn import spawn_main"),
        30: (1, "/venv/bin/python -c from multiprocessing.spawn import spawn_main"),             # an orphaned pool child
        40: (1, "/venv/bin/python -m planetgen.cli.migrate"),
        41: (1, "apache2 -k start"),
        42: (1, "/usr/sbin/mysqld"),
        50: (1, "/venv/bin/python -m planetgen.cli.warm planetgen-bulk"),                      # protected
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
    from planetgen.queue import lightjobs
    monkeypatch.setattr(lightjobs, "ask_supervisor", lambda connection=None, url=None: None)
    connection = redisqueue.connect(redis_server)
    for _ in range(2):
        lightjobs.add_job(connection, "math:pow", [2, 3])      # OPS.235: the jobs that wait are the site's light jobs
    worker = _spawn_worker("import time; time.sleep(60)")
    time.sleep(0.5)
    try:
        done = stopall.stop_all("boss", connection=connection)
        assert done["already"] is False
        assert done["queued"] >= 2
        assert done["killed"] >= 1
        assert lightjobs.live(connection) == []
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
    page_jobs.forget_all_jobs()
    job_id = jobs.start_job("galaxy", "Stuck", [{"label": "Runs", "argv": [PYTHON, "-c", step]}], spawn=False)
    return job_id, os.path.join(jobs.jobs_dir(), job_id)


def test_cancel_of_a_job_whose_runner_died_ends_its_step_and_frees_the_lock(tmp_path, monkeypatch):
    job_id, path = _job(tmp_path, monkeypatch, "pass")
    step = subprocess.Popen([PYTHON, "-c", "import time; time.sleep(60)"], start_new_session=True)
    try:
        page_jobs.running(job_id, pid=page_jobs.DEAD_PID)           # the runner of the step is gone, the step runs on
        with open(os.path.join(path, job_runner.STEP_PID_NAME), "w", encoding="utf-8") as handle:
            handle.write(str(step.pid))
        assert jobs.get_job(job_id)["status"] == "interrupted"      # what the page shows: the runner is gone
        assert jobs.cancel_job(job_id, who="boss") is True
        assert step.wait(timeout=10) is not None                   # its step was still running; now it is not
        assert jobs.get_job(job_id)["status"] == "canceled"
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
    monkeypatch.setattr(job_runner, "ROW_CHECK_SECONDS", 0.2)
    runner = threading.Thread(target=page_jobs.run_to_end, args=(job_id,), daemon=True)
    runner.start()
    deadline = time.time() + 30
    while time.time() < deadline and "up" not in (jobs.read_log(job_id) or ("",))[0]:
        time.sleep(0.1)
    started = time.time()
    assert jobs.cancel_job(job_id, who="boss") is True
    runner.join(timeout=60)
    assert not runner.is_alive() and time.time() - started < 60
    assert jobs.get_job(job_id)["status"] == "canceled"
