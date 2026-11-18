"""OPS.213: the supervisor with real processes: a real supervisor process, real Redis, and a real kill (SIGKILL) of the supervisor and of the process that started a run.

The stage 0 skeleton and the simulator (`tests/supervisor_sim.py`, OPS.212) cover every order of faults on a fake clock; these show the same promises hold with the operating system's own
processes, the host lock and the sign of life that runs out in real seconds. They start no workers and generate no galaxy, so they run on a small box; the kill of a unit's process in the
conclusion, and the comparison of a galaxy with an undisturbed run, are the next test (a galaxy run needs the gate machine).
"""

import os
import signal
import subprocess
import sys
import threading
import time

import pytest

from planetgen.queue import redisqueue, runline, runrecord, supervisor

SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


@pytest.fixture
def host(redis_server, tmp_path, monkeypatch):
    """A state folder of this test's own (so its supervisor's lock is not the machine's), and the processes it starts, which end with the test."""
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))
    monkeypatch.setattr(runrecord, "PREFIX", "planetgen")        # the processes started here know no other prefix than the real one: the suite's private prefix is not theirs
    env = {**os.environ, "PLANETGEN_STATE_DIR": str(tmp_path), "PYTHONPATH": os.pathsep.join([SRC, *sys.path])}
    processes = []

    class Host:
        url = redis_server
        conn = redisqueue.connect(redis_server)

        @staticmethod
        def start(argv):
            process = subprocess.Popen(argv, env=env, cwd=SRC)
            processes.append(process)
            token = supervisor._claim_token()

            def release():          # what `supervisor._spawn` does through `redisqueue.start_detached(on_exit=...)` (OPS.251): the claim to be starting ends with the process it was made for
                process.wait()
                supervisor._release_claim(Host.conn, token)

            threading.Thread(target=release, daemon=True).start()

        @staticmethod
        def ask():
            supervisor.ensure_started(Host.conn, redis_server, start=Host.start)

        @staticmethod
        def supervisor_pid():
            fields = Host.conn.hget(supervisor.alive_key(), "pid")
            return int(runrecord._s(fields)) if fields else None

    for key in Host.conn.scan_iter(match="planetgen:*"):
        Host.conn.delete(key)
    yield Host
    pid = Host.supervisor_pid()
    if pid:
        try:
            os.kill(pid, signal.SIGKILL)         # one a starter of the test started, detached: it is not in `processes`
        except ProcessLookupError:
            pass
    for process in processes:
        if process.poll() is None:
            process.kill()
        process.wait()
    for key in Host.conn.scan_iter(match="planetgen:*"):
        Host.conn.delete(key)


def _start(title):
    return {"title": title, "stages": [{"key": "mass", "label": "Mass"}]}


def _wait_until(check, seconds, what):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        found = check()
        if found:
            return found
        time.sleep(0.2)
    raise AssertionError(f"{what} did not happen within {seconds} seconds")


def test_a_real_supervisor_process_starts_the_run_that_asks_and_a_second_one_refuses_the_lock(host):
    run_id = runline.start_when_free(host.conn, _start("First"), 60, host.ask)
    assert runrecord.active_run(host.conn) == run_id
    pid = _wait_until(host.supervisor_pid, 10, "the supervisor's sign of life")
    assert subprocess.run([sys.executable, "-m", "planetgen.cli.supervisor", "--url", host.url], cwd=SRC, env={**os.environ, "PYTHONPATH": SRC},
                          capture_output=True, timeout=60, check=False).returncode == 0
    assert host.supervisor_pid() == pid, "a second supervisor leaves the first the host"
    runrecord.finish_run(host.conn, run_id, "done")


def test_a_run_waiting_for_the_machine_is_started_by_the_next_supervisor_when_the_first_is_killed(host):
    holder = runrecord.start_run(host.conn, "Holder", [runrecord.Stage("mass", "Mass")])
    host.ask()
    first = _wait_until(host.supervisor_pid, 10, "the first supervisor's sign of life")
    results = {}

    def waiter():
        results["run_id"] = runline.start_when_free(redisqueue.connect(host.url), _start("Waits"), 120, host.ask)

    thread = threading.Thread(target=waiter, daemon=True)
    thread.start()
    _wait_until(lambda: runline.of_type(host.conn, runline.RUN), 10, "the run's place in the line")
    os.kill(first, signal.SIGKILL)                                      # the supervisor dies with a run waiting, without a word
    runrecord.finish_run(host.conn, holder, "done")                     # the machine is free, and nobody is there to start the waiting run
    thread.join(timeout=90)
    assert not thread.is_alive() and results.get("run_id"), "the waiting run was started by a supervisor that took the place of the dead one"
    assert host.supervisor_pid() not in (None, first)
    assert runrecord.active_run(host.conn) == results["run_id"]
    runrecord.finish_run(host.conn, results["run_id"], "done")


def test_a_run_whose_starter_is_killed_is_carried_on_and_the_supervisor_says_so_once(host):
    code = ("import sys, time\nfrom planetgen.queue import runline, redisqueue, supervisor\nconn = redisqueue.connect(sys.argv[1])\n"
            "rid = runline.start_when_free(conn, {'title': 'Starter', 'stages': [{'key': 'mass', 'label': 'Mass'}]}, 60,\n"
            "                              lambda: supervisor.ensure_started(conn, sys.argv[1]))\nprint(rid, flush=True)\ntime.sleep(600)\n")
    starter = subprocess.Popen([sys.executable, "-c", code, host.url], stdout=subprocess.PIPE, text=True, cwd=SRC,
                               env={**os.environ, "PYTHONPATH": SRC, "PLANETGEN_STATE_DIR": os.environ["PLANETGEN_STATE_DIR"]})
    try:
        run_id = starter.stdout.readline().strip()
        assert run_id and runrecord.active_run(host.conn) == run_id
        os.kill(starter.pid, signal.SIGKILL)                            # the process that started the run is gone, and nothing follows the run any more
        starter.wait()
        host.conn.delete(runrecord.alive_key(run_id))
        said = _wait_until(lambda: [fields for _id, fields in runrecord.read_events(host.conn, run_id) if fields.get("kind") == "run.starter_stopped"], 30,
                           "the supervisor's line about the starter")
        time.sleep(6)
        said = [fields for _id, fields in runrecord.read_events(host.conn, run_id) if fields.get("kind") == "run.starter_stopped"]
        assert len(said) == 1 and "carries the run on" in said[0]["message"]
        assert runrecord.read_run(host.conn, run_id)["state"] == "running" and not runrecord.holder_is_dead(host.conn, run_id)
    finally:
        if starter.poll() is None:
            starter.kill()
        host.conn.hset(runrecord.run_key(run_id), "state", "done") if run_id else None
