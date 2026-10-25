# tests/test_stray_workers.py

"""Finding and ending the queue workers a test run left behind (TEST.393, `tests/stray_workers.py`)."""

import os
import signal
import subprocess
import sys
import time

from tests import stray_workers

LISTING = """\
  100   100     5400 /usr/bin/python3 -m planetgen.cli.worker --url redis://127.0.0.1:6379/3 planetgen-bulk
  101   100     5399 /usr/bin/python3 -m planetgen.cli.worker --url redis://127.0.0.1:6379/3 planetgen-bulk
  200   200       12 /usr/bin/python3 -m planetgen.cli.worker --url redis://127.0.0.1:6379/0 planetgen-bulk
  300   300     9000 /usr/bin/python3 -m planetgen.cli.worker --url redis://127.0.0.1:6379 planetgen-bulk
  400   400      800 /usr/bin/python3 -m planetgen.cli.worker --burst planetgen-bulk
  500   500      700 python -m pytest tests/test_workers.py -k planetgen.cli.worker --url redis://127.0.0.1:6379/4
  600   600       30 /usr/bin/python3 -m planetgen.cli.worker --url redis://127.0.0.1:6379/15 planetgen-interactive
  700   700       30 /usr/bin/mysqld --basedir=/usr
garbage line
"""


def test_a_worker_is_read_with_its_age_and_the_redis_database_it_names():
    found = {worker["pid"]: worker for worker in stray_workers.parse(LISTING)}
    assert sorted(found) == [100, 200, 300, 400, 600]                 # the horse (101), pytest, mysqld and noise are not
    assert found[100] == {"pid": 100, "seconds": 5400, "database": 3}
    assert found[200]["database"] == 0 and found[300]["database"] == 0     # a URL with no database is database 0
    assert found[400]["database"] is None


def test_only_workers_on_a_private_test_database_are_strays():
    assert [worker["pid"] for worker in stray_workers.strays(LISTING)] == [100, 600]
    assert [worker["pid"] for worker in stray_workers.strays(LISTING, databases=[3])] == [100]
    assert [worker["pid"] for worker in stray_workers.strays(LISTING, older_than=60)] == [100]


def test_reap_ends_the_group_and_kills_what_stays():
    sent = []
    alive = {100}
    found = stray_workers.strays(LISTING)
    count = stray_workers.reap(found, killpg=lambda pid, sig: sent.append((pid, sig)), sleep=lambda seconds: None,
                               exists=lambda pid: pid in alive)
    assert count == 2
    assert sent == [(100, signal.SIGTERM), (600, signal.SIGTERM), (100, signal.SIGKILL)]


def test_a_worker_that_is_already_gone_is_not_counted():
    def gone(pid, sig):
        raise ProcessLookupError

    assert stray_workers.reap(stray_workers.strays(LISTING), killpg=gone, sleep=lambda seconds: None) == 0


def test_a_real_worker_process_on_a_test_database_is_found_and_ended():
    code = "import time; time.sleep(60)"
    process = subprocess.Popen([sys.executable, "-c", code, "-m", "planetgen.cli.worker", "--url",
                                "redis://127.0.0.1:6379/14"], start_new_session=True,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        deadline = time.monotonic() + 10
        while process.pid not in [w["pid"] for w in stray_workers.strays(databases=[14])] and time.monotonic() < deadline:
            time.sleep(0.1)
        assert stray_workers.reap_database("redis://127.0.0.1:6379/14") >= 1
        assert process.wait(timeout=20) == -signal.SIGTERM
    finally:
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()


def test_the_machines_own_database_is_never_reaped():
    assert stray_workers.reap_database("redis://127.0.0.1:6379/0") == 0
    assert stray_workers.reap_database("redis://127.0.0.1:6379") == 0


def test_the_runner_reports_what_it_ended(monkeypatch):
    import importlib.util
    import pathlib

    path = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "test.py"
    spec = importlib.util.spec_from_file_location("runner_for_stray_test", path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    monkeypatch.setattr(stray_workers, "reap_all", lambda: [{"pid": 1, "seconds": 3700, "database": 4},
                                                            {"pid": 2, "seconds": 20, "database": 5}])
    sentence = runner.reap_strays()
    assert sentence.startswith("2 queue worker process(es) a test started and did not stop")
    assert "1h" in sentence or "61" in sentence                     # the oldest, 3,700 s
    monkeypatch.setattr(stray_workers, "reap_all", lambda: [])
    assert runner.reap_strays() == ""
