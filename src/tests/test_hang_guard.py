"""TEST.621: a test that hangs past its limit is stopped and named, and the run goes on."""
import os
import subprocess
import sys
import textwrap
import time
from pathlib import Path

from tests import hang_guard

SRC = str(Path(__file__).resolve().parents[1])

POOL_TEST = '''
import concurrent.futures
import time


def _sleep(seconds):
    time.sleep(seconds)


def test_pool_that_never_ends():
    with concurrent.futures.ProcessPoolExecutor(2) as pool:
        for future in [pool.submit(_sleep, 100000) for _ in range(2)]:
            future.result()


def test_after_it():
    assert True
'''

STUCK_TEST = '''
import signal
import time


def test_stuck_on_nothing_a_child_can_end():
    signal.signal(signal.SIGALRM, signal.SIG_IGN)
    while True:
        time.sleep(1)


def test_after_it():
    assert True
'''


def _run(tmp_path, source, extra=()):
    (tmp_path / "pytest.ini").write_text("[pytest]\ntimeout = 2\ntimeout_method = signal\nhang_guard_grace = 2\n")
    (tmp_path / "test_case.py").write_text(textwrap.dedent(source))
    started = time.monotonic()
    done = subprocess.run([sys.executable, "-m", "pytest", "-p", "tests.hang_guard", "-p", "no:cacheprovider", "-q", *extra],
                          cwd=tmp_path, env={**os.environ, "PYTHONPATH": SRC}, capture_output=True, text=True, timeout=120)
    return done, time.monotonic() - started


def test_a_pool_that_never_finishes_is_killed_and_the_next_test_runs(tmp_path):
    done, seconds = _run(tmp_path, POOL_TEST)
    output = done.stdout + done.stderr
    assert "HUNG TEST test_case.py::test_pool_that_never_ends" in output
    assert "killed 2 child process(es)" in output
    assert "1 failed, 1 passed" in output
    assert seconds < 60


def test_a_test_that_waits_on_no_child_ends_its_worker_under_xdist_and_the_run_goes_on(tmp_path):
    done, seconds = _run(tmp_path, STUCK_TEST, ("-n", "2"))
    output = done.stdout + done.stderr
    assert "crashed while running 'test_case.py::test_stuck_on_nothing_a_child_can_end'" in output
    assert "1 passed" in output and "1 failed" in output
    assert seconds < 90


def test_a_test_with_no_limit_is_not_watched(tmp_path):
    (tmp_path / "pytest.ini").write_text("[pytest]\ntimeout = 0\nhang_guard_grace = 0\n")
    (tmp_path / "test_case.py").write_text("import time\n\n\ndef test_slow():\n    time.sleep(1)\n")
    done = subprocess.run([sys.executable, "-m", "pytest", "-p", "tests.hang_guard", "-p", "no:cacheprovider", "-q"], cwd=tmp_path,
                          env={**os.environ, "PYTHONPATH": SRC}, capture_output=True, text=True, timeout=60)
    assert "1 passed" in done.stdout and "HUNG" not in done.stdout + done.stderr


def test_children_of_finds_a_grandchild_and_kill_children_ends_it():
    child = subprocess.Popen([sys.executable, "-c", "import subprocess, sys, time; "
                              "subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)']); time.sleep(60)"])
    try:
        deadline = time.monotonic() + 10
        while len(hang_guard.children_of(os.getpid())) < 2 and time.monotonic() < deadline:
            time.sleep(0.1)
        assert child.pid in hang_guard.children_of(os.getpid())
        assert len(hang_guard.children_of(os.getpid())) == 2
        assert hang_guard.kill_children() == 2
        assert child.wait(timeout=10) == -9
    finally:
        child.kill()
        child.wait()
