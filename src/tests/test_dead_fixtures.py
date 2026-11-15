"""TEST.190: a module-scoped fixture whose setup hung is not set up again, and so not waited for, by every later module."""
import os
import subprocess
import sys
import textwrap
import time
from pathlib import Path

SRC = str(Path(__file__).resolve().parents[1])

SHARED = '''
import time

import pytest


@pytest.fixture(scope="module")
def server():
    time.sleep(100000)
    yield "up"


@pytest.fixture
def per_test():
    time.sleep(100000)
    yield "up"
'''

MODULE = '''
from helper import server  # noqa: F401


def test_one(server):
    assert server == "up"


def test_two(server):
    assert server == "up"
'''

PER_TEST = '''
from helper import per_test  # noqa: F401


def test_one(per_test):
    assert per_test == "up"
'''


def _run(tmp_path, files, extra=(), env=None):
    (tmp_path / "pytest.ini").write_text("[pytest]\ntimeout = 3\ntimeout_method = signal\n")
    for name, source in files.items():
        (tmp_path / name).write_text(textwrap.dedent(source))
    started = time.monotonic()
    done = subprocess.run([sys.executable, "-m", "pytest", "-p", "tests.dead_fixtures", "-p", "no:cacheprovider", "-q", *extra], cwd=tmp_path,
                          env={**os.environ, "PYTHONPATH": SRC + os.pathsep + str(tmp_path), **(env or {})}, capture_output=True, text=True, timeout=120)
    return done.stdout + done.stderr, time.monotonic() - started


def test_a_server_fixture_that_hung_is_not_waited_for_by_each_later_module(tmp_path):
    output, seconds = _run(tmp_path, {"helper.py": SHARED, "test_a.py": MODULE, "test_b.py": MODULE, "test_c.py": MODULE})
    assert "6 errors" in output
    assert "TEST.190: fixture 'server' hung in the setup for test_a.py" in output, output
    assert seconds < 10, f"every module waited its own timeout: {seconds:.0f} s"


def test_a_function_fixture_that_timed_out_is_set_up_again_for_the_next_test(tmp_path):
    output, _ = _run(tmp_path, {"helper.py": SHARED, "test_a.py": PER_TEST, "test_b.py": PER_TEST})
    assert "TEST.190" not in output and "2 errors" in output


def test_the_switch_sets_every_module_up_as_pytest_does(tmp_path):
    output, _ = _run(tmp_path, {"helper.py": SHARED, "test_a.py": MODULE, "test_b.py": MODULE}, env={"PLANETGEN_TEST_RETRY_DEAD_FIXTURES": "1"})
    assert "TEST.190" not in output


def test_the_plugin_is_loaded_by_the_repo_settings():
    with open(os.path.join(os.path.dirname(SRC), "pytest.ini"), encoding="utf-8") as handle:
        assert "-p tests.dead_fixtures" in handle.read()
