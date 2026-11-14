# tests/test_queue_command_env.py

"""TEST.717: a test that runs a queue command in a process of its own gets a Redis database with no worker on it (`queue_command_env`, tests/conftest.py)."""

import ast
import os
import signal
import subprocess
import sys
import time
from urllib.parse import urlsplit

import pytest

from tests import stray_workers

TESTS = os.path.dirname(os.path.abspath(__file__))

QUEUE_COMMANDS = ("check-math", "status", "plan")
"""tuple: The commands of `planetgen.cli.generate` and `planetgen.cli.datasets` that run through the queue as check runs or plans (`wrapper.run`)."""


@pytest.fixture
def decoy_worker():
    """A worker process on this test worker's own database, like the one an earlier test leaves behind."""
    url = os.environ.get("PLANETGEN_REDIS_URL")
    path = urlsplit(url or "").path.strip("/")
    if not path.isdigit() or int(path) not in stray_workers.TEST_DATABASES:
        pytest.skip("no test Redis database for this worker (PLANETGEN_TEST_REDIS_URL is not set)")
    process = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)", "-m", stray_workers.WORKER_COMMAND, "--url", url],
                               start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    deadline = time.monotonic() + 10
    while process.pid not in [w["pid"] for w in stray_workers.strays(databases=[int(path)])] and time.monotonic() < deadline:
        time.sleep(0.1)
    yield process, int(path)
    if process.poll() is None:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()


def test_a_worker_an_earlier_test_left_on_the_database_is_ended_before_a_command_runs(decoy_worker, queue_command_env):
    process, database = decoy_worker
    assert process.wait(timeout=20) == -signal.SIGTERM
    assert not stray_workers.strays(databases=[database])
    assert queue_command_env["PLANETGEN_REDIS_URL"] == os.environ["PLANETGEN_REDIS_URL"]


def _queue_command_tests():
    """`(file, function, args)` of every test that runs `-m planetgen.cli.generate` or `-m planetgen.cli.datasets` with a queue command in a process."""
    found = []
    for name in sorted(os.listdir(TESTS)):
        if not (name.startswith("test_") and name.endswith(".py")) or name == os.path.basename(__file__):
            continue
        with open(os.path.join(TESTS, name), encoding="utf-8") as handle:
            source = handle.read()
        if "subprocess" not in source or "planetgen.cli." not in source:
            continue
        for node in ast.walk(ast.parse(source)):
            if not isinstance(node, ast.FunctionDef) or not node.name.startswith("test_"):
                continue
            text = ast.get_source_segment(source, node) or ""
            runs_cli = "subprocess." in text and ('"planetgen.cli.generate"' in text or '"planetgen.cli.datasets"' in text)
            if runs_cli and any(f'"{command}"' in text for command in QUEUE_COMMANDS):
                found.append((name, node.name, [arg.arg for arg in node.args.args]))
    return found


def test_the_guard_finds_the_tests_that_run_a_queue_command_in_a_process():
    names = {(file, function) for file, function, _ in _queue_command_tests()}
    assert ("test_math_gate.py", "test_check_math_runs_from_the_command_line") in names
    assert ("test_datasets_gate.py", "test_the_status_command_the_generate_page_runs_exits_1_while_data_is_missing") in names


def test_every_test_that_runs_a_queue_command_in_a_process_takes_the_queue_command_env():
    """A command in a process hands its work to any worker serving the Redis database, in that worker's environment (the machine's data folder, say)."""
    missing = [f"{file}::{function}" for file, function, args in _queue_command_tests() if "queue_command_env" not in args]
    assert not missing, f"these tests run a queue command in a process without the `queue_command_env` fixture: {missing}"
