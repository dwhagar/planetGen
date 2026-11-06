# tests/hang_guard.py

"""
A test that hangs past its limit is stopped and named, and the run goes on (TEST.621).

`pytest-timeout` (`timeout_method = signal`, `pytest.ini`) raises in the test when its limit is up, but only once: a test that
holds a `ProcessPoolExecutor` or a `multiprocessing` pool then waits for its children to finish as the `with` block closes
(`shutdown(wait=True)`), and a child that is busy, or that was forked holding a lock, never does. The alarm has been used, so
nothing stops the wait: the worker sits in a futex with the load at 0 for hours and the run never ends.

This plugin gives every test a second, independent clock that runs on a thread of its own, over all three phases of the test.
Past its limit plus `hang_guard_grace` seconds (120 by default):

1. it prints which test, on which worker, has run how long, and kills every child process of the worker (SIGKILL). The wait on
   the pool ends at once, the timeout `pytest-timeout` raised is reported as the test's failure, and the run goes on;
2. if the test is still running `hang_guard_grace` seconds after that (it waits on something other than a child), it prints so
   and ends the worker with `os._exit`. Under xdist the controller reports the worker as crashed while running that test and
   starts another; without xdist the run ends with the message.

A test with its own `@pytest.mark.timeout(n)` is held to n; `timeout = 0` (no limit) is not watched. Loaded with
`-p tests.hang_guard` (`pytest.ini`), like `tests.depmap`. Linux only (children are found in /proc); elsewhere the
first step finds none and the second still ends the worker.
"""

import os
import signal
import sys
import threading
import time

import pytest

DEFAULT_GRACE_S = 120.0
"""float: Seconds past a test's limit before its children are killed, and again before the worker is ended."""

EXIT_CODE = 70
"""int: The status of a worker ended for a test that would not stop (EX_SOFTWARE)."""


def pytest_addoption(parser):
    parser.addini("hang_guard_grace", "seconds past a test's timeout before its child processes are killed, and again before "
                                      "its worker is ended (TEST.621)", default=str(DEFAULT_GRACE_S))


def children_of(pid):
    """Every descendant process id of `pid`, from /proc, deepest last; empty where /proc is not there."""
    parents = {}
    try:
        names = os.listdir("/proc")
    except OSError:
        return []
    for name in names:
        if not name.isdigit():
            continue
        try:
            with open(f"/proc/{name}/stat", "rb") as handle:
                text = handle.read().decode("utf-8", "replace")
            parents[int(name)] = int(text[text.rindex(")") + 2:].split()[1])      # the fields after the command's closing bracket: state, ppid
        except (OSError, ValueError, IndexError):
            continue
    found = []
    frontier = [pid]
    while frontier:
        nxt = [child for child, parent in parents.items() if parent in frontier and child not in found]
        found.extend(nxt)
        frontier = nxt
    return found


def kill_children(pid=None):
    """Kills every descendant of this process (or `pid`); returns how many it signalled."""
    count = 0
    for child in reversed(children_of(os.getpid() if pid is None else pid)):
        try:
            os.kill(child, signal.SIGKILL)
            count += 1
        except OSError:
            pass
    return count


def say(text):
    sys.stderr.write(text + "\n")
    sys.stderr.flush()


class Watch:
    """The second clock of one test: kills the children when `limit_s` has passed, ends the process `grace_s` after that."""

    def __init__(self, nodeid, limit_s, grace_s, worker="main"):
        self.nodeid = nodeid
        self.limit_s = limit_s
        self.grace_s = grace_s
        self.worker = worker
        self.started = time.monotonic()
        self.killed = 0
        self.stage = 0
        self.timer = None

    def start(self):
        self._arm(self.limit_s + self.grace_s)

    def stop(self):
        if self.timer is not None:
            self.timer.cancel()

    def _arm(self, seconds):
        self.timer = threading.Timer(seconds, self._fire)
        self.timer.daemon = True
        self.timer.start()

    def _fire(self):
        elapsed = time.monotonic() - self.started
        if self.stage == 0:
            self.stage = 1
            self.killed = kill_children()
            say(f"\nHUNG TEST {self.nodeid} (worker {self.worker}) is still running after {elapsed:.0f} s, its limit being "
                f"{self.limit_s:.0f} s plus {self.grace_s:.0f} s: killed {self.killed} child process(es) of the worker (TEST.621).")
            self._arm(self.grace_s)
            return
        say(f"\nHUNG TEST {self.nodeid} (worker {self.worker}) is still running after {elapsed:.0f} s and waits on nothing a child "
            f"process can end: ending the worker (TEST.621).")
        os._exit(EXIT_CODE)


def limit_of(item):
    """The seconds the test may take: its own `timeout` mark, else the `timeout` option; 0 when it has none."""
    marker = item.get_closest_marker("timeout")
    if marker is not None and marker.args:
        return float(marker.args[0])
    if marker is not None and "timeout" in marker.kwargs:
        return float(marker.kwargs["timeout"])
    try:
        return float(item.config.getini("timeout") or 0)
    except (ValueError, KeyError):
        return 0.0


def grace_of(config):
    try:
        return max(float(config.getini("hang_guard_grace")), 0.0)
    except (ValueError, KeyError):
        return DEFAULT_GRACE_S


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_protocol(item, nextitem):
    limit = limit_of(item)
    watch = None
    if limit > 0:
        watch = Watch(item.nodeid, limit, grace_of(item.config), os.environ.get("PYTEST_XDIST_WORKER", "main"))
        watch.start()
    try:
        yield
    finally:
        if watch is not None:
            watch.stop()
