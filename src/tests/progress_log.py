# tests/progress_log.py

"""
A line for every test that finishes (TEST.150, step one of resuming a run, TEST.148).

Off unless `PLANETGEN_PROGRESS_LOG` names a file (`scripts/test.py` sets it for every run, to
`.test-progress/<build id>.jsonl`). When a test's last phase (teardown) is reported, one JSON line is appended:
the node id, its outcome, how long its phases took together, the worker and the time. A run that is killed loses at
most the tests that were running at that moment: they have no line, so a resume runs them again.

Outcomes, worst first: `failed` (any phase), `skipped`, `xfailed`, `xpassed`, `passed`. Only `passed` says the
test is done; the rest run again on a resume.

Loaded with `-p tests.progress_log` (`pytest.ini`), like `tests.depmap`. Under xdist only the workers write, one
`write` per line, so lines from several workers do not mix.
"""

import json
import os
import time

import pytest

ORDER = ("failed", "skipped", "xfailed", "xpassed", "passed")
"""tuple: Outcomes from the one a test is judged by first (a failing teardown fails a passed test) to the least."""


def phase_outcome(report):
    """The outcome one phase report stands for."""
    if hasattr(report, "wasxfail"):
        return "xfailed" if report.outcome == "skipped" else "xpassed" if report.outcome == "passed" else report.outcome
    return report.outcome


def worst(outcomes):
    return min(outcomes, key=lambda outcome: ORDER.index(outcome) if outcome in ORDER else 0)


class ProgressLog:
    def __init__(self, path, worker="main"):
        self.path = path
        self.worker = worker
        self.running = {}  # node id -> (outcomes so far, seconds so far)

    @pytest.hookimpl(trylast=True)
    def pytest_runtest_logreport(self, report):
        outcomes, seconds = self.running.get(report.nodeid, ([], 0.0))
        outcomes = outcomes + [phase_outcome(report)]
        seconds += report.duration
        if report.when != "teardown":
            self.running[report.nodeid] = (outcomes, seconds)
            return
        self.running.pop(report.nodeid, None)
        line = {"kind": "test", "nodeid": report.nodeid, "outcome": worst(outcomes), "seconds": round(seconds, 2),
                "worker": self.worker, "at": time.strftime("%H:%M:%SZ", time.gmtime())}
        data = (json.dumps(line) + "\n").encode("utf-8")
        fd = os.open(self.path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
        try:
            os.write(fd, data)
        finally:
            os.close(fd)


def pytest_configure(config):
    path = os.environ.get("PLANETGEN_PROGRESS_LOG")
    # Under xdist the controller sees every report the workers send; only the workers write.
    controller = bool(getattr(config.option, "numprocesses", None)) and not hasattr(config, "workerinput")
    if path and not controller:
        config.pluginmanager.register(ProgressLog(path, os.environ.get("PYTEST_XDIST_WORKER", "main")), "planetgen-progress-log")


def read(path):
    """(header, {node id: its last line}) of a progress log; ({}, {}) when it is missing. A torn last line (the
    process was killed in the middle of a write) is ignored."""
    header, tests = {}, {}
    try:
        with open(path, encoding="utf-8") as f:
            for text in f:
                try:
                    line = json.loads(text)
                except ValueError:
                    continue
                if line.get("kind") == "header":
                    header = line
                elif line.get("kind") == "test":
                    tests[line["nodeid"]] = line
    except OSError:
        pass
    return header, tests
