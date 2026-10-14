# tests/failure_capture.py

"""
Keeps every test failure's full report (TEST.139), so a failure that "passed
when I re-ran it alone" leaves something to read.

Off unless `PLANETGEN_FAILURE_LOG` names a file (`scripts/test.py` sets it for
every run). Each failing phase of a test appends one JSON line: the node id,
the phase, how long it ran, the machine's load average when it ended, the
worker, and the full text of the failure. The runner then re-runs the failed
tests alone, one worker, and tells a failure that only happens under load from
one that happens alone (`scripts/test.py`, "Failures").
"""

import json
import os
import time

import pytest


def _line(report, worker):
    try:
        load = round(os.getloadavg()[0], 1)
    except (AttributeError, OSError):
        load = None
    return {
        "nodeid": report.nodeid, "when": report.when, "outcome": report.outcome, "seconds": round(report.duration, 2),
        "load": load, "worker": worker, "at": time.strftime("%H:%M:%SZ", time.gmtime()),
        "text": str(report.longrepr) if report.longrepr is not None else "",
    }


class FailureCapture:
    def __init__(self, path):
        self.path = path

    @pytest.hookimpl(trylast=True)
    def pytest_runtest_logreport(self, report):
        if report.outcome != "failed":
            return
        record = _line(report, os.environ.get("PYTEST_XDIST_WORKER", "main"))
        data = (json.dumps(record) + "\n").encode("utf-8")
        fd = os.open(self.path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
        try:
            os.write(fd, data)  # one write per line: lines from several workers do not mix
        finally:
            os.close(fd)


def register(config):
    path = os.environ.get("PLANETGEN_FAILURE_LOG")
    # Under xdist the controller sees every report the workers send; only the workers write.
    controller = bool(getattr(config.option, "numprocesses", None)) and not hasattr(config, "workerinput")
    if path and not controller:
        config.pluginmanager.register(FailureCapture(path), "planetgen-failure-capture")
