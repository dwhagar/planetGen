"""OPS.167: a warning is logged as a warning. A line that says "Warning:" at normal level never reaches a warning reader (error.log, the summary's warning list)."""

import re
from pathlib import Path

import planetgen

SOURCE = Path(planetgen.__file__).parent
SEVERITY_WORD_AT_NORMAL = re.compile(r"""log\.(?:normal|info)\(\s*f?["'](?:Warning|WARNING|warning)\b""")


def test_no_line_says_warning_at_normal_level():
    offenders = []
    for path in SOURCE.rglob("*.py"):
        for number, line in enumerate(path.read_text().splitlines(), 1):
            if SEVERITY_WORD_AT_NORMAL.search(line):
                offenders.append(f"{path.relative_to(SOURCE)}:{number}")
    assert not offenders, f"use log.warning (the severity says it; leave the word out of the text): {offenders}"


class _Conn:
    def rollback(self):
        pass


class _Gate:
    def is_closed(self):
        return False

    def grant(self, need):
        return True

    def release(self):
        pass


def test_a_deadlock_retry_is_a_notice_and_a_granted_try_that_fails_too_is_a_warning(caplog, monkeypatch):
    import logging

    import pymysql
    import pytest

    from planetgen.db import contention
    from planetgen.util import logrecord
    monkeypatch.setattr(contention.time, "sleep", lambda seconds: None)

    def deadlock(conn):
        raise pymysql.err.OperationalError(1213, "Deadlock found")

    with caplog.at_level(1), pytest.raises(pymysql.err.OperationalError):
        contention.run_with_retries(_Conn(), deadlock, "gate", "Unit 4", attempts=2, gate_object=_Gate())
    retry = [r for r in caplog.records if "trying again" in r.getMessage()]
    final = [r for r in caplog.records if "granted try failed too" in r.getMessage()]
    assert retry and all(r.levelno == logrecord.LOGGING_LEVELS["notice"] for r in retry)
    assert final and all(r.levelno == logging.WARNING for r in final)
