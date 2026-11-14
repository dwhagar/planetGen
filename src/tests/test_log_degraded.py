"""OPS.165: a fallback that goes on without Redis leaves a trace: a debug line every time, a warning once in five minutes."""

import logging

import redis

from planetgen.db import countcache
from planetgen.util import log


def _warnings(caplog):
    return [r for r in caplog.records if r.levelno == logging.WARNING and "fallback" in r.getMessage()]


def test_a_fallback_warns_once_in_five_minutes_and_says_debug_every_time(caplog, monkeypatch):
    monkeypatch.setattr(log, "_degraded_said", {})
    now = [1000.0]
    monkeypatch.setattr(log.time, "monotonic", lambda: now[0])
    with caplog.at_level(1), caplog.at_level(1, logger=log._logger.name):
        for _ in range(3):
            log.degraded("pagecache", "reading a page from Redis failed", redis.ConnectionError("refused"))
        assert len(_warnings(caplog)) == 1
        assert sum("reading a page from Redis failed" in r.getMessage() and r.levelno == logging.DEBUG for r in caplog.records) == 3
        log.degraded("pagecache", "sharing a page failed", redis.ConnectionError("refused"))       # another place warns at once
        assert len(_warnings(caplog)) == 2
        now[0] += log.DEGRADED_WARN_SECONDS
        log.degraded("pagecache", "reading a page from Redis failed", redis.ConnectionError("refused"))
        assert len(_warnings(caplog)) == 3
    assert "ConnectionError: refused" in _warnings(caplog)[0].getMessage()


def test_a_silent_redis_fallback_in_the_count_cache_now_warns(caplog, monkeypatch):
    monkeypatch.setattr(log, "_degraded_said", {})

    class Dead:
        def get(self, *args, **kwargs):
            raise redis.ConnectionError("Error 111 connecting to localhost:6379")

    monkeypatch.setattr(countcache, "_redis", lambda: Dead())
    with caplog.at_level(1):
        assert countcache._recall("db1", "k") is None
    (warning,) = _warnings(caplog)
    assert "countcache" in warning.getMessage() and "connecting to localhost:6379" in warning.getMessage()
