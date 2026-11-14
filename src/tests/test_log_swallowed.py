"""OPS.168: `log.swallowed` writes one debug line and nothing at normal level."""

import logging

from planetgen.util import log


def test_a_swallowed_error_is_one_debug_line_and_never_a_warning(caplog):
    with caplog.at_level(1), caplog.at_level(1, logger=log._logger.name):
        log.swallowed("queue/scheduler", "_write_usage", RuntimeError("redis gone"))
    (record,) = [r for r in caplog.records if "went on past" in r.getMessage()]
    assert record.levelno == logging.DEBUG
    assert "queue/scheduler: went on past _write_usage: RuntimeError: redis gone" in record.getMessage()


def test_a_swallowed_error_costs_nothing_when_debug_is_off(caplog, monkeypatch):
    class Boom:
        def __str__(self):
            raise AssertionError("the text was built although nobody records it")

    monkeypatch.setattr(log._logger, "isEnabledFor", lambda level: False)
    log.swallowed("db/stats", "server_info", Boom())
