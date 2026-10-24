"""UX.121: the logging layer (`planetgen/util/logrecord.py`). The same records give a file line, a JSON line and a display
line that agree, and a file line reads back to the record it was written from."""

import datetime
import io
import json
import logging
import zoneinfo

import pytest

from planetgen.util import logrecord

NEW_YORK = zoneinfo.ZoneInfo("America/New_York")


def _ms(year, month, day, hour, minute, second, milli):
    moment = datetime.datetime(year, month, day, hour, minute, second, tzinfo=datetime.timezone.utc)
    return int(moment.timestamp()) * 1000 + milli


OCTOBER = _ms(2026, 10, 24, 18, 40, 0, 123)
JANUARY = _ms(2027, 1, 5, 18, 40, 0, 7)

RECORDS = [
    logrecord.Record(OCTOBER, "info", "generate", "Placed the sectors", {"count": 12, "elapsed": 1.5}),
    logrecord.Record(OCTOBER, "debug", "draw", "", {}),
    logrecord.Record(OCTOBER, "warning", "queue", 'He said "stop" = now \\ then', {"who": 'a "b" c', "n": 0}),
    logrecord.Record(OCTOBER, "error", "web", "two\nlines\tand a tab\x01", {"path": "/admin/log?x=1&y=2"}),
    logrecord.Record(OCTOBER, "notice", "update", "Zeta ζ and   separator", {"flag": True, "off": False}),
    logrecord.Record(OCTOBER, "critical", "db", "x", {"text_number": "12", "text_float": "1.5", "text_true": "true",
                                                      "text_empty": "", "text_dash": "-"}),
    logrecord.Record(JANUARY, "info", "generate", "A winter record", {"big": 10 ** 12, "tiny": 1e-07, "whole": 3.0}),
]


def test_a_file_line_is_logfmt_with_the_time_in_utc():
    record = logrecord.Record(OCTOBER, "info", "generate", "Stage 4 of 12: Scatter", {"stage": 4, "stages": 12})
    assert logrecord.file_line(record) == (
        'ts=2026-10-24T18:40:00.123Z level=info component=generate msg="Stage 4 of 12: Scatter" stage=4 stages=12')


def test_a_json_line_has_the_same_names_with_the_fields_beside_them():
    record = logrecord.Record(OCTOBER, "warning", "queue", "Slow", {"seconds": 2.5, "flag": True})
    assert json.loads(logrecord.json_line(record)) == {
        "ts": "2026-10-24T18:40:00.123Z", "level": "warning", "component": "queue", "msg": "Slow", "seconds": 2.5,
        "flag": True}


@pytest.mark.parametrize("record", RECORDS)
def test_every_line_is_one_line_and_reads_back_to_its_record(record):
    for text, parse in ((logrecord.file_line(record), logrecord.parse_file_line),
                        (logrecord.json_line(record), logrecord.parse_json_line)):
        assert "\n" not in text and "\r" not in text and " " not in text and "\x01" not in text
        assert parse(text) == record
        assert [type(value) for value in parse(text).fields.values()] == [type(value) for value in record.fields.values()]


@pytest.mark.parametrize("record", RECORDS)
def test_the_screen_shows_what_the_file_holds(record):
    from_file = logrecord.parse_file_line(logrecord.file_line(record))
    from_json = logrecord.parse_json_line(logrecord.json_line(record))
    assert logrecord.display_line(from_file, NEW_YORK) == logrecord.display_line(record, NEW_YORK)
    assert logrecord.display_line(from_json, NEW_YORK) == logrecord.display_line(record, NEW_YORK)
    assert logrecord.plain_line(from_file) == logrecord.plain_line(record)


def test_a_display_line_is_local_time_labelled_with_its_zone():
    record = logrecord.Record(OCTOBER, "info", "generate", "Placed the sectors", {"count": 12})
    assert logrecord.display_line(record, NEW_YORK) == (
        "2026-10-24 14:40:00 EDT INFO     generate: Placed the sectors (count=12)")
    winter = logrecord.Record(JANUARY, "error", "web", "Down")
    assert logrecord.display_line(winter, NEW_YORK) == "2027-01-05 13:40:00 EST ERROR    web: Down"
    assert logrecord.display_line(winter, datetime.timezone.utc) == "2027-01-05 18:40:00 UTC ERROR    web: Down"
    assert logrecord.display_line(record, NEW_YORK, times=False, component=False, fields=False) == "INFO     Placed the sectors"


def test_the_command_line_text_is_the_message_and_a_warning_says_what_it_is():
    assert logrecord.plain_line(logrecord.Record(OCTOBER, "info", "generate", "Placed 12")) == "Placed 12"
    assert logrecord.plain_line(logrecord.Record(OCTOBER, "notice", "generate", "Placed 12")) == "Placed 12"
    assert logrecord.plain_line(logrecord.Record(OCTOBER, "warning", "generate", "Slow")) == "warning: Slow"
    assert logrecord.plain_line(logrecord.Record(OCTOBER, "error", "generate", "Broke")) == "error: Broke"


def test_a_stage_has_one_wording_for_start_finish_and_failure_and_one_divider():
    start = logrecord.stage_record(logrecord.STAGE_START, "generate", 4, 12, "Scatter the massive stars", now=1.0)
    done = logrecord.stage_record(logrecord.STAGE_FINISH, "generate", 4, 12, seconds=30.6, now=2.0)
    failed = logrecord.stage_record(logrecord.STAGE_FAILED, "generate", 4, 12, seconds=31, reason="MySQL went away", now=3.0)
    run = logrecord.stage_record(logrecord.STAGE_START, "generate", 4, 12, "Settle", count=9, now=4.0)
    assert [item.message for item in (start, done, failed, run)] == [
        "Stage 4 of 12: Scatter the massive stars", "Stage 4 of 12 finished after 31 s",
        "Stage 4 of 12 failed after 31 s: MySQL went away", "Stages 4 to 12 of 12: Settle"]
    assert [item.severity for item in (start, done, failed)] == ["info", "info", "error"]
    assert start.fields == {"event": "stage.start", "stage": 4, "stages": 12}
    assert done.fields["seconds"] == 30.6 and run.fields["count"] == 9
    assert logrecord.plain_line(start) == "=== Stage 4 of 12: Scatter the massive stars ==="
    assert logrecord.display_line(failed, NEW_YORK, times=False) == (
        "ERROR    generate: === Stage 4 of 12 failed after 31 s: MySQL went away ===")
    assert logrecord.parse_file_line(logrecord.file_line(failed)) == failed
    with pytest.raises(ValueError, match="no stage event"):
        logrecord.stage_record("stage.begin", "generate", 1, 2)


def test_the_stage_wording_is_the_stage_markers():
    from planetgen.util import stagemark
    assert logrecord.stage_start(4, 12, "x").startswith(stagemark.marker(4, 12))
    assert logrecord.stage_start(4, 12, "x", count=9).startswith(stagemark.marker(4, 12, 9))


def test_the_time_is_whole_milliseconds_in_utc():
    record = logrecord.Record.at(1_792_867_200.1239, "info", "generate", "x")
    assert record.ms % 1000 == 124 and logrecord.file_line(record).startswith("ts=2026-10-24T18:40:00.124Z ")
    assert logrecord.Record.now("info", "generate", "x").ms > OCTOBER


def test_severities_are_the_syslog_levels_in_order_and_map_to_logging():
    assert logrecord.SEVERITIES == ("debug", "info", "notice", "warning", "error", "critical")
    assert [logrecord.SYSLOG_CODES[name] for name in logrecord.SEVERITIES] == [7, 6, 5, 4, 3, 2]
    levels = [logrecord.LOGGING_LEVELS[name] for name in logrecord.SEVERITIES]
    assert levels == sorted(levels) and logging.getLevelName(logrecord.NOTICE) == "NOTICE"


@pytest.mark.parametrize("make", [
    lambda: logrecord.Record(0, "loud", "generate", "x"),
    lambda: logrecord.Record(0, "info", "Generate", "x"),
    lambda: logrecord.Record(0, "info", "", "x"),
    lambda: logrecord.Record(0, "info", "generate", "x", {"msg": "y"}),
    lambda: logrecord.Record(0, "info", "generate", "x", {"Bad Name": 1}),
    lambda: logrecord.Record(0, "info", "generate", "x", {"n": float("nan")}),
    lambda: logrecord.Record(0, "info", "generate", "x", {"n": [1]}),
    lambda: logrecord.Record(0, "info", "generate", "x", {"n": None}),
])
def test_a_record_refuses_what_a_file_line_could_not_hold(make):
    with pytest.raises(ValueError):
        make()


@pytest.mark.parametrize("line", [
    "", "hello", "ts=2026-10-24T18:40:00.123Z level=info component=generate",
    'level=info ts=2026-10-24T18:40:00.123Z component=generate msg="x"',
    'ts=2026-10-24 level=info component=generate msg="x"'])
def test_a_line_that_is_not_a_record_is_refused(line):
    with pytest.raises(ValueError):
        logrecord.parse_file_line(line)


def _handler(formatter):
    stream = io.StringIO()
    handler = logging.StreamHandler(stream)
    handler.setFormatter(formatter)
    return handler, stream


def test_a_logger_takes_fields_and_a_handler_writes_the_file_line_the_screen_line_and_the_json_line():
    log = logrecord.get_logger("layer-test")
    file_handler, file_out = _handler(logrecord.FileFormatter())
    json_handler, json_out = _handler(logrecord.JsonFormatter())
    screen_handler, screen_out = _handler(logrecord.DisplayFormatter(NEW_YORK))
    plain_handler, plain_out = _handler(logrecord.DisplayFormatter(plain=True))
    logger = log.logger
    old_level = logger.level
    logger.setLevel(logging.DEBUG)
    handlers = (file_handler, json_handler, screen_handler, plain_handler)
    for handler in handlers:
        logger.addHandler(handler)
    try:
        log.info("Placed the sectors", count=12, elapsed=1.5)
        log.notice("A notice")
        log.warning("Slow disk", mb=3)
        log.stage(logrecord.STAGE_START, 4, 12, "Scatter")
    finally:
        for handler in handlers:
            logger.removeHandler(handler)
        logger.setLevel(old_level)
    records = [logrecord.parse_file_line(line) for line in file_out.getvalue().splitlines()]
    assert [(r.severity, r.component, r.message) for r in records] == [
        ("info", "layer-test", "Placed the sectors"), ("notice", "layer-test", "A notice"),
        ("warning", "layer-test", "Slow disk"), ("info", "layer-test", "Stage 4 of 12: Scatter")]
    assert records[0].fields == {"count": 12, "elapsed": 1.5} and records[3].fields["event"] == "stage.start"
    assert [logrecord.parse_json_line(line) for line in json_out.getvalue().splitlines()] == records
    assert screen_out.getvalue().splitlines() == [logrecord.display_line(r, NEW_YORK) for r in records]
    assert plain_out.getvalue().splitlines() == [
        "Placed the sectors", "A notice", "warning: Slow disk", "=== Stage 4 of 12: Scatter ==="]


def test_a_field_may_be_named_like_an_argument_of_the_logging_calls():
    record = logrecord.Record.at(1.0, "info", "generate", "x", severity="high", message="m", component_name="c")
    assert record.fields == {"severity": "high", "message": "m", "component_name": "c"}
    assert logrecord.parse_file_line(logrecord.file_line(record)) == record


def test_a_logger_refuses_a_field_a_file_line_could_not_hold_and_a_bad_component():
    log = logrecord.get_logger("layer-test")
    log.logger.setLevel(logging.DEBUG)
    with pytest.raises(ValueError):
        log.info("x", msg="y")
    with pytest.raises(ValueError):
        logrecord.get_logger("Bad Name")


def test_a_standard_log_record_from_any_logger_becomes_a_record():
    log_record = logging.LogRecord("planetgen.queue.units", logging.WARNING, "f.py", 1, "Unit %s slow", ("7",), None)
    record = logrecord.record_of(log_record)
    assert (record.severity, record.component, record.message) == ("warning", "queue.units", "Unit 7 slow")
    odd = logrecord.record_of(logging.LogRecord("urllib3.connectionpool", 5, "f.py", 1, "x", None, None))
    assert (odd.severity, odd.component) == ("debug", "urllib3.connectionpool")
