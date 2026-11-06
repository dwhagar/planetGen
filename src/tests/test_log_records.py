"""UX.122: the command line's log lines go through the logging layer, and a record file holds the same lines as records."""

import re
import argparse

import pytest

from planetgen.util import log, logrecord


@pytest.fixture(autouse=True)
def _reset():
    yield
    log.set_record_file(None)
    log.configure(log.NORMAL)


def test_the_console_text_is_what_the_program_has_always_printed(capsys):
    log.configure(log.NORMAL)
    log.normal("Placed 12 sectors")
    log.error("Error: no such sector")
    log.error("The seed found no sector")
    assert capsys.readouterr().out.splitlines() == [
        "Placed 12 sectors", "Error: no such sector", "error: The seed found no sector"]


def test_a_record_file_holds_each_line_as_a_record_the_console_text_is_the_plain_line_of(tmp_path, capsys):
    path = tmp_path / "records.log"
    log.set_record_file(str(path))
    log.configure(log.NORMAL)
    log.normal("Placed 12 sectors")
    log.error("Error: no such sector")
    log.debug("a decision")                                  # below the console's level: not in the file
    for handler in log._logger.handlers:
        handler.flush()
    records = [logrecord.parse_file_line(line) for line in path.read_text().splitlines()]
    assert [(r.severity, r.component, r.message) for r in records] == [
        ("info", "test_log_records", "Placed 12 sectors"), ("error", "test_log_records", "Error: no such sector")]
    assert capsys.readouterr().out.splitlines() == [logrecord.plain_line(r) for r in records]
    assert all(re.match(r"^ts=\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{3}[+-]\d{4} level=", line) for line in path.read_text().splitlines())


def test_a_debug_run_puts_its_debug_lines_in_the_record_file_too(tmp_path):
    path = tmp_path / "records.log"
    log.set_record_file(str(path))
    log.configure(log.DEBUG)
    log.debug("a decision", area="decision")
    for handler in log._logger.handlers:
        handler.flush()
    assert any(logrecord.parse_file_line(line).severity == "debug" and logrecord.parse_file_line(line).message == "a decision"
               for line in path.read_text().splitlines())


def test_the_record_file_can_be_json_lines(tmp_path):
    path = tmp_path / "records.jsonl"
    log.set_record_file(str(path), as_json=True)
    log.configure(log.NORMAL)
    log.normal("Placed 12 sectors")
    for handler in log._logger.handlers:
        handler.flush()
    record = logrecord.parse_json_line(path.read_text().splitlines()[0])
    assert record.message == "Placed 12 sectors" and record.severity == "info"


def test_a_command_takes_the_record_file_from_its_options_or_the_environment(tmp_path, monkeypatch):
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("go")
    log.add_debug_arguments(parser)
    args = parser.parse_args(["--record-file", str(tmp_path / "a.log"), "--record-json", "go"])
    assert args.record_file == str(tmp_path / "a.log") and args.record_json
    log.apply_record_arguments(args)
    assert log._record_file == (str(tmp_path / "a.log"), True)
    monkeypatch.setenv("PLANETGEN_RECORD_FILE", str(tmp_path / "env.log"))
    monkeypatch.setenv("PLANETGEN_RECORD_FORMAT", "json")
    log.apply_record_arguments(parser.parse_args(["go"]))
    assert log._record_file == (str(tmp_path / "env.log"), True)
    monkeypatch.delenv("PLANETGEN_RECORD_FILE")
    log.apply_record_arguments(parser.parse_args(["go"]))
    assert log._record_file is None


def test_a_run_log_file_is_file_lines_and_the_admin_run_log_shows_the_same_lines(redis_server, tmp_path):
    from planetgen.queue import redisqueue, runrecord
    conn = redisqueue.connect(redis_server)
    run_id = runrecord.start_run(conn, "Fill", [runrecord.Stage("mass", "Scatter the massive stars", units=1)], now=100.0)
    runrecord.begin_stage(conn, run_id, 0, now=101.0)
    runrecord.unit_done(conn, run_id, 0, 0, seconds=1.0, now=131.0)
    runrecord.finish_run(conn, run_id, "done", now=132.0)
    text = open(runrecord.write_log_file(conn, run_id, tmp_path), encoding="utf-8").read()
    assert text == runrecord.log_text(conn, run_id)
    records = [logrecord.parse_file_line(line) for line in text.splitlines()]
    assert [r.message for r in records] == [
        "run started: Fill", "Stage 1 of 1: Scatter the massive stars", "unit 0 done", "Stage 1 of 1 finished after 30 s",
        "run done"]
    assert [r.fields.get("event") for r in records] == ["run.start", "stage.start", "unit.done", "stage.finish", "run.done"]
    assert records[1].fields == {"event": "stage.start", "stage": 1, "stages": 1}
    assert records[3].fields == {"event": "stage.finish", "stage": 1, "stages": 1, "seconds": 30.0}
    assert logrecord.plain_line(records[1]) == "=== Stage 1 of 1: Scatter the massive stars ==="


TIMED = re.compile(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \w+ (INFO|WARNING|ERROR) +\w+: ")


def _command_args(*extra):
    parser = argparse.ArgumentParser()
    log.add_debug_arguments(parser)
    return parser, parser.parse_args(list(extra))


def test_a_command_whose_output_is_not_a_terminal_prints_every_line_with_its_time(capsys, monkeypatch):
    """OPS.98: piped into a file, `tee` or a log, a command's lines are display lines; capsys is not a terminal."""
    monkeypatch.delenv("PLANETGEN_LOG_TIMES", raising=False)
    parser, args = _command_args()
    log.apply_debug_arguments(args, parser)
    log.normal("Placed 12 sectors")
    log.warning("The seed found no sector")
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 2 and all(TIMED.match(line) for line in lines)
    assert lines[0].endswith("INFO     test_log_records: Placed 12 sectors")


def test_the_lines_stay_plain_on_a_terminal_and_the_switch_can_force_either(capsys, monkeypatch):
    parser, args = _command_args()
    monkeypatch.setenv("PLANETGEN_LOG_TIMES", "never")
    log.apply_debug_arguments(args, parser)
    log.normal("Placed 12 sectors")
    assert capsys.readouterr().out == "Placed 12 sectors\n"
    monkeypatch.setenv("PLANETGEN_LOG_TIMES", "always")
    log.apply_debug_arguments(args, parser)
    log.normal("Placed 12 sectors")
    assert TIMED.match(capsys.readouterr().out)
    monkeypatch.delenv("PLANETGEN_LOG_TIMES")
    monkeypatch.setattr("sys.stdout", type("Terminal", (), {"isatty": lambda self: True, "write": lambda self, text: len(text), "flush": lambda self: None})())
    assert log.times_on_console() is False


def test_a_plain_configure_keeps_the_plain_line_so_only_a_command_line_asks_for_times(capsys, monkeypatch):
    monkeypatch.setenv("PLANETGEN_LOG_TIMES", "always")
    log.configure(log.NORMAL)
    log.normal("Placed 12 sectors")
    assert capsys.readouterr().out == "Placed 12 sectors\n"
