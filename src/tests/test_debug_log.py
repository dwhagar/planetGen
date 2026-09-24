"""
Tests for the debug log (`config.json`'s `"debug"`, see `stellarObjects/log.py`):
it's written only when debug is on, carries timestamps, errors and every
random roll, never changes what a seeded run generates, and never breaks
the program when the file can't be opened.
"""

import logging
import os
import random
import re
import subprocess
import sys

import pytest

from stellarObjects import appconfig, log

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_HTML_LIB_DIR = os.path.join(_SRC_DIR, "html", "lib")

ROLL_CHANCE = 0.25


@pytest.fixture
def debug_log(tmp_path, monkeypatch):
    """Turns debug on with the log in a temp file; puts logging back as it was afterwards."""
    path = tmp_path / "planetgen.log"
    monkeypatch.setenv("PLANETGEN_DEBUG", "1")
    monkeypatch.setenv("PLANETGEN_LOG_FILE", str(path))
    log.configure(log.NORMAL)
    yield path
    monkeypatch.delenv("PLANETGEN_DEBUG")
    monkeypatch.delenv("PLANETGEN_LOG_FILE")
    log.configure(log.NORMAL)


def _flush():
    for handler in logging.getLogger("planetgen").handlers:
        handler.flush()


def test_debug_enabled_reads_env_then_config(monkeypatch):
    monkeypatch.delenv("PLANETGEN_DEBUG", raising=False)
    assert appconfig.debug_enabled({"debug": True}) is True
    assert appconfig.debug_enabled({"debug": False}) is False
    assert appconfig.debug_enabled({}) is False
    for value, expected in (("1", True), ("true", True), ("0", False), ("off", False), ("", False)):
        monkeypatch.setenv("PLANETGEN_DEBUG", value)
        assert appconfig.debug_enabled({"debug": not expected}) is expected


def test_log_file_path_defaults_to_var_log(monkeypatch):
    monkeypatch.delenv("PLANETGEN_LOG_FILE", raising=False)
    assert appconfig.log_file_path({}) == "/var/log/planetgen.log"
    assert appconfig.log_file_path({"log_file": "/tmp/x.log"}) == "/tmp/x.log"
    monkeypatch.setenv("PLANETGEN_LOG_FILE", "/srv/pg.log")
    assert appconfig.log_file_path({"log_file": "/tmp/x.log"}) == "/srv/pg.log"


def test_no_log_file_when_debug_is_off(tmp_path, monkeypatch):
    path = tmp_path / "planetgen.log"
    monkeypatch.setenv("PLANETGEN_DEBUG", "0")
    monkeypatch.setenv("PLANETGEN_LOG_FILE", str(path))
    log.configure(log.NORMAL)
    log.debug("nothing should record this")
    log.error("nor this")
    assert not log.debug_log_active()
    assert not path.exists()
    assert not getattr(random.random, "_planetgen_traced", False)


def test_debug_log_records_timestamped_lines_with_their_source(debug_log):
    assert log.debug_log_active()
    log.debug("a narrated decision")
    log.choice("Widget color", "blue", "the only color left")
    _flush()
    text = debug_log.read_text()
    assert "Debug log opened" in text
    line = next(l for l in text.splitlines() if "a narrated decision" in l)
    assert re.match(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{3}[+-]\d{4} DEBUG ", line)
    assert f"[{os.getpid()}]" in line
    assert "test_debug_log.py:" in line
    assert "test_debug_log_records_timestamped_lines_with_their_source" in line
    choice_line = next(l for l in text.splitlines() if "Widget color" in l)
    assert "test_debug_log.py:" in choice_line
    assert "chose 'blue' (the only color left)" in choice_line


def test_passwords_on_the_command_line_are_withheld(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["generate.py", "sector", "--mysql-password", "hunter2",
                                      "--mysql-password=hunter3", "--name", "X"])
    assert log._redacted_argv() == ["generate.py", "sector", "--mysql-password", "<withheld>",
                                    "--mysql-password=<withheld>", "--name", "X"]


def test_errors_and_tracebacks_reach_the_debug_log(debug_log):
    try:
        raise ValueError("boom")
    except ValueError:
        log.exception("while doing a thing")
    logging.getLogger("some.library").warning("library warning")
    _flush()
    text = debug_log.read_text()
    assert "ERROR" in text and "while doing a thing" in text
    assert "Traceback (most recent call last)" in text and "ValueError: boom" in text
    assert "library warning" in text


def test_random_rolls_are_traced_with_the_values_they_reference(debug_log):
    if random.random() < ROLL_CHANCE:
        pass
    random.choices(["a", "b"], weights=[1, 3], k=1)
    _flush()
    rolls = [l for l in debug_log.read_text().splitlines() if " roll random." in l]
    assert any("if random.random() < ROLL_CHANCE:" in l and "[ROLL_CHANCE=0.25]" in l for l in rolls)
    assert any("random.choices(['a', 'b'], weights=[1, 3], k=1)" in l for l in rolls)
    assert all("test_debug_log.py:" in l for l in rolls)


def test_tracing_does_not_change_a_seeded_run(tmp_path, monkeypatch):
    def draws():
        random.seed(1234)
        values = [random.random(), random.uniform(1, 9), random.randint(0, 100), random.gauss(0, 1),
                  random.choice("abcdef"), random.choices("xyz", weights=[1, 2, 3], k=4),
                  random.sample(range(50), 5)]
        items = list(range(10))
        random.shuffle(items)
        return values + [items]

    monkeypatch.setenv("PLANETGEN_DEBUG", "0")
    log.configure(log.NORMAL)
    plain = draws()
    monkeypatch.setenv("PLANETGEN_DEBUG", "1")
    monkeypatch.setenv("PLANETGEN_LOG_FILE", str(tmp_path / "planetgen.log"))
    log.configure(log.NORMAL)
    try:
        assert getattr(random.random, "_planetgen_traced", False)
        traced = draws()
    finally:
        monkeypatch.setenv("PLANETGEN_DEBUG", "0")
        log.configure(log.NORMAL)
    assert traced == plain
    assert not getattr(random.random, "_planetgen_traced", False)


def test_rolls_stay_off_the_console_under_cli_debug(debug_log, capsys):
    log.configure(log.DEBUG)
    try:
        random.random()
        log.debug("console-visible decision")
    finally:
        log.configure(log.NORMAL)
    out = capsys.readouterr().out
    assert "console-visible decision" in out
    assert " roll random." not in out
    _flush()
    assert " roll random." in debug_log.read_text()


def test_unwritable_log_file_does_not_break_the_program(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("PLANETGEN_DEBUG", "1")
    monkeypatch.setenv("PLANETGEN_LOG_FILE", str(tmp_path / "missing-dir" / "planetgen.log"))
    try:
        log.configure(log.NORMAL)
        assert not log.debug_log_active()
        log.normal("still works")
    finally:
        monkeypatch.setenv("PLANETGEN_DEBUG", "0")
        log.configure(log.NORMAL)
    captured = capsys.readouterr()
    assert "still works" in captured.out
    assert "can't be opened" in captured.err


def test_web_page_error_goes_to_the_log_not_the_page(tmp_path):
    path = tmp_path / "planetgen.log"
    env = dict(os.environ, PLANETGEN_DEBUG="1", PLANETGEN_LOG_FILE=str(path), REQUEST_METHOD="GET",
               QUERY_STRING="db=galaxy&password=hunter2", SCRIPT_FILENAME="/srv/html/broken.py",
               REMOTE_ADDR="203.0.113.9")
    code = (
        "import sys; sys.path[:0] = [sys.argv[1], sys.argv[2]]\n"
        "import page\n"
        "def handler():\n"
        "    raise ZeroDivisionError('division by zero in a page')\n"
        "page.run(handler)\n"
    )
    result = subprocess.run([sys.executable, "-c", code, _HTML_LIB_DIR, _SRC_DIR], env=env,
                            capture_output=True, text=True, timeout=60)
    assert "Status: 500" in result.stdout
    assert "ZeroDivisionError" not in result.stdout
    assert "The traceback is in the debug log" in result.stdout
    text = path.read_text()
    assert "web/broken.py[" in text
    assert "Request: GET" in text and "203.0.113.9" in text
    assert "ZeroDivisionError: division by zero in a page" in text
    assert "Response: 500 Internal Server Error" in text
    assert "hunter2" not in text
