# tests/test_run_progress_width.py

"""
PERF.179: a generation progress bar written to a job log (output that is not a terminal) keeps its whole line. Rich assumes 80 columns
there and cut the last columns off ("elapsed 0:0... remaining 0:00:..."); the bars are drawn `run_common.LOG_CONSOLE_WIDTH` wide instead.
"""

import io
import re
import sys
import time

import pytest

from planetgen.generation import run_common

LINE_END = re.compile(r"elapsed (\d+ days?, )?\d+:\d\d:\d\d remaining (\d+:\d\d:\d\d|-:--:--)$")


def _log_of(monkeypatch, bars):
    """What a run writes to a job log while `bars` (`[(description, completed, total, unit, hours_running)]`) are on screen."""
    monkeypatch.delenv("COLUMNS", raising=False)
    log = io.StringIO()
    monkeypatch.setattr(sys, "stdout", log)
    progress = run_common._generation_progress()
    progress.start()
    for description, completed, total, unit, hours in bars:
        task = progress.add_task(description, total=total, **({"unit": unit} if unit else {"percent": True}))
        progress.update(task, completed=completed)
        progress._tasks[task].start_time = time.monotonic() - hours * 3600.0
    progress.refresh()
    progress.stop()
    return [line for line in log.getvalue().splitlines() if line.strip()]


def test_a_bar_written_to_a_log_ends_in_its_whole_elapsed_and_remaining_times(monkeypatch):
    lines = _log_of(monkeypatch, [("Sectors", 1234, 5000, "sector", 0.01), ("Stars of the neighbourhood (12,345 of 98,765 shells)", 37.0, 100.0, None, 0.5)])
    assert len(lines) == 2
    for line in lines:
        assert "…" not in line and "..." not in line, line
        assert LINE_END.search(line), line


def test_the_longest_line_a_run_draws_still_fits_the_log_width(monkeypatch):
    """A long step name, a count in the hundreds of thousands and a run of over a day (rich writes it as days, then hours)."""
    bars = [("Neighbourhood sectors, systems and their stars by territory", 123_456, 987_654, "sector", 30.0)]
    [line] = _log_of(monkeypatch, bars)
    assert len(line) <= run_common.LOG_CONSOLE_WIDTH and LINE_END.search(line), line
    assert "1 day, 6:00:0" in line


def test_a_terminal_or_a_columns_setting_is_left_to_rich(monkeypatch):
    monkeypatch.setenv("COLUMNS", "100")
    assert run_common._log_console() is None
    monkeypatch.delenv("COLUMNS")
    assert run_common._log_console().width == run_common.LOG_CONSOLE_WIDTH
