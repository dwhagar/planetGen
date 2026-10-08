# tests/test_console_log.py

"""ADM.23: console log lines are never hard-wrapped, including while a
progress bar's `rich` console carries them (`log.set_console`)."""

import io

import pytest
from rich.console import Console

from planetgen.util import log

LONG = "Generated sector [3.40.7.0] " + "with a long message " * 20 + "[/done]"


@pytest.fixture
def console():
    out = io.StringIO()
    log.configure(log.NORMAL)
    log.set_console(Console(file=out, width=40))
    yield out
    log.reset_console()
    log.configure(log.NORMAL)


def test_a_long_line_stays_one_line_on_a_rich_console(console):
    log.normal(LONG)
    assert console.getvalue() == LONG + "\n"


def test_brackets_are_text_not_rich_markup(console):
    log.normal("[bold]not bold[/bold] [/x]")
    assert console.getvalue() == "[bold]not bold[/bold] [/x]\n"


def test_a_long_line_stays_one_line_on_stdout(capsys):
    log.configure(log.NORMAL)
    log.normal(LONG)
    assert capsys.readouterr().out == LONG + "\n"
