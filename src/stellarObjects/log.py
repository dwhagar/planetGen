"""
Unified CLI Logging
====================

A single, shared logging channel for `generate.py` and every module under
`stellarObjects` -- everything routes through here instead of bare
`print()` calls, so `--debug`/`--quiet`/`--silent` behave consistently
everywhere.

Three severities, matched one-to-one to stdlib `logging` levels so
threshold filtering is just `Logger.setLevel`, not reinvented:

    SILENT  (logging.ERROR)    -- only `error(...)` is shown.
    NORMAL  (logging.INFO)     -- the default; matches the plain, untimed
                                   text this program has always printed.
    DEBUG   (logging.DEBUG)    -- everything, always timestamped, and
                                   optionally mirrored to a file.

`configure()` is called exactly once, from `generate.py`'s `main()`, right
after argument parsing. Every other module just does
`from stellarObjects import log` and calls `log.debug(...)`/`log.normal(...)`/
`log.error(...)`.

`set_console`/`reset_console` let a caller with a live `rich.progress.Progress`
bar redirect this module's console output through `progress.console.print(...)`
instead of a raw stdout write, which is required while that bar is live (see
`generate.py`'s `_generation_progress()` docstring) -- a plain stdout write
fights the bar's own redraw.
"""

import logging
import sys

SILENT = "silent"
NORMAL = "normal"
DEBUG = "debug"

_LEVELS = {
    SILENT: logging.ERROR,
    NORMAL: logging.INFO,
    DEBUG: logging.DEBUG,
}

_TIMESTAMPED_FORMAT = "%(asctime)s [%(levelname)s] %(message)s"
_PLAIN_FORMAT = "%(message)s"

_logger = logging.getLogger("planetgen")
_logger.propagate = False


class _ConsoleHandler(logging.Handler):
    """
    Writes to stdout by default; `set_console`/`reset_console` retarget it
    to a `rich.Console` (a live `Progress`'s own `console`) so output stays
    redraw-safe while a progress bar is on screen.
    """

    def __init__(self):
        super().__init__()
        self._console = None

    def set_console(self, console):
        self._console = console

    def emit(self, record):
        message = self.format(record)
        if self._console is not None:
            self._console.print(message)
        else:
            print(message, file=sys.stdout)


_console_handler = _ConsoleHandler()


def configure(level=NORMAL, debug_file=None):
    """
    (Re)configures the shared logger for the given severity `level`
    (`SILENT`/`NORMAL`/`DEBUG`). `debug_file`, if given, is opened as an
    additional handler that mirrors everything the console handler shows --
    only meaningful (and only ever passed) when `level == DEBUG`.

    Safe to call more than once (handlers are cleared and rebuilt each
    time), which matters for tests that configure the logger repeatedly.

    Args:
        level (str): One of `SILENT`, `NORMAL`, `DEBUG`.
        debug_file (str | None): Path to also write every logged line to.
    """
    for handler in list(_logger.handlers):
        _logger.removeHandler(handler)
        if handler is not _console_handler:
            handler.close()

    _logger.setLevel(_LEVELS[level])

    formatter = logging.Formatter(_TIMESTAMPED_FORMAT if level == DEBUG else _PLAIN_FORMAT)
    _console_handler.setFormatter(formatter)
    _logger.addHandler(_console_handler)

    if debug_file:
        file_handler = logging.FileHandler(debug_file)
        file_handler.setFormatter(formatter)
        _logger.addHandler(file_handler)


def set_console(console):
    """Routes console output through `console.print(...)` (a live `rich.progress.Progress`'s own console)."""
    _console_handler.set_console(console)


def reset_console():
    """Reverts console output to a plain stdout write."""
    _console_handler.set_console(None)


def debug(message, *args, **kwargs):
    _logger.debug(message, *args, **kwargs)


def normal(message, *args, **kwargs):
    _logger.info(message, *args, **kwargs)


def error(message, *args, **kwargs):
    _logger.error(message, *args, **kwargs)


def choice(description, chosen, reason):
    """Debug-only convenience for narrating a randomized/algorithmic decision: what was chosen and why."""
    debug(f"{description}: chose {chosen!r} ({reason})")


# Sensible default so anything that imports stellarObjects modules directly
# (tests, `python -c`, etc.) without going through `generate.py`'s `main()`
# still gets today's plain, untimed output rather than silence.
configure(NORMAL)
