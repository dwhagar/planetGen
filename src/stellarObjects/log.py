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

Separately from the console, `config.json`'s `"debug": true` (or
`PLANETGEN_DEBUG=1`) turns on the **debug log**: every entry point --
`generate.py`, the maintenance scripts, the `html/` CGI pages and the
Flask API -- appends everything at DEBUG severity to `/var/log/
planetgen.log` (`"log_file"`/`PLANETGEN_LOG_FILE` to move it), whatever
the console's own level is. Each line carries a millisecond timestamp,
the process (`generate.py[1234]`, `web/system.py[88]`, `api[77]`) and the
source line that logged it. On top of the narrated decisions below, the
debug log also gets:

  - every draw from the global `random` module made by planetGen's own
    code (`_trace_random`): the call, its arguments (weights included),
    the result, and the source line that asked for it -- so a roll like
    `if random.random() < FLAVOR_CHANCE_PLANET:` shows both the number
    drawn and the probability it was compared against;
  - every SQL statement `_db.Connection` runs, with its timing;
  - uncaught exceptions (`sys.excepthook`/`threading.excepthook`),
    Python warnings, and anything other libraries log (Flask, urllib3, ...).

The tracer only wraps calls, it never draws extra numbers, so a seeded run
produces the same result with the debug log on or off. With debug off
none of this is installed and nothing is written.

`set_console`/`reset_console` let a caller with a live `rich.progress.Progress`
bar redirect this module's console output through `progress.console.print(...)`
instead of a raw stdout write, which is required while that bar is live (see
`generate.py`'s `_generation_progress()` docstring) -- a plain stdout write
fights the bar's own redraw.
"""

import functools
import keyword
import linecache
import logging
import logging.handlers
import os
import random
import re
import sys
import threading
import time
from contextlib import contextmanager

from stellarObjects import appconfig

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

# High-volume detail (every random draw, every SQL statement) that belongs
# in the debug log only -- never on the console, even under --debug.
_trace_logger = logging.getLogger("planetgen.trace")
_trace_logger.propagate = False
_trace_logger.setLevel(logging.DEBUG)

_component = None
"""str | None: Short name for this process in the debug log (`set_component`);
defaults to the running script's file name."""

_debug_log_handler = None
_root_stderr_handler = None
_warned_unwritable = set()


def set_component(name):
    """Names this process in the debug log (e.g. `"api"`, `"web/system.py"`)."""
    global _component
    _component = name


def _component_name():
    if _component:
        return _component
    return os.path.basename(sys.argv[0]) if sys.argv and sys.argv[0] else "python"


class _OriginFilter(logging.Filter):
    """Stamps each record with `origin`: `component[pid]`, plus the thread
    name when it isn't the main thread (the API serves requests on threads)."""

    def filter(self, record):
        origin = f"{_component_name()}[{record.process}]"
        if record.threadName != "MainThread":
            origin += f"/{record.threadName}"
        record.origin = origin
        return True


class _DebugLogFormatter(logging.Formatter):
    """Local time with milliseconds and UTC offset, e.g.
    `2026-09-24T02:57:37.123-0400`, so lines from several processes (and
    across a DST change) still sort and compare correctly."""

    def formatTime(self, record, datefmt=None):
        local = time.localtime(record.created)
        return (time.strftime("%Y-%m-%dT%H:%M:%S", local)
                + f".{int(record.msecs):03d}" + time.strftime("%z", local))


_DEBUG_LOG_FORMAT = ("%(asctime)s %(levelname)-8s %(origin)s %(name)s "
                     "%(filename)s:%(lineno)d %(funcName)s | %(message)s")


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


def configure(level=NORMAL, debug_file=None, console=True):
    """
    (Re)configures the shared logger for the given console severity `level`
    (`SILENT`/`NORMAL`/`DEBUG`). `debug_file`, if given, is opened as an
    additional handler that mirrors everything the console handler shows --
    only meaningful (and only ever passed) when `level == DEBUG`.

    Independently of `level`, when debug mode is on (`appconfig.
    debug_enabled()`) the debug log (see this module's docstring) is
    attached too, at DEBUG severity; the logger's own threshold is the
    lower of the two, and each handler filters to its own.

    Safe to call more than once (handlers are cleared and rebuilt each
    time), which matters for tests that configure the logger repeatedly.

    Args:
        level (str): One of `SILENT`, `NORMAL`, `DEBUG`.
        debug_file (str | None): Path to also write every console line to.
        console (bool): False leaves the console handler off -- for the
            web interface, where stdout is the HTTP response (CGI) or not
            ours to write to (mod_wsgi).
    """
    for handler in list(_logger.handlers):
        _logger.removeHandler(handler)
        if handler is not _console_handler and handler is not _debug_log_handler:
            handler.close()

    console_level = _LEVELS[level]
    formatter = logging.Formatter(_TIMESTAMPED_FORMAT if level == DEBUG else _PLAIN_FORMAT)
    _console_handler.setFormatter(formatter)
    _console_handler.setLevel(console_level)
    if console:
        _logger.addHandler(_console_handler)
    else:
        # Without any handler, logging's "last resort" would print errors
        # to stderr anyway; the web interface already reports its own.
        _logger.addHandler(logging.NullHandler())

    if debug_file:
        file_handler = logging.FileHandler(debug_file)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(console_level)
        _logger.addHandler(file_handler)

    was_active = _debug_log_handler is not None
    debug_log = _configure_debug_log()
    _logger.setLevel(logging.DEBUG if debug_log else console_level)
    if debug_log and not was_active:
        _logger.debug("Debug log opened (planetGen %s, Python %s, argv=%r, cwd=%s, uid=%s)",
                      _package_version(), sys.version.split()[0], _redacted_argv(), os.getcwd(),
                      os.getuid() if hasattr(os, "getuid") else "n/a")


def _redacted_argv():
    """`sys.argv` with any password option's value (`--mysql-password X`,
    `--mysql-password=X`) withheld."""
    shown, hide_next = [], False
    for arg in sys.argv:
        if hide_next:
            shown.append("<withheld>")
            hide_next = False
        elif arg.startswith("-") and "password" in arg.lower():
            if "=" in arg:
                shown.append(arg.split("=", 1)[0] + "=<withheld>")
            else:
                shown.append(arg)
                hide_next = True
        else:
            shown.append(arg)
    return shown


def _package_version():
    try:
        from stellarObjects._version import __version__
        return __version__
    except Exception:  # noqa: BLE001
        return "unknown"


def _configure_debug_log():
    """
    Attaches (or, with debug mode off, detaches) the debug log and the
    extras that go with it. Never raises: a broken `config.json` or a log
    file this process can't write to just leaves the debug log off, with
    one warning on stderr (the Apache error log, for the web interface).

    Returns:
        bool: Whether the debug log is attached.
    """
    global _debug_log_handler
    try:
        config = appconfig.load_config()
        enabled = appconfig.debug_enabled(config)
        path = appconfig.log_file_path(config)
    except Exception as exc:  # noqa: BLE001 -- logging must never stop the program
        print(f"planetgen: debug log disabled, could not read config.json: {exc}", file=sys.stderr)
        enabled, path = False, None

    handler = _debug_log_handler
    if handler is not None and (not enabled or handler.baseFilename != os.path.abspath(path)):
        _detach_debug_log()
        handler = None

    if not enabled:
        return False

    if handler is None:
        try:
            # WatchedFileHandler reopens the file when logrotate moves it
            # away, so a long-running process (the API, a galaxy run)
            # follows the rotation instead of writing into the old file.
            handler = logging.handlers.WatchedFileHandler(path, encoding="utf-8")
        except OSError as exc:
            if path not in _warned_unwritable:
                _warned_unwritable.add(path)
                print(f"planetgen: debug is on but the debug log {path} can't be opened ({exc}); "
                      f"run install.sh/update.sh as root to create it.", file=sys.stderr)
            return False
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(_DebugLogFormatter(_DEBUG_LOG_FORMAT))
        handler.addFilter(_OriginFilter())
        _debug_log_handler = handler
        _attach_extras(handler)

    _logger.addHandler(handler)
    if handler not in _trace_logger.handlers:
        _trace_logger.addHandler(handler)
    return True


def _attach_extras(handler):
    """Routes other libraries' logging, uncaught exceptions and warnings to
    the debug log, and starts tracing random draws."""
    global _root_stderr_handler
    root = logging.getLogger()
    if not root.handlers:
        # Adding a handler to the root logger switches off logging's
        # built-in "last resort" stderr output for WARNING and up; keep
        # that output going (it's how library warnings reach Apache's
        # error log today).
        _root_stderr_handler = logging.StreamHandler(sys.stderr)
        _root_stderr_handler.setLevel(logging.WARNING)
        root.addHandler(_root_stderr_handler)
    root.addHandler(handler)
    root.setLevel(logging.DEBUG)
    logging.captureWarnings(True)
    _install_excepthooks()
    _trace_random(True)


def _detach_debug_log():
    global _debug_log_handler, _root_stderr_handler
    handler = _debug_log_handler
    if handler is None:
        return
    _logger.removeHandler(handler)
    _trace_logger.removeHandler(handler)
    root = logging.getLogger()
    root.removeHandler(handler)
    if _root_stderr_handler is not None:
        root.removeHandler(_root_stderr_handler)
        _root_stderr_handler = None
    root.setLevel(logging.WARNING)
    logging.captureWarnings(False)
    _trace_random(False)
    handler.close()
    _debug_log_handler = None


def debug_log_active():
    """Whether the debug log is attached -- lets a caller skip building
    expensive debug-only detail when nothing would record it."""
    return _debug_log_handler is not None


_previous_excepthook = None
_previous_threading_excepthook = None


def _install_excepthooks():
    """Logs uncaught exceptions (main thread and other threads) to the
    debug log, then hands them to whatever hook was there before, so the
    usual traceback still reaches stderr."""
    global _previous_excepthook, _previous_threading_excepthook
    if _previous_excepthook is not None:
        return
    _previous_excepthook = sys.excepthook

    def _excepthook(exc_type, exc, tb):
        if _debug_log_handler is not None and not issubclass(exc_type, KeyboardInterrupt):
            _logger.critical("Uncaught exception", exc_info=(exc_type, exc, tb))
        _previous_excepthook(exc_type, exc, tb)

    sys.excepthook = _excepthook

    if hasattr(threading, "excepthook"):
        _previous_threading_excepthook = threading.excepthook

        def _threading_excepthook(hook_args):
            if _debug_log_handler is not None:
                _logger.critical("Uncaught exception in thread %s",
                                 getattr(hook_args.thread, "name", "?"),
                                 exc_info=(hook_args.exc_type, hook_args.exc_value, hook_args.exc_traceback))
            _previous_threading_excepthook(hook_args)

        threading.excepthook = _threading_excepthook


# -- Random-draw tracing ------------------------------------------------------

_RANDOM_FUNCTIONS = (
    "random", "uniform", "randint", "randrange", "choice", "choices", "sample", "shuffle",
    "gauss", "normalvariate", "lognormvariate", "triangular", "expovariate", "betavariate",
    "gammavariate", "paretovariate", "weibullvariate", "vonmisesvariate", "binomialvariate",
    "getrandbits",
)
_original_random_functions = {}

# planetGen's own code lives under the repo (a source checkout, including
# src/html/) or, once pip-installed, in the stellarObjects package and the
# top-level generate module.
_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_DIR = os.path.dirname(os.path.dirname(_PACKAGE_DIR))
_OWN_MODULES = ("generate", "__main__")

_REPR_LIMIT = 400


def _short_repr(value):
    text = repr(value)
    if len(text) > _REPR_LIMIT:
        text = text[:_REPR_LIMIT] + f"... ({len(text)} chars)"
    return text


def _is_own_frame(frame):
    filename = frame.f_code.co_filename
    if filename.startswith(_PACKAGE_DIR) or filename.startswith(_REPO_DIR + os.sep):
        return True
    return frame.f_globals.get("__name__") in _OWN_MODULES


_IDENTIFIER = re.compile(r"(?<![\w.'\"])([A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*)")
_ASSIGNMENT_TARGET = re.compile(r"^\s*([\w.\[\], ]+?)\s*(?:[-+*/]?=)(?!=)")
_MISSING = object()


def _is_simple_value(value):
    if isinstance(value, (bool, int, float)):
        return True
    if isinstance(value, (tuple, list)) and 0 < len(value) <= 8:
        return all(isinstance(v, (bool, int, float, str)) for v in value)
    return isinstance(value, str) and len(value) <= 40


def _referenced_values(frame, source):
    """
    The plain values (numbers, short strings, small tuples) of the names
    a roll's source line refers to, e.g. `FLAVOR_CHANCE_PLANET=0.15` for
    `if random.random() < FLAVOR_CHANCE_PLANET:` -- the probability or
    threshold the drawn number is about to be compared against. Names
    being assigned on that line are skipped (they still hold the old
    value).
    """
    assigned = set()
    target = _ASSIGNMENT_TARGET.match(source)
    if target:
        assigned = set(_IDENTIFIER.findall(target.group(1)))
    found = []
    seen = set()
    for dotted in _IDENTIFIER.findall(source):
        if dotted in seen or dotted in assigned or dotted.startswith("random.") or keyword.iskeyword(dotted):
            continue
        seen.add(dotted)
        head, *rest = dotted.split(".")
        value = frame.f_locals.get(head, _MISSING)
        if value is _MISSING:
            value = frame.f_globals.get(head, _MISSING)
        for attr in rest:
            if value is _MISSING:
                break
            try:
                value = getattr(value, attr, _MISSING)
            except Exception:  # noqa: BLE001 -- a property may raise
                value = _MISSING
        if value is not _MISSING and _is_simple_value(value):
            found.append(f"{dotted}={_short_repr(value)}")
    return found


def _make_tracer(name, original):
    @functools.wraps(original)
    def traced(*args, **kwargs):
        result = original(*args, **kwargs)
        frame = sys._getframe(1)
        if _is_own_frame(frame):
            arg_text = ", ".join([_short_repr(a) for a in args]
                                 + [f"{k}={_short_repr(v)}" for k, v in kwargs.items()])
            shown = args[0] if name == "shuffle" and args else result
            source = linecache.getline(frame.f_code.co_filename, frame.f_lineno).strip()
            values = _referenced_values(frame, source)
            _trace_logger.debug("roll random.%s(%s) -> %s  <- %s%s", name, arg_text, _short_repr(shown), source,
                          f"  [{', '.join(values)}]" if values else "", stacklevel=2)
        return result

    traced._planetgen_traced = True
    return traced


def _trace_random(enable):
    """
    Wraps (or unwraps) the global `random` module's functions so each call
    from planetGen's own code is logged at DEBUG severity. Only the
    module-level functions are wrapped: they're what every generator
    module calls (`random.uniform(...)`), while the methods they use
    internally belong to the hidden `random.Random` instance and stay
    untouched -- so nothing is logged twice and no extra numbers are drawn.
    """
    if enable:
        for name in _RANDOM_FUNCTIONS:
            original = getattr(random, name, None)
            if original is None or getattr(original, "_planetgen_traced", False):
                continue
            _original_random_functions[name] = original
            setattr(random, name, _make_tracer(name, original))
    else:
        for name, original in _original_random_functions.items():
            setattr(random, name, original)
        _original_random_functions.clear()


def set_console(console):
    """Routes console output through `console.print(...)` (a live `rich.progress.Progress`'s own console)."""
    _console_handler.set_console(console)


def reset_console():
    """Reverts console output to a plain stdout write."""
    _console_handler.set_console(None)


# stacklevel=2 on each wrapper so the debug log's file:line/function is the
# caller's, not this module's.
def debug(message, *args, **kwargs):
    kwargs.setdefault("stacklevel", 2)
    _logger.debug(message, *args, **kwargs)


def normal(message, *args, **kwargs):
    kwargs.setdefault("stacklevel", 2)
    _logger.info(message, *args, **kwargs)


def error(message, *args, **kwargs):
    kwargs.setdefault("stacklevel", 2)
    _logger.error(message, *args, **kwargs)


def trace(message, *args, **kwargs):
    """Debug-log-only detail (never shown on the console, even with
    `--debug`); a no-op unless the debug log is on."""
    if _debug_log_handler is not None:
        kwargs.setdefault("stacklevel", 2)
        _trace_logger.debug(message, *args, **kwargs)


def exception(message, *args, **kwargs):
    """`error`, plus the traceback of the exception being handled."""
    kwargs.setdefault("stacklevel", 2)
    _logger.error(message, *args, exc_info=True, **kwargs)


def choice(description, chosen, reason):
    """Debug-only convenience for narrating a randomized/algorithmic decision: what was chosen and why."""
    _logger.debug(f"{description}: chose {chosen!r} ({reason})", stacklevel=2)


@contextmanager
def timed_phase(label):
    """
    Debug-only convenience for timing one named phase of generation:
    logs `"{label}: {elapsed}ms"` at DEBUG severity when the `with` block
    exits, so `--debug`'s own timestamped output doubles as a per-phase
    profile with no separate instrumentation needed -- a benchmark can
    either read the timestamps directly or just parse the elapsed value
    already printed in the message.

    Formats the elapsed time whether or not DEBUG is enabled (matching
    every other call in this module, which never checks the level
    itself) -- a `time.perf_counter()` delta and an f-string are cheap
    enough that gating on `isEnabledFor` isn't worth the extra branch.

    Args:
        label (str): A short description of the phase being timed.
    """
    start = time.perf_counter()
    try:
        yield
    finally:
        elapsed_ms = (time.perf_counter() - start) * 1000
        # stacklevel 3: past contextlib's __exit__ to the `with` statement.
        _logger.debug(f"{label}: {elapsed_ms:.3f}ms", stacklevel=3)


# Sensible default so anything that imports stellarObjects modules directly
# (tests, `python -c`, etc.) without going through `generate.py`'s `main()`
# still gets today's plain, untimed output rather than silence.
configure(NORMAL)
