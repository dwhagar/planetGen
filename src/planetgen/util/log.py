"""
Unified CLI Logging
====================

A single, shared logging channel for `planetgen` and every planetGen
module -- everything routes through here instead of bare
`print()` calls, so `--debug`/`--quiet`/`--silent` behave consistently
everywhere.

Three severities, matched one-to-one to stdlib `logging` levels so
threshold filtering is just `Logger.setLevel`, not reinvented:

    SILENT  (logging.ERROR)    -- only `error(...)` is shown.
    NORMAL  (logging.INFO)     -- the default; matches the plain, untimed
                                   text this program has always printed.
    DEBUG   (logging.DEBUG)    -- everything, always timestamped, and
                                   optionally mirrored to a file.

`configure()` is called exactly once, from `planetgen`'s `main()`, right
after argument parsing. Every other module just does
`from planetgen.util import log` and calls `log.debug(...)`/`log.normal(...)`/
`log.error(...)`.

**Debug mode** is `--debug` on the command line, or `config.json`'s
`"debug": true` (or `PLANETGEN_DEBUG=1`), which makes every `planetgen`
command a debug run. A debug run writes every decision it makes to
**stdout**, one line each, in one format:

    2026-10-23T17:29:42.729-0400 DEBUG    worker.py[2142] draw.random() ... | roll draw.random() -> 0.31

(ISO time to the millisecond with the UTC offset, severity, the process
and, for a thread other than the main one, its name, the logger, the
source line that logged it, and the message). The lines of the worker
processes that do a run's work reach the command's stdout too, as they
are written: a unit that is running in a work-queue worker (whose own
stdout goes nowhere) relays its lines through Redis (`relayed`,
`read_relayed`) to the process that started the run, which prints them
(`planetgen.queue.units`, `planetgen.queue.wrapper`). A debug run writes
to stderr instead when the command's own result goes to stdout
(`--output -`).

What a debug run writes is split into areas (`AREAS`), so that a galaxy
run's volume can be cut down: `--debug-areas draw,queue`, `"debug_areas"`
in `config.json` or `PLANETGEN_DEBUG_AREAS` pick the areas shown, and a
name with a leading `-` leaves one out (`-draw`). Errors and warnings
are always shown, and so are the normal (INFO) lines.

Separately from the console, `config.json`'s `"debug": true` (or
`PLANETGEN_DEBUG=1`) turns on the **debug log**: every entry point --
`planetgen`, the maintenance scripts, the `html/` CGI pages and the
Flask API -- appends everything at DEBUG severity to `/var/log/
planetgen.log` (`"log_file"`/`PLANETGEN_LOG_FILE` to move it), whatever
the console's own level is, in the same line format. Besides the
narrated decisions, the debug log (and a debug run's stdout) gets:

  - every draw from `planetgen.util.draw` made by planetGen's own
    code (`_trace_random`): the call, its arguments (weights included),
    the result, and the source line that asked for it -- so a roll like
    `if draw.random() < FLAVOR_CHANCE_PLANET:` shows both the number
    drawn and the probability it was compared against; draws on a
    `draw.Stream` a generator holds are traced too, a bulk draw
    (`random_array`, `randbits_list`) as one line;
  - every SQL statement `store.Connection` runs, with its timing;
  - the run's events: stages, units, workers, the machine's lease;
  - cache hits and misses;
  - uncaught exceptions (`sys.excepthook`/`threading.excepthook`),
    Python warnings, and anything other libraries log (Flask, urllib3, ...).

The tracer only wraps calls, it never draws extra numbers, so a seeded run
produces the same result with debug on or off. With debug off none of this
is installed and nothing is written.

`set_console`/`reset_console` let a caller with a live `rich.progress.Progress`
bar redirect this module's console output through `progress.console.print(...)`
instead of a raw stdout write, which is required while that bar is live (see
`planetgen`'s `_generation_progress()` docstring) -- a plain stdout write
fights the bar's own redraw.
"""

import functools
import keyword
import linecache
import logging
import logging.handlers
import os
import re
import sys
import threading
import time
from contextlib import contextmanager

from planetgen.util import draw, logpaths

SILENT = "silent"
NORMAL = "normal"
DEBUG = "debug"

_LEVELS = {
    SILENT: logging.ERROR,
    NORMAL: logging.INFO,
    DEBUG: logging.DEBUG,
}

_PLAIN_FORMAT = "%(message)s"

AREAS = {
    "decision": "narrated decisions and choices, stage timings",
    "draw": "every random draw: the call, its arguments, the result, the source line and the probability it met",
    "sql": "every SQL statement, with its timing",
    "queue": "runs, stages, units, workers and the machine's lease",
    "cache": "cache hits, misses and rebuilds",
    "web": "web requests, API calls and authorization decisions",
    "lib": "what other libraries log (Flask, RQ, Redis, urllib3)",
}
"""dict: The areas of debug output (`--debug-areas`), each with what it holds."""

_logger = logging.getLogger("planetgen")
_logger.propagate = False

# High-volume detail (every random draw, every SQL statement). It goes to
# the same places as everything else a debug run writes (stdout, the debug
# log, a worker's relay); only its area (`draw`, `sql`, ...) tells it apart.
_trace_logger = logging.getLogger("planetgen.trace")
_trace_logger.setLevel(logging.DEBUG)

_component = None
"""str | None: Short name for this process in the debug log (`set_component`);
defaults to the running script's file name."""

_debug_log_handler = None
_relay_handler = None
_console_debug = False
_tracing = False
_root_stderr_handler = None
_warned_unwritable = set()

_areas = None
"""frozenset | None: The areas debug output is limited to (`set_areas`); `None` is every area."""

_areas_spec = ""
"""str: The text `_areas` came from, as `--debug-areas` takes it; empty is every area."""


def parse_areas(spec):
    """
    The areas a `--debug-areas` text names: names separated by commas or spaces, `all` (or `*`) for every area, a
    leading `-` to leave one out (`-draw`; with only exclusions, every other area is shown). Empty is every area.

    Returns:
        frozenset | None: The areas shown, or `None` for all of them.

    Raises:
        ValueError: A name that is not an area.
    """
    names = [part for part in re.split(r"[,\s]+", spec or "") if part]
    chosen, left_out = set(), set()
    for name in names:
        if name in ("all", "*"):
            chosen.update(AREAS)
            continue
        target = left_out if name.startswith("-") else chosen
        bare = name.lstrip("-")
        if bare not in AREAS:
            raise ValueError(f"{bare!r} is not a debug area; the areas are {', '.join(AREAS)}")
        target.add(bare)
    if not names or (not chosen and not left_out):
        return None
    shown = (chosen or set(AREAS)) - left_out
    return None if shown == set(AREAS) else frozenset(shown)


def set_areas(spec):
    """Limits a debug run's output to the areas `spec` names (`parse_areas`); empty shows them all."""
    global _areas, _areas_spec
    _areas = parse_areas(spec)
    _areas_spec = (spec or "").strip()


def areas_spec():
    """The text of the areas in force (`""` for all), as `set_areas` takes it."""
    return _areas_spec


def area_enabled(area):
    """Whether output of `area` is shown: lets a caller skip building a line nobody would see."""
    return _areas is None or area in _areas


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


_WEB_DIRS = tuple(os.sep + "planetgen" + os.sep + name + os.sep for name in ("web", "api"))
_QUEUE_DIR = os.sep + "planetgen" + os.sep + "queue" + os.sep


def _area_of(record):
    """The area of a record that did not name one: another library's is `lib`; planetGen's own is `web`
    or `queue` by the folder it was logged from, else `decision`."""
    if not record.name.startswith("planetgen"):
        return "lib"
    if any(part in record.pathname for part in _WEB_DIRS):
        return "web"
    if _QUEUE_DIR in record.pathname:
        return "queue"
    return "decision"


class _AreaFilter(logging.Filter):
    """Stamps each record with its `area` and drops a DEBUG record of an area `--debug-areas` left out. Normal (INFO)
    lines, warnings and errors always pass: a filter that hid an error would hide the one line that matters."""

    def filter(self, record):
        area = getattr(record, "area", None)
        if area is None:
            area = record.area = _area_of(record)
        if record.levelno != logging.DEBUG:
            return True
        return _areas is None or area in _areas


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
_origin_filter = _OriginFilter()
_area_filter = _AreaFilter()


def _debug_formatter():
    return _DebugLogFormatter(_DEBUG_LOG_FORMAT)


class _ConsoleHandler(logging.Handler):
    """
    Writes to stdout by default; `set_console`/`reset_console` retarget it
    to a `rich.Console` (a live `Progress`'s own `console`) so output stays
    redraw-safe while a progress bar is on screen. A debug run's lines are
    flushed as they are written, so they can be watched as they happen even
    when stdout is a pipe.
    """

    def __init__(self):
        super().__init__()
        self._console = None
        self.to_stderr = False

    def set_console(self, console):
        self._console = console

    def emit(self, record):
        self.write(self.format(record))

    def write(self, message):
        """One finished line to the console (a line a worker formatted comes through here too)."""
        if self._console is not None:
            # ADM.23: one line, wrapped by whatever shows it (a terminal,
            # the web job log), never hard-broken at the console's width
            # (80 when piped); and brackets are text, not rich markup.
            self._console.print(message, soft_wrap=True, markup=False)
        else:
            stream = sys.stderr if self.to_stderr else sys.stdout
            print(message, file=stream)
        try:
            (sys.stderr if self.to_stderr else sys.stdout).flush()
        except (OSError, ValueError):
            pass


_console_handler = _ConsoleHandler()


def configure(level=NORMAL, debug_file=None, console=True, to_stderr=False, areas=None):
    """
    (Re)configures the shared logger for the given console severity `level`
    (`SILENT`/`NORMAL`/`DEBUG`). `debug_file`, if given, is opened as an
    additional handler that mirrors everything the console handler shows --
    only meaningful (and only ever passed) when `level == DEBUG`.

    Independently of `level`, when debug mode is on (`logpaths.
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
        to_stderr (bool): Write the console's lines to stderr, for a command whose own result goes to stdout.
        areas (str | None): The areas debug output is limited to (`parse_areas`); `None` takes `PLANETGEN_DEBUG_AREAS`,
            then `config.json`'s `"debug_areas"`, then every area.
    """
    global _console_debug
    for handler in list(_logger.handlers):
        _logger.removeHandler(handler)
        if (handler is not _console_handler and handler is not _debug_log_handler
                and handler is not _relay_handler):
            handler.close()

    set_areas(areas if areas is not None else _configured_areas())
    console_level = _LEVELS[level]
    formatter = _debug_formatter() if level == DEBUG else logging.Formatter(_PLAIN_FORMAT)
    _console_handler.setFormatter(formatter)
    _console_handler.setLevel(console_level)
    _console_handler.to_stderr = to_stderr
    _console_handler.filters.clear()
    if level == DEBUG:
        _console_handler.addFilter(_origin_filter)
        _console_handler.addFilter(_area_filter)
    if console:
        _logger.addHandler(_console_handler)
    else:
        # Without any handler, logging's "last resort" would print errors
        # to stderr anyway; the web interface already reports its own.
        _logger.addHandler(logging.NullHandler())
    _console_debug = bool(console and level == DEBUG)

    if debug_file:
        file_handler = logging.FileHandler(debug_file)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(console_level)
        if level == DEBUG:
            file_handler.addFilter(_origin_filter)
            file_handler.addFilter(_area_filter)
        _logger.addHandler(file_handler)

    if _relay_handler is not None:
        _logger.addHandler(_relay_handler)
    was_active = _debug_log_handler is not None
    debug_log = _configure_debug_log()
    _logger.setLevel(logging.DEBUG if (debug_log or _relay_handler is not None) else console_level)
    _sync_tracing()
    if debug_log and not was_active:
        _logger.debug("Debug log opened (planetGen %s, Python %s, argv=%r, cwd=%s, uid=%s)",
                      _package_version(), sys.version.split()[0], _redacted_argv(), os.getcwd(),
                      os.getuid() if hasattr(os, "getuid") else "n/a")


def _configured_areas():
    """The areas `PLANETGEN_DEBUG_AREAS` or `config.json` name; every area when neither does or the text is wrong
    (one warning on stderr: logging must never stop the program)."""
    try:
        spec = logpaths.debug_areas()
        parse_areas(spec)
        return spec
    except Exception as exc:  # noqa: BLE001
        print(f"planetgen: debug_areas ignored, every area is shown: {exc}", file=sys.stderr)
        return ""


def config_debug():
    """Whether debug mode is on in `config.json` or the environment (`False` when the file cannot be read)."""
    try:
        return logpaths.debug_enabled()
    except Exception:  # noqa: BLE001 -- logging must never stop the program
        return False


def console_debug():
    """Whether this process writes a debug run's lines to its console: what a run hands to its workers."""
    return _console_debug


def add_debug_arguments(parser):
    """Adds `--debug` (and `--debug-areas`) to a standalone command's `parser` and to each of its subcommands, so
    every command in the project takes it. Call it once the subcommands are added; `apply_debug_arguments` reads it."""
    import argparse

    def add(target, default):
        target.add_argument("--debug", action="store_true", default=default,
                            help="Write every decision, random draw, SQL statement and queue event, with ISO "
                                 "timestamps, to stdout (and to the debug log when config.json turns debug on).")
        target.add_argument("--debug-areas", metavar="AREAS", default=None if default is False else default,
                            help="Limit --debug to some areas: " + ", ".join(AREAS) + " (a leading - leaves one out).")

    add(parser, False)
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for sub in action.choices.values():
                add(sub, argparse.SUPPRESS)        # a value given before the subcommand is not reset by the subparser


def apply_debug_arguments(args, parser=None):
    """
    Configures logging for a standalone command from the options `add_debug_arguments` added, and from debug mode in
    `config.json`. Call it right after parsing. Debug output goes to stdout.
    """
    wanted = getattr(args, "debug", False) or config_debug()
    areas = getattr(args, "debug_areas", None)
    if areas is not None:
        try:
            parse_areas(areas)
        except ValueError as exc:
            if parser is None:
                raise
            parser.error(f"--debug-areas: {exc}")
    configure(DEBUG if wanted else NORMAL, areas=areas)


_SECRET_OPTIONS = ("--mysql-password",)
"""Every command-line option in this project whose value is a secret."""

_SECRET_OPTION_WORDS = ("password", "passwd", "token", "secret")


def _is_secret_option(arg):
    """Whether `arg` (an option, possibly `--opt=value`) names a secret:
    any option mentioning a secret word, or any abbreviation argparse
    would accept for one of `_SECRET_OPTIONS` (`allow_abbrev` is on by
    default, so `--mysql-pass X` sets `--mysql-password` too). Prefixes
    shorter than `--mysql-p` are left alone: they're ambiguous with the
    other `--mysql-*` options, so argparse rejects them anyway."""
    if not arg.startswith("-"):
        return False
    option = arg.split("=", 1)[0].lower()
    if any(word in option for word in _SECRET_OPTION_WORDS):
        return True
    return option.startswith("--") and any(
        secret.startswith(option) and len(option) > len("--mysql-") for secret in _SECRET_OPTIONS)


def _redacted_argv():
    """`sys.argv` with any secret option's value (`--mysql-password X`,
    `--mysql-password=X`, or an abbreviation such as `--mysql-pass X`)
    withheld."""
    shown, hide_next = [], False
    for arg in sys.argv:
        if hide_next:
            shown.append("<withheld>")
            hide_next = False
        elif _is_secret_option(arg):
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
        from planetgen._version import __version__
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
        config = logpaths.read_config_file()
        enabled = logpaths.debug_enabled(config)
        path = logpaths.log_file_path(config)
        if not isinstance(path, str):
            # A non-string "log_file" in config.json would otherwise raise
            # a TypeError from the file handler below, past the OSError
            # handling there.
            raise TypeError(f'"log_file" must be a string, not {type(path).__name__}')
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
        handler.setFormatter(_debug_formatter())
        handler.addFilter(_origin_filter)
        handler.addFilter(_area_filter)
        _debug_log_handler = handler
        _attach_extras(handler)

    _logger.addHandler(handler)
    return True


def _attach_extras(handler):
    """Routes other libraries' logging, uncaught exceptions and warnings to
    the debug log."""
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


def _detach_debug_log():
    global _debug_log_handler, _root_stderr_handler
    handler = _debug_log_handler
    if handler is None:
        return
    _logger.removeHandler(handler)
    root = logging.getLogger()
    root.removeHandler(handler)
    if _root_stderr_handler is not None:
        root.removeHandler(_root_stderr_handler)
        _root_stderr_handler = None
    root.setLevel(logging.WARNING)
    logging.captureWarnings(False)
    handler.close()
    _debug_log_handler = None


def debug_log_active():
    """Whether anything records debug detail (the debug log, a debug run's console, a worker's relay) -- lets a
    caller skip building expensive debug-only detail when nothing would record it."""
    return _tracing


def _sync_tracing():
    """Starts tracing random draws while anything records debug detail, and stops when nothing does."""
    global _tracing
    _tracing = _debug_log_handler is not None or _console_debug or _relay_handler is not None
    _trace_random(_tracing)
    if _tracing:
        _install_excepthooks()


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
        if _tracing and not issubclass(exc_type, KeyboardInterrupt):
            _logger.critical("Uncaught exception", exc_info=(exc_type, exc, tb))
        _previous_excepthook(exc_type, exc, tb)

    sys.excepthook = _excepthook

    if hasattr(threading, "excepthook"):
        _previous_threading_excepthook = threading.excepthook

        def _threading_excepthook(hook_args):
            if _tracing:
                _logger.critical("Uncaught exception in thread %s",
                                 getattr(hook_args.thread, "name", "?"),
                                 exc_info=(hook_args.exc_type, hook_args.exc_value, hook_args.exc_traceback))
            _previous_threading_excepthook(hook_args)

        threading.excepthook = _threading_excepthook


# -- Random-draw tracing ------------------------------------------------------

_RANDOM_FUNCTIONS = (
    "random", "uniform", "randint", "randrange", "choice", "choices", "sample", "shuffle",
    "gauss", "normal", "getrandbits",
)
_original_random_functions = {}

# planetGen's own code lives under the repo (a source checkout, including
# src/html/) or, once pip-installed, in the planetgen package and the
# top-level generate module. This file is src/planetgen/util/log.py.
_PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
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
    `if draw.random() < FLAVOR_CHANCE_PLANET:` -- the probability or
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
        if dotted in seen or dotted in assigned or dotted.startswith(("random.", "draw.")) or keyword.iskeyword(dotted):
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


_STREAM_METHODS = ("random", "uniform", "randrange", "randint", "choice", "choices", "shuffle", "sample", "gauss",
                   "normal", "gamma", "dirichlet", "getrandbits", "random_array", "randbits_list")
"""tuple: Every draw a `draw.Stream` makes. The module-level functions (`draw.uniform`) pass straight to these, so
tracing the class traces them, and the draws on a stream a generator holds itself, in one place."""

_draw_state = threading.local()
_DRAW_STEM = os.path.splitext(os.path.abspath(draw.__file__))[0]
_SKIPPED_FILES = (_DRAW_STEM, os.path.join(os.path.dirname(os.path.abspath(__file__)), "random"))


def _is_wrapper_frame(frame):
    """Whether `frame` is the draw module or a sampling helper (`planetgen.util.random`) that draws for its caller:
    the line worth showing is the one that called them."""
    return os.path.splitext(os.path.abspath(frame.f_code.co_filename))[0] in _SKIPPED_FILES


def _describe_result(value):
    if hasattr(value, "shape") and hasattr(value, "tolist"):
        size = getattr(value, "size", 0)
        return f"<array shape={value.shape}: {_short_repr(value.ravel()[:3].tolist())}{' ...' if size > 3 else ''}>"
    if isinstance(value, list) and len(value) > 12:
        return f"<list of {len(value)}: {_short_repr(value[:3])[:-1]}, ...]>"
    return _short_repr(value)


def _make_tracer(name, original):
    @functools.wraps(original)
    def traced(self, *args, **kwargs):
        if getattr(_draw_state, "inside", False):
            return original(self, *args, **kwargs)     # a draw built on other draws (choice on randrange): one line
        _draw_state.inside = True
        try:
            result = original(self, *args, **kwargs)
        finally:
            _draw_state.inside = False
        if _areas is not None and "draw" not in _areas:
            return result
        frame = sys._getframe(1)
        skipped, helpers = 0, []
        while frame is not None and _is_wrapper_frame(frame):
            if os.path.splitext(os.path.abspath(frame.f_code.co_filename))[0] != _DRAW_STEM:
                helpers.append(frame.f_code.co_name)
            frame, skipped = frame.f_back, skipped + 1
        if frame is not None and _is_own_frame(frame):
            if name == "__init__":
                text = f"new stream seeded with {_short_repr(args[0] if args else kwargs.get('seed'))}"
                source = ""
            else:
                arg_text = ", ".join([_short_repr(a) for a in args]
                                     + [f"{k}={_short_repr(v)}" for k, v in kwargs.items()])
                shown = args[0] if name == "shuffle" and args else result
                text = f"roll draw.{name}({arg_text}) -> {_describe_result(shown)}"
                source = None
            if source is None:
                source = linecache.getline(frame.f_code.co_filename, frame.f_lineno).strip()
                values = _referenced_values(frame, source)
                if helpers:
                    text += f"  (via {', '.join(helpers)})"
                text += f"  <- {source}" + (f"  [{', '.join(values)}]" if values else "")
            _trace_logger.debug("%s", text, stacklevel=2 + skipped, extra={"area": "draw"})
        return result

    traced._planetgen_traced = True
    return traced


def _trace_random(enable):
    """
    Wraps (or unwraps) the draws of `planetgen.util.draw.Stream`, the one
    place every random number comes from, so each draw from planetGen's own
    code is logged at DEBUG severity in the `draw` area. A draw built on
    others (`choice` on `randrange`, `gamma` on `gauss`) is one line, a
    bulk draw (`random_array`) one line for all of it, and a draw made
    through a sampling helper names the helper and the line that called
    it. The wrappers only observe: no extra numbers are drawn, so a seeded
    run is the same with tracing on or off.
    """
    if enable:
        for name in ("__init__", *_STREAM_METHODS):
            original = getattr(draw.Stream, name, None)
            if original is None or getattr(original, "_planetgen_traced", False):
                continue
            _original_random_functions[name] = original
            setattr(draw.Stream, name, _make_tracer(name, original))
    else:
        for name, original in _original_random_functions.items():
            setattr(draw.Stream, name, original)
        _original_random_functions.clear()


# -- Worker lines to the starter's stdout -------------------------------------

RELAY_MAXLEN = 200_000
"""int: The most lines a run keeps in Redis for its starter to print (the oldest go first)."""

RELAY_KEEP_SECONDS = 3600
"""int: How long a run's relayed lines stay in Redis once nothing is adding to them."""

_RELAY_BATCH = 500
_RELAY_SECONDS = 0.2


class _RelayHandler(logging.Handler):
    """
    Collects a unit's formatted lines and sends them to the run's stream in Redis (`XADD`, a batch at a time, at
    least every `_RELAY_SECONDS`), where the process that started the run reads them (`read_relayed`) and prints
    them as they come. A thread flushes on the clock, so the line before a long stretch of work is not held back
    until the stretch ends.
    """

    def __init__(self, conn, key):
        super().__init__(logging.DEBUG)
        self.conn, self.key = conn, key
        self._lines = []
        self._guard = threading.Lock()
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, name="planetgen-relay", daemon=True)
        self._thread.start()

    def emit(self, record):
        line = self.format(record)
        with self._guard:
            self._lines.append(line)
            full = len(self._lines) >= _RELAY_BATCH
        if full:
            self.flush()

    def _run(self):
        while not self._stop.wait(_RELAY_SECONDS):
            self.flush()

    def flush(self):
        with self._guard:
            lines, self._lines = self._lines, []
        if not lines:
            return
        try:
            pipe = self.conn.pipeline(transaction=False)
            for line in lines:
                pipe.xadd(self.key, {"l": line}, maxlen=RELAY_MAXLEN, approximate=True)
            pipe.expire(self.key, RELAY_KEEP_SECONDS)
            pipe.execute()
        except Exception:  # noqa: BLE001 -- debug output must never fail the unit it describes
            pass

    def close(self):
        self._stop.set()
        self._thread.join(timeout=2)
        self.flush()
        super().close()


@contextmanager
def relayed(conn, key, areas=""):
    """
    For a unit of work in a queue worker, whose own stdout goes nowhere: while the block runs, everything it logs at
    DEBUG goes to the stream `key` in Redis (`conn`) for the run's starter to print (`read_relayed`), in the areas
    `areas` names (`parse_areas`; empty is every area). A worker that already writes the debug log keeps doing so.
    """
    global _relay_handler, _areas, _areas_spec
    handler = _RelayHandler(conn, key)
    handler.setFormatter(_debug_formatter())
    handler.addFilter(_origin_filter)
    handler.addFilter(_area_filter)
    before = (_areas, _areas_spec, _logger.level, _relay_handler)
    set_areas(areas)
    _relay_handler = handler
    _logger.addHandler(handler)
    _logger.setLevel(logging.DEBUG)
    _sync_tracing()
    try:
        yield handler
    finally:
        _logger.removeHandler(handler)
        handler.close()
        _relay_handler = before[3]
        _areas, _areas_spec = before[0], before[1]
        _logger.setLevel(before[2])
        _sync_tracing()


def flush_relayed():
    """Sends the lines a unit has written so far to the run's starter now (no-op outside `relayed`)."""
    if _relay_handler is not None:
        _relay_handler.flush()


def read_relayed(conn, key, after="0", count=2000):
    """The lines workers relayed to `key` after stream id `after`: `(lines, last id)`."""
    rows = conn.xrange(key, min="(" + after if after != "0" else "-", max="+", count=count)
    if not rows:
        return [], after
    text = lambda value: value.decode("utf-8", "replace") if isinstance(value, bytes) else value   # noqa: E731
    return [text(fields.get(b"l", fields.get("l", ""))) for _entry, fields in rows], text(rows[-1][0])


class RelayReader:
    """
    The starter's side of `relayed`: `drain()` prints, on this process's console (and in the `--debug FILE` it was
    given), the lines the run's workers have relayed since the last call. A run calls it on every poll and once more
    when it ends.
    """

    def __init__(self, conn, key):
        self.conn, self.key, self.after = conn, key, "0"

    def finish(self, quiet=0.4, limit=3.0):
        """The run is over: prints what is left, and keeps looking until no line has come for `quiet` seconds (the
        unit that ended the run is still writing its last ones), but no longer than `limit`."""
        deadline = time.monotonic() + limit
        last = time.monotonic()
        while True:
            if self.drain():
                last = time.monotonic()
            if time.monotonic() - last >= quiet or time.monotonic() >= deadline:
                return
            time.sleep(0.05)

    def drain(self):
        """Prints what is waiting; returns how many lines."""
        total = 0
        while True:
            try:
                lines, self.after = read_relayed(self.conn, self.key, self.after)
            except Exception:  # noqa: BLE001 -- a lost line is not a reason to stop the run
                return total
            if not lines:
                return total
            for line in lines:
                write_line(line)
            total += len(lines)


def write_line(text):
    """Writes a line another process already formatted to the console and to the `--debug FILE` mirror."""
    if _console_debug:
        _console_handler.write(text)
    for handler in _logger.handlers:
        if isinstance(handler, logging.FileHandler) and handler is not _debug_log_handler:
            try:
                handler.stream.write(text + "\n")
                handler.flush()
            except (OSError, ValueError):
                pass


def set_console(console):
    """Routes console output through `console.print(...)` (a live `rich.progress.Progress`'s own console)."""
    _console_handler.set_console(console)


def reset_console():
    """Reverts console output to a plain stdout write."""
    _console_handler.set_console(None)


# stacklevel=2 on each wrapper so the debug log's file:line/function is the
# caller's, not this module's.
def debug(message, *args, area=None, **kwargs):
    """A decision, narrated. `area` (one of `AREAS`) files it where `--debug-areas` can pick it; the default is by
    the folder of the module that logs it (`web`, `queue`, else `decision`)."""
    kwargs.setdefault("stacklevel", 2)
    if area is not None:
        kwargs["extra"] = {"area": area}
    _logger.debug(message, *args, **kwargs)


def normal(message, *args, **kwargs):
    kwargs.setdefault("stacklevel", 2)
    _logger.info(message, *args, **kwargs)


def error(message, *args, **kwargs):
    kwargs.setdefault("stacklevel", 2)
    _logger.error(message, *args, **kwargs)


def trace(message, *args, area="decision", **kwargs):
    """High-volume detail (an SQL statement, a cache hit), in `area` (one of `AREAS`). It goes where every other
    debug line goes (stdout, the debug log, a worker's relay); a no-op unless one of them is on, so the caller
    may build `message` only when `debug_log_active()` says something would record it."""
    if _tracing:
        kwargs.setdefault("stacklevel", 2)
        kwargs["extra"] = {"area": area}
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


# Sensible default so anything that imports planetgen modules directly
# (tests, `python -c`, etc.) without going through `planetgen`'s `main()`
# still gets today's plain, untimed output rather than silence.
configure(NORMAL)
