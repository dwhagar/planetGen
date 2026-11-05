"""
No code under `planetgen/` formats its own log or report lines (TEST.183, UX.102).

A line an admin or the user reads in a log, a report or a console stream goes
through the one logging layer (`planetgen.util.logrecord`, `util/log.py`,
`util/stagemark.py`, `queue/runrecord.py`), so it carries the same RFC 3339
timestamp, severity, component, fields, dividers and "Stage N of M" markers as
every other line (docs/design/operations/logging-standard.md, rule 5). This test
reads every module of the package as a syntax tree and fails when a call that
writes (`print`, `.write`, `.writelines`) or logs (`.debug`, `.info`, `.notice`,
`.warning`, `.error`, `.critical`, `.exception`, `.log`, `.normal`) builds its
text with

* a time call (`strftime`, `isoformat`, `ctime`, `asctime`, `now`, `utcnow`,
  `fromtimestamp`, `utcfromtimestamp`), directly, through a `:%H:%M` style
  format spec, or through a variable the same function set from one; or
* a divider: a string that starts or ends with three or more of `= - * # _ ~ +`,
  or one such character multiplied by a count.

The result is the set of places that format their own lines. Each must be in
`ALLOWED` under its file and function with the reason it is not a log line.
`ALLOWED` is empty of anything else: a new output path either goes through the
layer or is named here, with its reason, in the change that adds it.

Migrations (`db/migrations`) are frozen history and are not scanned.
"""
import ast
import os
import re

import pytest

PACKAGE = os.path.join(os.path.dirname(__file__), "..", "planetgen")

SKIPPED = ("db/migrations/",)
"""tuple: Path prefixes (relative to the package, with `/`) that are not scanned."""

TIME_CALLS = {"strftime", "isoformat", "ctime", "asctime", "now", "utcnow", "fromtimestamp", "utcfromtimestamp"}
"""set: Calls that turn a moment into text."""

WRITERS = {"print", "write", "writelines"}
"""set: Calls that put text on a screen or in a file."""

LOGGERS = {"debug", "info", "notice", "warning", "warn", "error", "critical", "exception", "log", "normal"}
"""set: Logger methods whose text is a log line."""

DIVIDER_CHARS = "=-*#_~+"
"""str: The characters a hand-made divider is drawn with."""

EDGE = re.compile(r"^\s*[=\-*#_~+]{3,}|[=\-*#_~+]{3,}\s*$")
"""re.Pattern: A string that starts or ends with a run of divider characters."""

TIME_SPEC = re.compile(r"%[YymdHIMSpfjaAbBzZcxX]")
"""re.Pattern: A time directive in a format spec (`:%H:%M`), not a percentage (`:.1%`)."""

ALLOWED = {
    ("cli/migrate.py", "print_initial_admin_login"): (
        "The first admin password is shown once on the console, boxed so it is not missed; it must never reach a log file, "
        "so it cannot go through the logging layer (an interactive display, not a record)."
    ),
    ("cli/lockouts.py", "main"): (
        "`planetgen lockouts` lists the locked subjects as a table for the person at the console; a data listing like "
        "`query`, `config` and `queue` (logging standard, rule 5), not a log line."
    ),
    ("cli/queue.py", "run_queue"): (
        "`planetgen queue` prints its status sentence for the person at the console; a data listing, not a log line."
    ),
}
"""dict: `(file, function)` to the reason that function's own formatting is not a log line."""


def _name(func):
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return None


def _has_time(node):
    for part in ast.walk(node):
        if isinstance(part, ast.Call) and _name(part.func) in TIME_CALLS:
            return True
        if isinstance(part, ast.FormattedValue) and part.format_spec is not None:
            spec = "".join(s.value for s in ast.walk(part.format_spec) if isinstance(s, ast.Constant) and isinstance(s.value, str))
            if TIME_SPEC.search(spec):
                return True
    return False


def _divider(node):
    for part in ast.walk(node):
        if isinstance(part, ast.Constant) and isinstance(part.value, str):
            if len(part.value.strip()) >= 3 and EDGE.search(part.value):
                return part.value
        if isinstance(part, ast.BinOp) and isinstance(part.op, ast.Mult):
            for text, count in ((part.left, part.right), (part.right, part.left)):
                if isinstance(text, ast.Constant) and isinstance(text.value, str) and len(text.value) == 1 and text.value in DIVIDER_CHARS:
                    return ast.unparse(part)
    return None


def find_own_formatting(source, filename="<source>"):
    """
    Returns `[(function, line, kind)]` for every writing or logging call in
    `source` that builds its text with a time call or a divider. `function` is
    the innermost function that holds the call, or `<module>`.

    Arguments:
        source (str): Python source text.
        filename (str): Its name, for syntax errors.

    Returns:
        list: Tuples of function name, line number and what was found.
    """
    tree = ast.parse(source, filename)
    found = []

    def visit(scope, stamped, scope_name):
        # Variables this scope set from a time call, so `when = ...strftime(...); print(when)` is caught.
        stamped = set(stamped)
        body = list(ast.iter_child_nodes(scope))
        stack = list(body)
        while stack:
            node = stack.pop()
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                visit(node, stamped, node.name)
                continue
            if isinstance(node, ast.Assign) and _has_time(node.value):
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        stamped.add(target.id)
            stack.extend(ast.iter_child_nodes(node))
        stack = list(body)
        while stack:
            node = stack.pop()
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            stack.extend(ast.iter_child_nodes(node))
            if not isinstance(node, ast.Call) or _name(node.func) not in WRITERS | LOGGERS:
                continue
            for arg in list(node.args) + [k.value for k in node.keywords]:
                if _has_time(arg):
                    found.append((scope_name, node.lineno, "time"))
                elif any(isinstance(n, ast.Name) and n.id in stamped for n in ast.walk(arg)):
                    found.append((scope_name, node.lineno, "time"))
                if _divider(arg):
                    found.append((scope_name, node.lineno, "divider"))

    visit(tree, set(), "<module>")
    return sorted(set(found), key=lambda hit: (hit[1], hit[0], hit[2]))


def _modules():
    for folder, _dirs, names in os.walk(PACKAGE):
        for name in sorted(names):
            if not name.endswith(".py"):
                continue
            path = os.path.join(folder, name)
            relative = os.path.relpath(path, PACKAGE).replace(os.sep, "/")
            if not relative.startswith(SKIPPED):
                yield relative, path


def _all_hits():
    hits = {}
    for relative, path in _modules():
        with open(path, encoding="utf-8") as handle:
            for function, line, kind in find_own_formatting(handle.read(), path):
                hits.setdefault((relative, function), []).append((line, kind))
    return hits


def test_no_module_formats_its_own_log_lines():
    stray = {key: where for key, where in _all_hits().items() if key not in ALLOWED}
    lines = [f"{path}:{line} in {function} ({kind})" for (path, function), where in sorted(stray.items()) for line, kind in where]
    assert not stray, (
        "These calls write or log a line they format themselves (their own timestamp or divider). Route them through "
        "planetgen.util.logrecord / util.log, or add the function to ALLOWED with the reason it is not a log line "
        "(docs/design/operations/logging-standard.md, rule 5):\n  " + "\n  ".join(lines)
    )


def test_every_allowed_exception_is_still_needed_and_has_a_reason():
    hits = _all_hits()
    for key, reason in ALLOWED.items():
        assert len(reason.split()) >= 8, f"{key} needs a real reason"
        assert key in hits, f"{key} no longer formats its own lines; remove it from ALLOWED"


@pytest.mark.parametrize("source, kinds", [
    ("import time\ndef f(log):\n    log.info(time.strftime('%H:%M'))\n", ["time"]),
    ("def f():\n    print(f'{when:%H:%M} done')\n", ["time"]),
    ("from datetime import datetime\ndef f(out):\n    stamp = datetime.now().isoformat()\n    out.write(stamp + ' x')\n", ["time"]),
    ("def f():\n    print('=' * 72)\n", ["divider"]),
    ("def f(log):\n    log.error('=== failed ===')\n", ["divider"]),
    ("def f(out):\n    out.write('----------')\n", ["divider"]),
    ("def f():\n    print('plain words', 3)\n", []),
    ("def f(x, log):\n    log.info(f'{x:.1%} done, {x:,} rows')\n", []),
    ("def f():\n    stamp = time.strftime('%H')\n    return stamp\n", []),
    ("def f():\n    return 'a' + '=' * 3\n", []),
    ("def f(log):\n    log.info('a - b')\n", []),
])
def test_the_scan_sees_what_it_should(source, kinds):
    assert [kind for _function, _line, kind in find_own_formatting(source)] == kinds


# --- one logging module (OPS.104) -------------------------------------------------------------------------------------

THE_MODULE = "util/log.py"
"""str: The one module that holds a formatter class, a handler class or a `setFormatter` of the logging layer."""

HANDLER_WORDS = ("Formatter", "Handler")

HANDLER_CALLS = {"Formatter", "StreamHandler", "FileHandler", "WatchedFileHandler", "RotatingFileHandler", "TimedRotatingFileHandler",
                 "SocketHandler", "SysLogHandler", "NullHandler", "MemoryHandler", "QueueHandler", "setFormatter", "getLogger"}
"""set: Calls that make a formatter or a handler of the standard library, give a handler its formatter or fetch a logger by name (a module logs through
`log.normal` / `log.get`, OPS.105)."""

OTHER_HANDLERS = {}
"""dict: `(file, what)` to the reason a formatter, handler or `getLogger` outside the module is allowed. It is empty since OPS.105 moved the
activity log, the queue run log and the direct users onto `util/log.py`; a new entry needs a reason a person can check."""


def find_handlers(source, filename="<source>"):
    """`[(what, line)]`: every class that derives from a formatter or handler, and every call that makes one or sets one, in `source`."""
    found = []
    for node in ast.walk(ast.parse(source, filename)):
        if isinstance(node, ast.ClassDef):
            for base in node.bases:
                name = _name(base) or ""
                if name.endswith(HANDLER_WORDS) and (ast.unparse(base).startswith("logging") or isinstance(base, ast.Name)):
                    found.append((f"class {node.name}", node.lineno))
                    break
        elif isinstance(node, ast.Call) and _name(node.func) in HANDLER_CALLS:
            found.append((_name(node.func), node.lineno))
    return sorted(set(found), key=lambda hit: hit[1])


def _handler_hits():
    hits = {}
    for relative, path in _modules():
        if relative == THE_MODULE:
            continue
        with open(path, encoding="utf-8") as handle:
            for what, line in find_handlers(handle.read(), path):
                hits.setdefault((relative, what), []).append(line)
    return hits


def test_no_formatter_class_or_handler_is_made_outside_the_logging_module():
    stray = {key: lines for key, lines in _handler_hits().items() if key not in OTHER_HANDLERS}
    assert not stray, (
        "A formatter or handler is made outside planetgen/util/log.py. The one module writes every line (logging standard, rule 5): "
        "add the form to log.LineFormatter or the handler to log.py, or name the straggler in OTHER_HANDLERS with its reason:\n  "
        + "\n  ".join(f"{path}:{lines[0]} {what}" for (path, what), lines in sorted(stray.items())))


def test_every_straggler_is_still_there_and_has_a_reason():
    hits = _handler_hits()
    for key, reason in OTHER_HANDLERS.items():
        assert len(reason.split()) >= 8, f"{key} needs a real reason"
        assert key in hits, f"{key} is gone; remove it from OTHER_HANDLERS"


def test_the_logging_module_has_one_formatter_class_and_the_record_layer_has_none():
    with open(os.path.join(PACKAGE, "util", "log.py"), encoding="utf-8") as handle:
        formatters = [node.name for node in ast.walk(ast.parse(handle.read())) if isinstance(node, ast.ClassDef)
                      and any((_name(base) or "").endswith("Formatter") for base in node.bases)]
    assert formatters == ["LineFormatter"]
    with open(os.path.join(PACKAGE, "util", "logrecord.py"), encoding="utf-8") as handle:
        assert not find_handlers(handle.read())


@pytest.mark.parametrize("source, found", [
    ("import logging\nclass F(logging.Formatter):\n    pass\n", ["class F"]),
    ("import logging\nh = logging.StreamHandler()\nh.setFormatter(x)\n", ["StreamHandler", "setFormatter"]),
    ("import logging\nclass H(logging.handlers.WatchedFileHandler):\n    pass\n", ["class H"]),
    ("import logging\nlogging.getLogger('x').addHandler(h)\n", ["getLogger"]),
])
def test_the_scan_finds_a_formatter_a_handler_and_what_sets_one(source, found):
    assert [what for what, _line in find_handlers(source)] == found


# --- no print for a log line (OPS.105) --------------------------------------------------------------------------------

PRINT_ALLOWED = {
    "_version.py": "Prints the banner a person asked for with the version command: an answer, not a record.",
    "cli/accounts.py": "Prints who is the Owner after the command ran: the answer the admin asked for.",
    "cli/backup.py": "Lists, makes, checks and restores backups; every line is the answer to the command, not a record of a run.",
    "cli/benchmark.py": "A benchmark report to the terminal, a result table the person asked for, not a log.",
    "cli/benchmark_suite.py": "The benchmark suite's own report and step lines for the person at the terminal; it has no log file or severity.",
    "cli/config.py": "Shows, checks and changes settings; each line is the answer to the command.",
    "cli/forget_devices.py": "The count of trusted devices revoked: the answer to the command.",
    "cli/generate.py": "Prints the fingerprint lines the person asked for.",
    "cli/lift_lockout.py": "The count of lockouts lifted: the answer to the command.",
    "cli/lockouts.py": "Lists the lockouts: the answer to the command.",
    "cli/maintenance.py": "Shows whether maintenance mode is on: the answer to the command.",
    "cli/migrate.py": "The database state answer, and the one-time admin password box that must never reach a log file.",
    "cli/query.py": "Lists the sectors and systems the person asked for: the data, not a log.",
    "cli/queue.py": "Lists what is running or queued and whether work is on hold: the answer to the command.",
    "cli/reset.py": "The confirmation prompt before a wipe: an interactive display, not a record.",
    "cli/reset_two_factor.py": "Says whether two-factor sign-in is off for the account: the answer to the command.",
    "cli/smoke.py": "The smoke report table the person asked for.",
    "cli/species_check.py": "Lists the species found outside their range: the data the check exists to show.",
    "cli/update_guard.py": "Says whether an update is running and what is queued: the answer to the command.",
    "cli/version_history.py": "Says which version key was recorded: the answer to the command.",
    "galaxy/dials.py": "Prints the generated documentation block when run to rebuild the docs.",
    "physics/mathcheck.py": "Prints the math check report the person asked for.",
    "util/log.py": "The console writer itself, and the three notes it prints before the logger exists.",
    "util/secret_key.py": "Says the key was set: the answer to the command.",
}
"""dict: `file` to the reason it may still call `print()`. Everything a person reads as a record of a run (a progress line, a warning, an error, a
status line) goes through `log.normal` / `log.warning` / `log.error` so it carries the severity, time and stderr rule of the standard; what is
allowed here is only the answer a command was asked for (a listing, a report, a count, a prompt), which is data on stdout and not a log line.
A new file needs a reason a person can check."""


def find_prints(source, filename="<source>"):
    """`[(line, to_stderr)]`: every call of the builtin `print` in `source`, and whether it writes to a stderr stream."""
    found = []
    for node in ast.walk(ast.parse(source, filename)):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "print":
            to_stderr = any(keyword.arg == "file" and "stderr" in ast.unparse(keyword.value) for keyword in node.keywords)
            found.append((node.lineno, to_stderr))
    return sorted(found)


def _print_hits():
    hits = {}
    for relative, path in _modules():
        with open(path, encoding="utf-8") as handle:
            found = find_prints(handle.read(), path)
        if found:
            hits[relative] = found
    return hits


def test_no_module_prints_a_log_line():
    stray = {path: found for path, found in _print_hits().items() if path not in PRINT_ALLOWED}
    assert not stray, (
        "print() writes a line outside the logging layer. A progress, status, warning or error line goes through log.normal / log.warning / "
        "log.error (logging standard, rules 1 and 5); only the answer a command was asked for may print, and its file is named in PRINT_ALLOWED "
        "with the reason:\n  " + "\n  ".join(f"{path}:{found[0][0]}" for path, found in sorted(stray.items())))


def test_nothing_prints_to_stderr_but_the_logging_module():
    stray = {path: [line for line, to_stderr in found if to_stderr] for path, found in _print_hits().items() if path != THE_MODULE}
    stray = {path: lines for path, lines in stray.items() if lines}
    assert not stray, "print(..., file=sys.stderr) is a log line: use log.warning or log.error:\n  " + "\n  ".join(
        f"{path}:{lines[0]}" for path, lines in sorted(stray.items()))


def test_every_printing_file_is_still_printing_and_has_a_reason():
    hits = _print_hits()
    for path, reason in PRINT_ALLOWED.items():
        assert len(reason.split()) >= 8, f"{path} needs a real reason"
        assert path in hits, f"{path} no longer prints; remove it from PRINT_ALLOWED"


@pytest.mark.parametrize("source, found", [
    ("print('x')\n", [(1, False)]),
    ("import sys\nprint('x', file=sys.stderr)\n", [(2, True)]),
    ("console.print('x')\nlog.normal('x')\n", []),
])
def test_the_print_scan_finds_a_print_and_where_it_writes(source, found):
    assert find_prints(source) == found
