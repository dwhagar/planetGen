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
