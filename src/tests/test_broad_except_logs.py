"""
A broad `except` says what it went past (OPS.168, TEST.183).

`except Exception` (and a bare `except`, and `except BaseException`) is how this code chooses to degrade: a cache that cannot be written, a
note that cannot be kept, a probe that fails. Each of those used to leave no trace, so a debug run was blind exactly where a fault hid. This test
reads every module of the package as a syntax tree and fails on a broad handler that neither re-raises nor logs. The way to satisfy it is
`log.swallowed(area, what, exc)` (one debug line, nothing at normal level), or a real log line when the failure matters.

`ALLOWED` names the handlers that cannot log, by file and function, with the reason. Migrations (`db/migrations`) are frozen history and are not scanned.
"""
import ast
import os

PACKAGE = os.path.join(os.path.dirname(__file__), "..", "planetgen")

SKIPPED = ("db/migrations/",)

LOGGING_MACHINERY = "The logging layer itself: a handler that logged here would call itself, or lose the line it is writing."

ALLOWED = {
    ("util/log.py", "_attach_error_log"): LOGGING_MACHINERY,
    ("util/log.py", "_configured_areas"): LOGGING_MACHINERY,
    ("util/log.py", "config_debug"): LOGGING_MACHINERY,
    ("util/log.py", "_package_version"): LOGGING_MACHINERY,
    ("util/log.py", "_configure_debug_log"): LOGGING_MACHINERY,
    ("util/log.py", "_CallbackHandler.emit"): LOGGING_MACHINERY,
    ("util/log.py", "_ErrorLogHandler.emit"): LOGGING_MACHINERY,
    ("util/log.py", "_RelayHandler.flush"): LOGGING_MACHINERY,
    ("util/log.py", "_referenced_values"): LOGGING_MACHINERY,
    ("util/log.py", "RelayReader.drain"): LOGGING_MACHINERY,
    ("util/log.py", "_masked_record"): LOGGING_MACHINERY,
    ("util/secrets_mask.py", "_settings_secrets"): LOGGING_MACHINERY,
}

LOG_NAMES = ("log.", "planetgen_log.", "logging.")
LOG_ENDINGS = ("record_error", "swallowed", "degraded", "explain", "explain_lost")


def _is_broad(handler):
    if handler.type is None:
        return True
    kinds = handler.type.elts if isinstance(handler.type, ast.Tuple) else [handler.type]
    return any(getattr(kind, "id", getattr(kind, "attr", "")) in ("Exception", "BaseException") for kind in kinds)


def _is_quiet(handler):
    for node in ast.walk(handler):
        if isinstance(node, ast.Raise):
            return False
        if isinstance(node, ast.Call):
            name = ast.unparse(node.func)
            if name.startswith(LOG_NAMES) or "logger" in name.lower() or name.endswith(LOG_ENDINGS):
                return False
    return True


def quiet_handlers():
    """`[(file, function)]` of every broad handler that neither re-raises nor logs."""
    found = []
    for folder, _dirs, names in os.walk(PACKAGE):
        for name in sorted(names):
            if not name.endswith(".py"):
                continue
            path = os.path.join(folder, name)
            relative = os.path.relpath(path, PACKAGE).replace(os.sep, "/")
            if relative.startswith(SKIPPED):
                continue
            with open(path, encoding="utf-8") as source:
                tree = ast.parse(source.read())
            parents = {child: parent for parent in ast.walk(tree) for child in ast.iter_child_nodes(parent)}
            for node in ast.walk(tree):
                if isinstance(node, ast.ExceptHandler) and _is_broad(node) and _is_quiet(node):
                    names_up, scope = [], node
                    while scope in parents:
                        scope = parents[scope]
                        if isinstance(scope, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                            names_up.append(scope.name)
                    found.append((relative, ".".join(reversed(names_up)) or "<module>"))
    return found


def test_every_broad_except_reraises_or_logs():
    offenders = sorted({site for site in quiet_handlers() if site not in ALLOWED})
    assert not offenders, (f"a broad `except` that neither re-raises nor logs: {offenders}. Add `log.swallowed(area, what, exc)` "
                           f"(`except Exception as exc:`), or name it in ALLOWED with the reason it cannot log.")


def test_allowed_names_only_handlers_that_exist():
    stale = sorted(site for site in ALLOWED if site not in set(quiet_handlers()))
    assert not stale, f"in ALLOWED but no longer a quiet broad handler: {stale}"
