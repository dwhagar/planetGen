# planetgen/util/settings_view.py

"""
What the settings page shows (ADM.84)
=====================================

`describe` turns the settings model and `config.json` into plain data for the Admin settings page: every option
with its help text, type, limits, default, the value in the file, where the value in use comes from, whether the
caller's role may change it and whether a change needs a restart. `pending_restart` lists the options that need a
restart and were saved after this process started.

The page, the API (`/api/admin/settings`) and the command line all read the same model, so a new option written
with `settings.opt(...)` shows up here with no code of its own.

A secret option's value never leaves this module: it is `None`, with `is_set` telling whether one is stored.

This module imports no Flask.
"""

import os
import threading
from typing import Literal, get_args, get_origin

from planetgen.util import settings, settings_writer

GROUPS = {
    "site": "Site identity",
    "look": "Look",
    "search": "Search engines",
    "previews": "Link previews",
    "speed": "Speed and disk",
    "protection": "Abuse protection",
    "housekeeping": "Housekeeping",
    "database": "Database connection",
    "queue": "Queue server",
    "logging": "Logging",
    "proxy": "Reverse proxy",
    "security": "Security",
    "paths": "Folders and addresses",
    "wiki": "Wiki upload",
    "file": "Settings file",
}
"""dict: Each `x-group` and its heading, in the order the page lists them. The first seven are the groups an admin
works in; the rest hold what only the server's file can change (or the Owner), listed after them. Every option's
group must be named here (a test checks it)."""

KINDS = ("bool", "int", "float", "text", "list", "choice")

CONTROLS = {"bool": "switch", "int": "number", "float": "number", "text": "text", "list": "chips", "choice": "choice"}
"""dict: The control an option's type gives when its `x-control` names none."""

CONTROL_NAMES = frozenset(CONTROLS.values()) | {"textarea", "secret", "slider", "duration", "rate", "color", "asset", "language"}
"""frozenset: Every control the page may draw."""


def kind_of(annotation):
    """The control an option needs: `"bool"`, `"int"`, `"float"`, `"text"`, `"list"` (a list of text) or `"choice"`."""
    origin = get_origin(annotation)
    if origin is Literal:
        return "choice"
    if origin is list or annotation is list:
        return "list"
    if annotation is bool:
        return "bool"
    if annotation is int:
        return "int"
    if annotation is float:
        return "float"
    return "text"


def _limits(info):
    """`{"ge", "gt", "le", "lt", "min_length"}` for the bounds the field declares."""
    found = {}
    for item in info.metadata:
        for key in ("ge", "gt", "le", "lt", "min_length", "max_length"):
            value = getattr(item, key, None)
            if value is not None:
                found[key] = value
    return found


def _source(path, meta, document, environ):
    var = meta.get("x-env")
    if var and (environ.get(var) not in (None, "") or (meta.get("x-env-empty") and var in environ)):
        return "environment"
    return "config.json" if settings_writer._get(document, path)[0] else "default"


def describe(role="admin", environ=None):
    """
    Every option as the page lists it.

    Args:
        role (str): The caller's role (`settings_writer.allowed`).
        environ (Mapping): Where to look for overriding variables (default `os.environ`).

    Returns:
        dict: `file_hash` (what a save sends back as `expect_hash`), `problem` (why the page can't save, or
            `""`), `file_error` (why config.json can't be read, or `""`), `groups` (`[{"key", "title",
            "options"}]`, in page order) and `options`, the same options in one flat list. Each option is
            `{"path", "section", "name", "title", "summary", "description", "group", "control", "advanced", "kind", "choices",
            "limits", "unit", "multiline",
            "default", "value", "effective", "is_set", "source", "env", "restart", "secret", "editable",
            "min_role", "allowed", "reason"}`. `value` is what the file holds (the default where it holds
            nothing): what the form edits. `effective` is the value in use, which differs from `value` when an
            environment variable overrides the file. A secret's `value`, `default` and `effective` are `None`.
    """
    environ = os.environ if environ is None else environ
    file_error = ""
    document = {}
    try:
        document = settings_writer.read_document()
    except settings_writer.SettingsWriteError as exc:
        file_error = str(exc)
    problem = settings_writer.write_problem() or ""
    file_hash = settings_writer.file_hash() or ""
    from_file = settings.field_errors(document)[0]
    if from_file is None:
        file_error = file_error or "config.json has values that do not fit; the page shows the defaults for them."
        from_file = settings.field_errors({})[0]
    try:
        in_use = settings.get_settings(environ)
    except settings.SettingsError:
        in_use = from_file
    options = []
    for path, info, meta in settings.iter_fields():
        permitted, why = settings_writer.allowed(meta, role)
        annotation = meta["annotation"]
        kind = kind_of(annotation)
        secret = bool(meta["x-secret"])
        value = settings._lookup(from_file, path)
        effective = settings._lookup(in_use, path)
        default = settings.field_default(info)
        section, _dot, name = path.rpartition(".")
        options.append({
            "path": path,
            "section": section,
            "name": name,
            "title": meta["x-label"],
            "summary": meta["x-summary"],
            "description": meta["description"],
            "group": meta["x-group"],
            "control": meta["x-control"] or ("secret" if secret else "textarea" if meta["x-multiline"] else CONTROLS[kind]),
            "advanced": bool(meta["x-advanced"]),
            "kind": kind,
            "choices": list(get_args(annotation)) if kind == "choice" else [],
            "limits": _limits(info),
            "unit": meta.get("x-unit") or "",
            "multiline": bool(meta.get("x-multiline")),
            "default": None if secret else default,
            "value": None if secret else (list(value) if kind == "list" else value),
            "effective": None if secret else (list(effective) if kind == "list" else effective),
            "is_set": bool(value) if secret else None,
            "source": _source(path, meta, document, environ),
            "env": meta.get("x-env") or "",
            "restart": bool(meta["x-restart"]),
            "secret": secret,
            "editable": bool(meta["x-editable"]),
            "min_role": meta["x-min_role"],
            "allowed": permitted,
            "reason": "" if permitted else why,
        })
    order = list(GROUPS)
    present = sorted({option["group"] for option in options}, key=order.index)
    groups = [{"key": key, "title": GROUPS[key], "options": [o for o in options if o["group"] == key]}
              for key in present]
    return {"file_hash": file_hash, "problem": problem, "file_error": file_error, "groups": groups,
            "options": options}


# ---------------------------------------------------------------------
# Options that need a restart
# ---------------------------------------------------------------------

_lock = threading.Lock()
_startup = {"values": None}


def _restart_values():
    current = settings.get_settings()
    return {path: _plain(settings._lookup(current, path))
            for path, _info, meta in settings.iter_fields() if meta["x-restart"]}


def _plain(value):
    return list(value) if isinstance(value, list) else value


def remember_startup():
    """Records the restart-only options as they are now, for `pending_restart`; the web app calls it as it
    starts. A settings file that doesn't load leaves nothing recorded."""
    try:
        values = _restart_values()
    except settings.SettingsError:
        return
    with _lock:
        _startup["values"] = values


def forget_startup():
    """Forgets the record (tests)."""
    with _lock:
        _startup["values"] = None


def pending_restart():
    """
    The options that need a restart and now hold a different value than they did when this process started, as
    `[(path, title)]` in model order. Empty when nothing was recorded, or the file doesn't load (the process
    keeps the last good settings; `config check` says what is wrong with the file).
    """
    with _lock:
        started = _startup["values"]
    if started is None:
        return []
    try:
        now = _restart_values()
    except settings.SettingsError:
        return []
    titles = {path: meta["x-label"] for path, _info, meta in settings.iter_fields()}
    return [(path, titles[path]) for path in now if now[path] != started.get(path)]
