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

CATEGORY_TITLES = {
    "site": "Site",
    "security": "Security",
    "database": "Database",
    "cache": "Caches",
    "jobs": "Background jobs",
    "backups": "Backups",
    "wiki": "Wiki upload",
    "logging": "Logging",
    "proxy": "Reverse proxy",
    "file": "This file",
    "general": "General",
}
"""dict: The heading of each `x-category`, in the order the page lists them (a category not named here follows,
under its own name)."""

KINDS = ("bool", "int", "float", "text", "list", "choice")


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
        for key in ("ge", "gt", "le", "lt", "min_length"):
            value = getattr(item, key, None)
            if value is not None:
                found[key] = value
    return found


def _title(path):
    return path.replace("_", " ").replace(".", " > ")


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
            `""`), `file_error` (why config.json can't be read, or `""`), `groups` (`[{"category", "title",
            "options"}]`, in page order) and `options`, the same options in one flat list. Each option is
            `{"path", "section", "name", "title", "description", "category", "kind", "choices", "limits", "unit",
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
            "title": _title(path),
            "description": meta["description"],
            "category": meta["x-category"],
            "kind": kind,
            "choices": list(get_args(annotation)) if kind == "choice" else [],
            "limits": _limits(info),
            "unit": meta.get("x-unit") or "",
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
    order = list(CATEGORY_TITLES)
    present = []
    for option in options:
        if option["category"] not in present:
            present.append(option["category"])
    present.sort(key=lambda c: order.index(c) if c in order else len(order))
    groups = [{"category": c, "title": CATEGORY_TITLES.get(c, c.capitalize()),
               "options": [o for o in options if o["category"] == c]} for c in present]
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
    return [(path, _title(path)) for path in now if now[path] != started.get(path)]
