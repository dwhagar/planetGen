# planetgen/util/settings_writer.py

"""
The settings writer (ADM.83)
============================

The write side of `planetgen.util.settings`, the reader. It saves the user's own `config.json`, the one file that
holds every program setting, for the configuration page (ADM.84), the backup, restore, import and export controls
(ADM.85) and `python3 -m planetgen.cli.config set`.

`config.json` is the user's copy of `config.json.example`: the example is only ever read by people, never by this
module, and the installers never overwrite the copy. If the copy does not exist the writer says so and creates
nothing; copying the example is a deployment step (INSTALL.md).

One save is:

1. **Check each change** against the settings model: the option exists, the caller's role may change it
   (`x-editable`, `x-min_role`), and the value fits (type, range, enum). Every problem is reported at its
   option's dotted path (`tile_cache.max_mb`), not as one block of text, and one problem refuses the whole save.
   The value written is the one the model parsed (`"30"` becomes `30`), never the raw input. Roles, lowest to
   highest: `"admin"` (the web), `"owner"` (the Owner of USR.2, from 8.2) and `"host"` (the command line:
   whoever can run it owns the server, so it may change the database, the secrets and the file-only options).
2. **Back up**: a settings-only archive through `planetgen.backup` (the old file, who is making the change,
   when; one `.tar.bz2` in the `backups` folder, `settings.backups_kept` of them kept). The archive holds the
   whole file, secrets included, which is why archives are 0600 and their download is for the Owner.
3. **Write atomically**: a temporary file in the same folder, `fsync`, the old file's mode and owner, `os.replace`,
   `fsync` of the folder. A reader sees the old file or the new one, never half of either. The document is edited
   in place as parsed JSON: every key the model does not know, and the order of the keys, are kept; only the
   changed options move.
4. **Audit**: one row per changed option (`config.update`, target `config:<option>`, detail
   `{"old":..,"new":..}`; a secret's values are written `(changed)`, never the text) through the caller's `audit`
   function, and one activity-log line (`DB config.update`) whatever the caller gave.

The whole read-modify-write holds a lock file next to the file, so two saves (two admins, or the page and the
command line) never interleave, and a form that carries the file's `file_hash` (`expect_hash`) is refused when
the file changed after the form was drawn.

Only keys that were set are stored: an option equal to the default is left out, and `RESET` removes a key so the
option falls back to the default. `config.json` gets a `config_version` on its first save.

Secrets (`x-secret`) are never returned by anything here that a page, a log or a download may show: `Change.shown`
and `redacted` mask them. `Plan.document_before` and `document_after` hold the raw documents and are for the
writer only.

This module imports no Flask.
"""

import copy
import hashlib
import json
import os
import threading
from dataclasses import dataclass, field

try:  # POSIX; elsewhere the in-process lock below is all there is
    import fcntl
except ImportError:  # pragma: no cover
    fcntl = None

from planetgen import backup
from planetgen.admin import activity_log
from planetgen.util import logpaths, settings

RESET = object()
"""A change value that removes the option from the file, so the default applies again."""

SECRET_MARK = "(changed)"
"""str: What an audit row, a diff and a log line show for a secret instead of its text."""

PART = "settings"
FILE = "config.json"
"""The backup part's name and the file inside it."""

LOCK_SUFFIX = ".lock"

MISSING = ("config.json does not exist. It is your own copy of config.json.example: copy the example to "
           "config.json and edit it (INSTALL.md, 'config.json'), then try again.")

FIX_HINT = ("config.json and its folder must belong to the web server's user: run ./update.sh on the server "
            "(or `sudo examples/apache/set-permissions.sh`) and try again.")

_ROLES = {"admin": 1, "owner": 2, "host": 3}
_process_lock = threading.Lock()


class SettingsWriteError(backup.BackupError):
    """config.json could not be read or written; the message says where and what to do."""


class StaleSettingsError(SettingsWriteError):
    """The file changed after the caller read it (`expect_hash`)."""


@dataclass
class Change:
    """One option that a save changes."""
    path: str
    old: object
    new: object
    secret: bool = False
    restart: bool = False
    reset: bool = False

    def shown(self, which):
        """`old` or `new` as a page, a log or an audit row may show it."""
        value = getattr(self, which)
        return SECRET_MARK if self.secret else value


@dataclass
class Plan:
    """What a save would do (the diff before saving), or why it can't."""
    changes: list = field(default_factory=list)
    errors: dict = field(default_factory=dict)
    shadowed: dict = field(default_factory=dict)
    document_before: dict = field(default_factory=dict)
    document_after: dict = field(default_factory=dict)
    file_hash: str = ""

    @property
    def ok(self):
        return not self.errors

    @property
    def restart(self):
        return [c.path for c in self.changes if c.restart]


@dataclass
class Result(Plan):
    """A `Plan` that was carried out (`saved`), with its backup."""
    saved: bool = False
    backup: str = None
    audit_errors: list = field(default_factory=list)


# ---------------------------------------------------------------------
# The file
# ---------------------------------------------------------------------

def config_path():
    return logpaths.CONFIG_PATH


def _folder(path):
    return os.path.dirname(os.path.abspath(path))


def write_problem():
    """Why this process can't write config.json (missing, or the file or its folder not writable), or `None`: a
    check made without writing, for the page."""
    path = config_path()
    if not os.path.isfile(path):
        return MISSING
    if not os.access(path, os.W_OK):
        return f"{path} is not writable by this user. {FIX_HINT}"
    if not os.access(_folder(path), os.W_OK | os.X_OK):
        return f"{_folder(path)} is not writable by this user, so the file can't be replaced safely. {FIX_HINT}"
    return None


def read_bytes(path):
    """The file's bytes, or `None` when it doesn't exist."""
    try:
        with open(path, "rb") as handle:
            return handle.read()
    except FileNotFoundError:
        return None
    except OSError as exc:
        raise SettingsWriteError(f"{path} can't be read: {exc.strerror or exc}") from exc


def file_hash():
    """A short hash of config.json's bytes (empty for none), for a form to detect a concurrent edit."""
    data = read_bytes(config_path())
    return hashlib.sha256(data).hexdigest()[:16] if data is not None else ""


def _parse(data, label):
    try:
        parsed = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as exc:
        raise SettingsWriteError(f"{label} is not valid JSON ({exc}); fix it by hand or restore a backup") from exc
    if not isinstance(parsed, dict):
        raise SettingsWriteError(f"{label} must contain a JSON object, not {type(parsed).__name__}")
    return parsed


def read_document():
    """config.json as a dict. Raises `SettingsWriteError` when it is missing or damaged."""
    path = config_path()
    data = read_bytes(path)
    if data is None:
        raise SettingsWriteError(MISSING)
    return _parse(data, path)


class _Lock:
    """An exclusive lock for the whole read-modify-write: a lock file between processes, a thread lock inside one."""

    def __init__(self, path):
        self.path = path + LOCK_SUFFIX
        self.handle = None

    def __enter__(self):
        _process_lock.acquire()
        try:
            if fcntl is not None:
                try:
                    fresh = not os.path.exists(self.path)
                    self.handle = open(self.path, "a")
                    if fresh:
                        backup.match_owner(self.path, _folder(self.path))
                    fcntl.flock(self.handle, fcntl.LOCK_EX)
                except OSError as exc:
                    self._close()
                    raise SettingsWriteError(f"The lock file {self.path} can't be used ({exc.strerror or exc}). "
                                             f"{FIX_HINT}") from exc
        except BaseException:
            _process_lock.release()
            raise
        return self

    def _close(self):
        if self.handle is not None:
            try:
                self.handle.close()
            finally:
                self.handle = None

    def __exit__(self, *exc_info):
        try:
            if self.handle is not None and fcntl is not None:
                fcntl.flock(self.handle, fcntl.LOCK_UN)
            self._close()
        finally:
            _process_lock.release()


def _atomic_write(path, data):
    try:
        backup.atomic_write(path, data)
    except backup.BackupError as exc:
        raise SettingsWriteError(str(exc)) from exc


def _dump(document):
    return (json.dumps(document, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


# ---------------------------------------------------------------------
# Reading the model
# ---------------------------------------------------------------------

def _fields():
    return {path: (info, meta) for path, info, meta in settings.iter_fields()}


def _allowed(meta, role):
    """Whether `role` may change the option `meta` describes; the second value says why not."""
    rank = _ROLES.get(role, 0)
    if not meta["x-editable"]:
        return rank >= _ROLES["host"], "can only be changed on the server, with python3 -m planetgen.cli.config"
    return rank >= _ROLES[meta["x-min_role"]], "can only be changed by the Owner"


def _get(tree, path):
    node = tree
    for part in path.split("."):
        if not isinstance(node, dict) or part not in node:
            return False, None
        node = node[part]
    return True, node


def _put(tree, path, value):
    parts = path.split(".")
    for part in parts[:-1]:
        if not isinstance(tree.get(part), dict):
            tree[part] = {}
        tree = tree[part]
    tree[parts[-1]] = value


def _drop(tree, path):
    parts = path.split(".")
    chain = [tree]
    for part in parts[:-1]:
        nxt = chain[-1].get(part)
        if not isinstance(nxt, dict):
            return
        chain.append(nxt)
    chain[-1].pop(parts[-1], None)
    for part, node in zip(reversed(parts[:-1]), reversed(chain[:-1])):  # a section left empty goes too
        if node.get(part) == {}:
            del node[part]


def _effective(settings_obj, path):
    value = settings._lookup(settings_obj, path)
    return list(value) if isinstance(value, list) else value


def redacted(document):
    """A copy of `document` with every secret option's value masked: `"(set)"` for a non-empty one, `""` for an
    empty one. For a page or an export; it never shows a secret."""
    shown = copy.deepcopy(document)
    for path, (_info, meta) in _fields().items():
        if meta["x-secret"]:
            found, value = _get(shown, path)
            if found:
                _put(shown, path, "(set)" if value else "")
    return shown


def _stamp_version(document):
    if "config_version" in document:
        return document
    return {"config_version": settings.CONFIG_VERSION, **document}


def _diff(before, after, options, known, environ, result):
    """Appends a `Change` to `result` for each of `options` whose effective value differs between the documents."""
    old_settings, _problems = settings.field_errors(before)
    new_settings, _problems = settings.field_errors(after)
    old_settings = old_settings or settings.Settings()
    for option in options:
        old, new = _effective(old_settings, option), _effective(new_settings, option)
        if old == new:
            continue
        _info, meta = known[option]
        result.changes.append(Change(option, old, new, secret=bool(meta["x-secret"]),
                                     restart=bool(meta["x-restart"]), reset=not _get(after, option)[0]))
        var = meta.get("x-env")
        if var and (environ.get(var) not in (None, "") or (meta.get("x-env-empty") and var in environ)):
            result.shadowed[option] = var


# ---------------------------------------------------------------------
# Planning a save
# ---------------------------------------------------------------------

def plan(changes, role="admin", environ=None):
    """
    Checks `changes` (`{dotted option path: value or RESET}`) against the model and the file, writing nothing: the
    diff a page shows before saving.

    Args:
        role (str): `"admin"` (the web), `"owner"` or `"host"` (the command line); see the module docstring.
        environ (Mapping): The environment to look for overriding variables in (default `os.environ`); they
            only produce `shadowed` notes, the value written never depends on them.

    Returns:
        Plan: `errors` maps an option path (or `""` for the file as a whole) to a message; `changes` lists the
        options that would really change (a value equal to what applies already is left out).
    """
    path = config_path()
    environ = os.environ if environ is None else environ
    result = Plan()
    data = read_bytes(path)
    if data is None:
        result.errors[""] = MISSING
        return result
    result.file_hash = hashlib.sha256(data).hexdigest()[:16]
    try:
        before = _parse(data, path)
    except SettingsWriteError as exc:
        result.errors[""] = str(exc)
        return result
    result.document_before = before
    known = _fields()
    after = copy.deepcopy(before)
    wanted = {}
    for option, value in changes.items():
        if option not in known:
            result.errors[option] = "is not a setting"
            continue
        _info, meta = known[option]
        permitted, why = _allowed(meta, role)
        if not permitted:
            result.errors[option] = why
        elif value is RESET:
            _drop(after, option)
            wanted[option] = RESET
        else:
            _put(after, option, value)
            wanted[option] = value
    checked, problems = settings.field_errors(after)
    if checked is None:
        for option, message in problems:
            if option in known:
                result.errors.setdefault(option, message)
            else:
                result.errors.setdefault("", f"{option}: {message}" if option else message)
    if result.errors:
        return result
    bare = copy.deepcopy(after)  # the same file without the options being set: what applies if they are left out
    for option in wanted:
        _drop(bare, option)
    baseline = settings.field_errors(bare)[0] or settings.Settings()
    for option, value in wanted.items():
        if value is RESET:
            continue
        normalized = _effective(checked, option)  # the model's reading of the text or number given
        if normalized == _effective(baseline, option):
            _drop(after, option)  # the default again: not worth storing
        else:
            _put(after, option, normalized)
    _diff(before, after, wanted, known, environ, result)
    result.document_after = _stamp_version(after) if result.changes else after
    return result


# ---------------------------------------------------------------------
# Saving
# ---------------------------------------------------------------------

def _audit(result, actor, audit, action):
    for change in result.changes:
        detail = json.dumps({"old": change.shown("old"), "new": change.shown("new")}, sort_keys=True)
        target = f"config:{change.path}"
        activity_log.event("DB", action, user=actor, target=target, detail=detail)
        if audit is not None:
            try:
                audit(action, target, detail)
            except Exception as exc:  # noqa: BLE001 -- the file is already written; say so, don't undo it
                result.audit_errors.append(f"{target}: {exc}")
                activity_log.event("DB", "config.audit_failed", user=actor, target=target, error=str(exc)[:200])


def _write(path, document, actor, action, result, archive):
    """Takes the settings-only backup archive (unless the caller already did) and writes `document` as the new
    file (the lock is held)."""
    try:
        if archive:
            result.backup = backup.take([PART], who=actor, reason=action).name
        _atomic_write(path, _dump(document))
        result.saved = True
    except backup.BackupError as exc:
        if isinstance(exc, SettingsWriteError):
            raise
        raise SettingsWriteError(str(exc)) from exc


def save(changes, actor, role="admin", audit=None, expect_hash=None, action="config.update", environ=None):
    """
    Plans and, when it is clean, writes `changes` (see `plan`). Nothing is written when any option has an error.

    Args:
        actor (str): Who is saving (a user name, `cli:<user>`), kept in the backup archive and the audit rows.
        audit (callable): `audit(action, target, detail)`, called once per changed option after the file is
            written; the web passes `planetgen.api.authz.audit`. A failure is listed in `Result.audit_errors`.
        expect_hash (str): The `file_hash` the caller's form was drawn from; a different one raises
            `StaleSettingsError`.

    Returns:
        Result: `saved` is false (and `errors` set) for a refused save, or when there was nothing to change.

    Raises:
        StaleSettingsError: The file changed since `expect_hash`.
        SettingsWriteError: The file or its folder can't be written.
    """
    path = config_path()
    with _Lock(path):
        planned = plan(changes, role=role, environ=environ)
        result = Result(**planned.__dict__)
        if expect_hash is not None and expect_hash != planned.file_hash:
            raise StaleSettingsError("The settings were changed by someone else since this form was opened; "
                                     "reload it and make the change again.")
        if result.errors or not result.changes:
            return result
        _write(path, result.document_after, actor, action, result, archive=True)
    _audit(result, actor, audit, action)
    return result


# ---------------------------------------------------------------------
# Restoring
# ---------------------------------------------------------------------

def _restore_plan(data, label, role, environ):
    """The `Plan` that restoring the backed-up file `data` would be, from its raw bytes (the lock is held).

    The command line (`host`) puts the whole file back. Anyone else gets only the options their role may change
    from the backup, and the rest of today's file stays as it is, so a web restore can never swap the database
    password or the secret key."""
    path = config_path()
    environ = os.environ if environ is None else environ
    result = Plan()
    archived = _parse(data, label)
    current_bytes = read_bytes(path)
    if current_bytes is None:
        result.errors[""] = MISSING
        return result, None
    result.file_hash = hashlib.sha256(current_bytes).hexdigest()[:16]
    current = _parse(current_bytes, path)
    result.document_before = current
    known = _fields()
    if _ROLES.get(role, 0) >= _ROLES["host"]:
        wanted = archived
    else:
        wanted = copy.deepcopy(current)
        for option, (_info, meta) in known.items():
            if _allowed(meta, role)[0]:
                found, value = _get(archived, option)
                if found:
                    _put(wanted, option, value)
                else:
                    _drop(wanted, option)
    checked, problems = settings.field_errors(wanted)
    if checked is None:
        for option, message in problems:
            result.errors.setdefault(option if option in known else "", message if option in known
                                     else f"{option}: {message}" if option else message)
        return result, None
    _diff(current, wanted, list(known), known, environ, result)
    result.document_after = wanted
    return result, wanted


def plan_restore(name, role="admin", environ=None):
    """What restoring the settings part of the archive `name` would change (the diff before a restore); nothing is
    written. Raises `SettingsWriteError` for an unreadable or altered archive."""
    try:
        data = backup.read_file(name, PART, FILE)
    except backup.BackupError as exc:
        raise SettingsWriteError(str(exc)) from exc
    with _Lock(config_path()):
        return _restore_plan(data, f"{name}/{PART}/{FILE}", role, environ)[0]


def restore_bytes(data, context):
    """
    The settings part's restore (`planetgen.backup.parts.SettingsPart`): `data` is the backed-up config.json. It
    must fit today's model; the options that differ get audit rows (`config.restore`). The caller
    (`backup.restore`) has already taken the safety archive of the file this replaces.

    Returns:
        Result: The options changed.

    Raises:
        SettingsWriteError: The backup does not fit today's model, or config.json is missing or can't be written.
    """
    path = config_path()
    with _Lock(path):
        planned, wanted = _restore_plan(data, "the backup's configuration", context.role, context.environ)
        result = Result(**planned.__dict__)
        if result.errors:
            raise SettingsWriteError("The backup's configuration can't be restored: "
                                     + "; ".join(f"{o + ': ' if o else ''}{m}" for o, m in result.errors.items()))
        if wanted != result.document_before:
            _write(path, _stamp_version(wanted), context.who, "config.restore", result, archive=False)
    _audit(result, context.who, context.audit, "config.restore")
    return result


def restore(name, actor, role="admin", audit=None, expect_hash=None, environ=None):
    """
    Restores the settings part of the archive `name` (the archive's other parts are left alone). The backup must
    fit today's model; the file being replaced is first backed up in a safety archive (`Result.backup`), and each
    option that differs gets an audit row (`config.restore`). Below the command line only the options the role
    may change are restored (see `_restore_plan`).

    Returns:
        Result: `errors` set and nothing written when the backup doesn't fit.

    Raises:
        SettingsWriteError: There is no such archive, or it is damaged or holds no configuration.
        StaleSettingsError: The file changed since `expect_hash`.
    """
    preview = plan_restore(name, role=role, environ=environ)
    if preview.errors:
        return Result(**preview.__dict__)
    if expect_hash is not None and expect_hash != preview.file_hash:
        raise StaleSettingsError("The settings were changed by someone else since this page was opened; "
                                 "reload it and try again.")
    try:
        report = backup.restore(name, [PART], backup.RestoreContext(who=actor, role=role, audit=audit,
                                                                    environ=environ))
    except backup.BackupError as exc:
        if isinstance(exc, SettingsWriteError):
            raise
        raise SettingsWriteError(str(exc)) from exc
    result = report.results[PART]
    result.backup = report.safety
    return result
