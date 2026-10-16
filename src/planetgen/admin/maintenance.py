# planetgen/admin/maintenance.py

"""
Maintenance mode (UX.112)
=========================

One switch that makes the whole site, pages and API alike, unavailable to
everyone but a signed-in admin (or an admin's API key): used while the
sector maps are rebuilt. A visitor gets a "the site is being updated" page
(HTTP 503, `Retry-After`) or, from the API, a JSON 503 saying the same,
with the admin's expected finish time when one was given.

The state is one small JSON file, `maintenance.json`, in the state folder
(`/var/lib/planetGen/state`, `PLANETGEN_STATE_DIR` overrides it, and
`PLANETGEN_MAINTENANCE_FILE` names the file itself). The web user writes it
(the admin page, and `warm_map` when the map is rebuilt), so the installer
gives that folder to Apache's user (`examples/apache/create-cache-dir.sh`, run
by `install.sh` and `update.sh`): the install folder above it is root's, and a
file written there from the web app fails with "Permission denied". It needs
no database, so the switch works while the galaxy database is being
migrated, and the web app, the command line (`python3 -m
planetgen.cli.maintenance`) and the admin page all read and write the same
file. A missing file means maintenance is off. A file that exists but can't
be read counts as on, because the safe mistake is a closed site.

```json
{"enabled": true, "eta": "about two hours", "reason": "", "cause": "", "since": "2026-10-15T08:00:00Z", "by": "boss"}
```

The site is open in its normal state. Two things close it for a time:
an admin (the admin page or the command), and the sector-grid switch
(`SECTOR_SWITCH`): `planetgen.cli.migrate` turns it on before it converts a
database from the cylinder grid to the spherical one, with that `cause`, and
the map rebuild turns it off again with `turn_off_for(SECTOR_SWITCH)` once the
new Galaxy Map is built (an admin's own maintenance, with another or no
cause, is left alone).

This module imports no Flask; `web/maintenance_gate.py` is the request hook.
"""

import json
import os
import tempfile
import threading
from datetime import datetime, timezone

FILE_ENV = "PLANETGEN_MAINTENANCE_FILE"
"""str: The environment variable that names the state file."""

STATE_DIR_ENV = "PLANETGEN_STATE_DIR"
"""str: The environment variable that names the state folder."""

DEFAULT_STATE_DIR = "/var/lib/planetGen/state"
"""str: The folder the installer gives to Apache's user (`examples/apache/deploy-paths.py` says the same)."""

FIX_HINT = ("The folder must belong to the web server's user: run ./update.sh on the server "
            "(or `sudo examples/apache/create-cache-dir.sh`) and try again.")
"""str: What to do about a switch that can't be written."""


class MaintenanceWriteError(OSError):
    """The maintenance switch could not be written; the message says where and what to do."""

SECTOR_SWITCH = "sector-switch"
"""str: The `cause` of the maintenance the move to the spherical sector grid starts."""

SECTOR_SWITCH_REASON = "The Galaxy Map is being rebuilt for the new sector grid."
"""str: What visitors are told during the sector-grid switch (the page's lead line)."""

MAX_ETA_CHARS = 200
"""int: The longest expected-finish text accepted."""

_lock = threading.Lock()
_cache = {"key": None, "state": None}


def state_dir():
    """The state folder: `PLANETGEN_STATE_DIR`, else `DEFAULT_STATE_DIR`."""
    return os.environ.get(STATE_DIR_ENV) or DEFAULT_STATE_DIR


def state_path():
    """The state file: `PLANETGEN_MAINTENANCE_FILE`, else `maintenance.json` in the state folder."""
    return os.environ.get(FILE_ENV) or os.path.join(state_dir(), "maintenance.json")


def write_problem():
    """
    Why the switch can't be written by this process, or `None` when it can: the file's folder is missing (and
    its parent is not writable) or not writable. A check made without writing, for the admin page.
    """
    folder = os.path.dirname(state_path()) or "."
    target = folder
    while not os.path.isdir(target):
        parent = os.path.dirname(target)
        if not parent or parent == target:
            break
        target = parent
    if os.access(target, os.W_OK | os.X_OK, effective_ids=os.access in os.supports_effective_ids):
        return None
    return f"The maintenance switch can't be written: {target} is not writable by the web server's user. {FIX_HINT}"


def _off():
    return {"enabled": False, "eta": "", "reason": "", "cause": "", "since": None, "by": None}


def read():
    """
    The current state: `{"enabled", "eta", "since", "by"}`. Re-read only
    when the file changes (a stat per call). Unreadable means on.
    """
    path = state_path()
    try:
        stat = os.stat(path)
    except FileNotFoundError:
        return _off()
    except OSError:
        return dict(_off(), enabled=True)
    key = (path, stat.st_mtime_ns, stat.st_size)
    with _lock:
        if _cache["key"] == key:
            return dict(_cache["state"])
        try:
            with open(path, "r", encoding="utf-8") as handle:
                data = json.load(handle)
            if not isinstance(data, dict):
                raise ValueError("not a JSON object")
            state = {"enabled": data.get("enabled") is True,
                     "eta": str(data.get("eta") or "")[:MAX_ETA_CHARS],
                     "reason": str(data.get("reason") or "")[:MAX_ETA_CHARS],
                     "cause": str(data.get("cause") or "")[:40],
                     "since": data.get("since") if isinstance(data.get("since"), str) else None,
                     "by": data.get("by") if isinstance(data.get("by"), str) else None}
        except (OSError, ValueError):
            state = dict(_off(), enabled=True)
        _cache["key"], _cache["state"] = key, state
        return dict(state)


def is_on():
    """Whether maintenance mode is on."""
    return read()["enabled"]


def _write(state):
    path = state_path()
    folder = os.path.dirname(path) or "."
    try:
        os.makedirs(folder, exist_ok=True)
        handle, temp = tempfile.mkstemp(prefix=".maintenance-", dir=folder)
    except OSError as exc:
        raise MaintenanceWriteError(
            f"The maintenance switch {path} can't be written ({exc.strerror or exc}). {FIX_HINT}") from exc
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as out:
            json.dump(state, out)
            out.write("\n")
        os.chmod(temp, 0o644)
        os.replace(temp, path)
    except BaseException as exc:
        if os.path.exists(temp):
            os.unlink(temp)
        if isinstance(exc, OSError):
            raise MaintenanceWriteError(
                f"The maintenance switch {path} can't be written ({exc.strerror or exc}). {FIX_HINT}") from exc
        raise
    with _lock:
        _cache["key"] = None


def turn_on(eta="", by=None, reason=None, cause=None):
    """
    Turns maintenance mode on.

    Args:
        eta (str): A free-text expected finish ("about two hours", "by
            Thursday evening"); empty when there is no estimate.
        by (str, optional): Who switched it on.
        reason (str, optional): The lead line visitors see instead of the
            standard one. `None` keeps the one already set while it is on.
        cause (str, optional): What switched it on (`SECTOR_SWITCH`), so
            `turn_off_for` clears only its own. `None` keeps the current one
            while it is on.

    Raises:
        ValueError: `eta` is longer than `MAX_ETA_CHARS`.
    """
    eta = " ".join((eta or "").split())
    if len(eta) > MAX_ETA_CHARS:
        raise ValueError(f"The expected finish must be at most {MAX_ETA_CHARS} characters.")
    current = read()
    if reason is None:
        reason = current["reason"] if current["enabled"] else ""
    if cause is None:
        cause = current["cause"] if current["enabled"] else ""
    _write({"enabled": True, "eta": eta, "reason": reason, "cause": cause,
            "since": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "by": by})


def turn_off(by=None):
    """Turns maintenance mode off."""
    _write(dict(_off(), by=by))


def turn_off_for(cause, by=None):
    """
    Turns maintenance mode off only if `cause` is what turned it on (an
    admin's own maintenance stays).

    Returns:
        bool: Whether it was turned off.
    """
    state = read()
    if not state["enabled"] or state["cause"] != cause:
        return False
    turn_off(by=by)
    return True
