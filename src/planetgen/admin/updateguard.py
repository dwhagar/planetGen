# planetgen/admin/updateguard.py

"""
The update guard (OPS.41)
=========================

While `update.sh` (or a manual `planetgen.cli.migrate` that has something to migrate) changes the checkout and the
databases, the site is closed and the update holds a lock:

* **The site is closed** (OPS.41): maintenance mode (`admin/maintenance.py`) is turned on with the cause
  `UPDATE_CAUSE` before anything changes and turned off only when the update has finished. An update that fails
  leaves it on, with a message naming the step and the command that turns it off, because the database may be
  half migrated. An admin's own maintenance stays theirs; the sector-grid switch takes the closing over
  (`cli/migrate.py`) and the Galaxy Map warm-up opens the site.
* **The update holds a lock**, `update.lock` in the state folder, for as long as its process lives. `is_updating`
  is true while the lock is held by a live process or the maintenance state still has the cause `UPDATE_CAUSE`.

This module imports no Flask and nothing outside the standard library and `admin/maintenance.py`.
"""

import contextlib
import json
import os
from datetime import datetime, timezone

from planetgen.admin import maintenance

UPDATE_CAUSE = "update"
"""str: The maintenance `cause` of an update."""

UPDATE_REASON = "The site is being updated and will be back in a moment."
"""str: What visitors are told while an update runs."""

LOCK_ENV = "PLANETGEN_UPDATE_LOCK_FILE"
"""str: The environment variable that names the lock file."""

OWNER_ENV = "PLANETGEN_UPDATE_PID"
"""str: The process id of the update that runs a command (set by `update.sh`), so a `migrate` it starts knows the
update already holds the lock."""

LOCK_NAME = "update.lock"

TURN_OFF_HINT = "python -m planetgen.cli.maintenance off"
"""str: The command that opens the site again after a failed update."""


def _now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def pid_alive(pid):
    """Whether a process with this id is running here (one of another user counts)."""
    if not isinstance(pid, int) or pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    try:
        with open(f"/proc/{pid}/stat", "rb") as handle:
            return handle.read().rpartition(b")")[2].split()[0] != b"Z"      # a zombie has ended
    except (OSError, IndexError):
        return True


def lock_path():
    """The lock file: `PLANETGEN_UPDATE_LOCK_FILE`, else `update.lock` in the state folder."""
    return os.environ.get(LOCK_ENV) or os.path.join(maintenance.state_dir(), LOCK_NAME)


def read_lock():
    """The update's lock (`{"pid", "since", "by", "kept"}`) while its process is alive, else `None`. The lock of a
    process that died is ignored (and replaced by the next update)."""
    try:
        with open(lock_path(), "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict) or not pid_alive(data.get("pid")):
        return None
    return data


def is_updating():
    """Whether an update runs: its lock is held by a live process, or maintenance mode still has the update's cause
    (an update that failed leaves the site closed, and nothing may run until an admin opens it)."""
    state = maintenance.read()
    return read_lock() is not None or (state["enabled"] and state["cause"] == UPDATE_CAUSE)


def inside_update():
    """Whether this process was started by the update that holds the lock (`OWNER_ENV`)."""
    lock = read_lock()
    try:
        return lock is not None and lock["pid"] == int(os.environ.get(OWNER_ENV, ""))
    except ValueError:
        return False


def _write_lock(pid, by, kept):
    path = lock_path()
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    temp = f"{path}.{os.getpid()}.tmp"
    with open(temp, "w", encoding="utf-8") as handle:
        json.dump({"pid": pid, "since": _now(), "by": by, "kept": kept}, handle)
        handle.write("\n")
    os.chmod(temp, 0o644)
    os.replace(temp, path)


def _remove_lock():
    with contextlib.suppress(OSError):
        os.unlink(lock_path())


# ---------------------------------------------------------------------
# Beginning and ending an update
# ---------------------------------------------------------------------

def begin(pid, by="update.sh", write=print):
    """
    Starts an update: takes the lock and closes the site, before the first change.

    Args:
        pid (int): The process that holds the lock for as long as it lives (`update.sh`'s shell).

    Returns:
        bool: Whether the update may go on. `False` changes nothing: another update runs, or the site could not be
        closed.
    """
    existing = read_lock()
    if existing is not None and existing["pid"] != pid:
        write(f"error: an update is already running (process {existing['pid']}, since {existing['since']}).")
        return False
    state = maintenance.read()
    already_closed = state["enabled"] and state["cause"] != UPDATE_CAUSE    # an admin's own closing stays theirs
    try:
        _write_lock(pid, by, kept=already_closed)
        if not already_closed:
            maintenance.turn_on(by=by, reason=UPDATE_REASON, cause=UPDATE_CAUSE)
    except OSError as exc:
        _remove_lock()
        write(f"error: could not close the site for the update ({exc}); nothing was changed.")
        return False
    return True


def end(pid, ok, step=None, write=print):
    """
    Ends the update `pid` holds the lock for. A good run opens the site again (unless an admin had closed it
    before, or the sector-grid switch took the closing over: the Galaxy Map warm-up opens it then). A failed one
    leaves it closed, says which step failed and how to open it, and still releases the lock so the operator can run
    the update again.

    Returns:
        bool: Whether the site is open afterwards.
    """
    lock = read_lock()
    if lock is not None and lock["pid"] != pid:
        write(f"error: the update lock belongs to process {lock['pid']}, not {pid}; left alone.")
        return not maintenance.is_on()
    _remove_lock()
    state = maintenance.read()
    if ok:
        if lock is not None and not lock.get("kept"):
            maintenance.turn_off_for(UPDATE_CAUSE, by="update")
        return not maintenance.is_on()
    if state["enabled"] and state["cause"] == UPDATE_CAUSE:
        write(f"error: the update failed{f' at {step}' if step else ''}. The site stays closed (maintenance mode) "
              f"because the database may be half migrated. Fix the problem and run ./update.sh again; to open the "
              f"site without updating: {TURN_OFF_HINT}")
    else:
        write(f"error: the update failed{f' at {step}' if step else ''}.")
    return False
