# planetgen/admin/webupdate.py

"""
The Update button of the admin pages (SEC.43): how a request to run `update.sh` gets from the web user to root with no
shell and no argument.

The hand-off is two empty files, nothing else:

* The web process (`www-data`) creates `update-request` in the state folder, and `update-request-cancel` as well when
  the admin agreed to cancel running work (what `update.sh --yes` means). Their content is never read. A root-owned
  systemd path unit (`planetgen-web-update.path`, installed by `examples/maintenance/install-web-update.sh`) starts
  `planetgen-web-update.service` when `update-request` exists, and that service runs a root-owned copy of
  `examples/maintenance/planetgen-web-update.sh`. The script removes both request files and runs `update.sh`, with
  `--yes` only when the cancel file was there. Nothing the browser sends reaches the command line.
* The script writes what happened to a folder only root can write (`status_dir`, default `/var/lib/planetgen-update`,
  readable by all): `status.json` (state, start and finish times, exit status) and `update.log`. The web process only
  reads them, so a web user who plants a symlink cannot make root write somewhere else.

This module imports no Flask. `view` is what the page and the API show; `request` is what the button does.
"""

import json
import os
import time

from planetgen.admin import maintenance, updateguard

STATUS_DIR_ENV = "PLANETGEN_UPDATE_STATUS_DIR"
"""str: The environment variable that names the folder root writes the update's status and log to."""

DEFAULT_STATUS_DIR = "/var/lib/planetgen-update"
"""str: Where `install-web-update.sh` makes it."""

REQUEST_NAME = "update-request"
"""str: The file in the state folder whose existence starts an update."""

CANCEL_NAME = "update-request-cancel"
"""str: The file whose existence, beside the request, adds `--yes` (cancel the work that is running)."""

STATUS_NAME = "status.json"
LOG_NAME = "update.log"

LOG_TAIL_BYTES = 24 * 1024
"""int: How much of the end of the log the page shows."""

STALE_REQUEST_SECONDS = 120
"""int: A request file this old that nothing picked up means the path unit is not running."""

STATES = ("not_installed", "idle", "requested", "running", "succeeded", "failed")
"""tuple: What `view` calls the state of the update machinery."""


class UpdateRequestError(Exception):
    """A request that cannot be made; the message says why, for the admin."""


def status_dir():
    """The folder root writes `status.json` and `update.log` to: `PLANETGEN_UPDATE_STATUS_DIR`, else the default."""
    return os.environ.get(STATUS_DIR_ENV) or DEFAULT_STATUS_DIR


def request_path():
    return os.path.join(maintenance.state_dir(), REQUEST_NAME)


def cancel_path():
    return os.path.join(maintenance.state_dir(), CANCEL_NAME)


def installed():
    """Whether `install-web-update.sh` has been run here: its status file is there."""
    return os.path.isfile(os.path.join(status_dir(), STATUS_NAME))


def _read_status():
    """`status.json` as a dict (empty when missing or unreadable)."""
    try:
        with open(os.path.join(status_dir(), STATUS_NAME), "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _text(value, limit=64):
    return value[:limit] if isinstance(value, str) else None


def log_tail(limit=LOG_TAIL_BYTES):
    """The last `limit` bytes of the update's log as text (empty when there is none). A cut first line is dropped."""
    path = os.path.join(status_dir(), LOG_NAME)
    try:
        with open(path, "rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - limit))
            data = handle.read(limit)
    except OSError:
        return ""
    text = data.decode("utf-8", errors="replace")
    if size > limit and "\n" in text:
        text = text.split("\n", 1)[1]
    return text


def _request_age():
    """Seconds since the request file was made, or `None` when there is none."""
    try:
        return max(0.0, time.time() - os.lstat(request_path()).st_mtime)
    except OSError:
        return None


def view():
    """
    What the Update page and `GET /api/admin/update` show:
    `{"state", "installed", "can_request", "reason", "started_at", "finished_at", "exit_code", "cancel_running",
    "log", "stale_request"}`. `state` is one of `STATES`; `can_request` is whether a request would be
    accepted now and `reason` says why not when it would not.
    """
    status = _read_status()
    on = installed()
    age = _request_age()
    lock = updateguard.read_lock()
    recorded = status.get("state") if status.get("state") in ("idle", "running", "succeeded", "failed") else "idle"
    if not on:
        state = "not_installed"
    elif lock is not None or recorded == "running":
        state = "running"
    elif age is not None:
        state = "requested"
    else:
        state = recorded
    reason = ""
    if state == "not_installed":
        reason = ("The Update button is not set up on this server. Run examples/maintenance/install-web-update.sh "
                  "once as root.")
    elif state == "running":
        reason = "An update is running."
    elif state == "requested":
        reason = "An update has been requested and is about to start."
    else:
        problem = maintenance.write_problem()
        if problem:
            reason = problem
    exit_code = status.get("exit_code")
    return {
        "state": state,
        "installed": on,
        "can_request": not reason,
        "reason": reason,
        "started_at": _text(status.get("started_at")),
        "finished_at": _text(status.get("finished_at")),
        "exit_code": exit_code if isinstance(exit_code, int) and not isinstance(exit_code, bool) else None,
        "cancel_running": status.get("cancel_running") is True,
        "log": log_tail() if on else "",
        "stale_request": state == "requested" and age is not None and age > STALE_REQUEST_SECONDS,
    }


def _make(path):
    """Creates the empty file `path`, failing when it already exists (a symlink counts)."""
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0), 0o660)
    os.close(fd)


def request(cancel_running=False):
    """
    Asks root to run `update.sh`: makes the request file (and the cancel file when `cancel_running`). The caller has
    already checked who may ask (the Owner, with a recent password).

    Raises:
        UpdateRequestError: Not set up here, an update is running or already requested, or the state folder cannot be
            written.
    """
    current = view()
    if not current["can_request"]:
        raise UpdateRequestError(current["reason"])
    try:
        # The cancel file first: root looks for it only once the request file is there.
        if cancel_running:
            _make(cancel_path())
        _make(request_path())
    except FileExistsError as exc:
        _forget(cancel_path())
        raise UpdateRequestError("An update has been requested and is about to start.") from exc
    except OSError as exc:
        _forget(cancel_path())
        raise UpdateRequestError(f"The request could not be written ({exc.strerror or exc}); is the state "
                                 f"folder {maintenance.state_dir()} writable by the web server's user?") from exc


def _forget(path):
    try:
        os.unlink(path)
    except OSError:
        pass


def withdraw():
    """Takes back a request root has not started (both files); returns whether there was one."""
    found = _request_age() is not None
    _forget(request_path())
    _forget(cancel_path())
    return found

