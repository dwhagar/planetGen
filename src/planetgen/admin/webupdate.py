# planetgen/admin/webupdate.py

"""
The Update button of the admin pages (SEC.43): how a request to run `update.sh` gets from the web user to root with no
shell and no argument.

The hand-off is one empty file. The web process (`www-data`) creates `update.request` in the state folder and never
writes anything else; its content is never read. The root-owned `planetgen-update.path` (`examples/maintenance`,
installed by `install.sh` and `update.sh`) sees it and starts `planetgen-update.service`, which removes the request and
runs `update.sh` as root with no argument. A request made while an update runs stays in place and runs after it. Whether
running work is canceled is the update guard's rule and the service's `--yes` choice, not the browser's.

Everything this module reads is public: the units' presence (`/etc/systemd/system`), the update's lock and the
maintenance switch, and `update.log` in the log folder, which `update.sh` appends every line to. This module imports no
Flask. `view` is what the page and the API show; `request` is what the button does.
"""

import os
import subprocess
import time

from planetgen.admin import maintenance, updateguard
from planetgen.util import logpaths, logrecord

REQUEST_NAME = "update.request"
"""str: The file in the state folder whose existence asks for an update."""

UNIT_DIR_ENV = "PLANETGEN_SYSTEMD_UNIT_DIR"
"""str: Where the units are looked for (default `/etc/systemd/system`; tests point it elsewhere)."""

PATH_UNIT = "planetgen-update.path"
SERVICE_UNIT = "planetgen-update.service"

LOG_TAIL_BYTES = 24 * 1024
"""int: How much of the end of the log the page shows."""

STALE_REQUEST_SECONDS = 120
"""int: A request this old that nothing picked up while no update runs means the path unit is not running."""

SYSTEMCTL_TIMEOUT_SECONDS = 3
"""int: How long `systemctl show` may take; a slow answer is no answer."""

STATES = ("not_installed", "idle", "requested", "running")
"""tuple: What `view` calls the state."""


class UpdateRequestError(Exception):
    """A request that cannot be made; the message says why, for the admin."""


def request_path():
    return os.path.join(maintenance.state_dir(), REQUEST_NAME)


def unit_dir():
    return os.environ.get(UNIT_DIR_ENV) or "/etc/systemd/system"


def installed():
    """Whether both root-owned units are in place (`install_maintenance_units` in `scripts/deploy-common.sh`)."""
    return all(os.path.isfile(os.path.join(unit_dir(), name)) for name in (PATH_UNIT, SERVICE_UNIT))


def log_path():
    try:
        return logpaths.job_log_path(logpaths.UPDATE_LOG_NAME)
    except (OSError, ValueError):
        return None


def log_tail(limit=LOG_TAIL_BYTES):
    """`(text, modified)`: the last `limit` bytes of `update.log` (a cut first line dropped) and its modification time
    (Unix seconds); `("", None)` when there is no log."""
    path = log_path()
    if not path:
        return "", None
    try:
        with open(path, "rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - limit))
            data = handle.read(limit)
        modified = os.stat(path).st_mtime
    except OSError:
        return "", None
    text = data.decode("utf-8", errors="replace")
    if size > limit and "\n" in text:
        text = text.split("\n", 1)[1]
    return text, modified


def _request_age():
    """Seconds since the request file was made, or `None` when there is none."""
    try:
        return max(0.0, time.time() - os.lstat(request_path()).st_mtime)
    except OSError:
        return None


def service_status():
    """
    What systemd says about the update service, from `systemctl show` (it needs no privilege): `{"active", "result",
    "status"}` (`ActiveState`, `Result` and `ExecMainStatus`), or `None` when systemctl is missing or does not answer.
    A unit that failed stays `failed` (`result` `exit-code`, or `start-limit-hit` when ten starts came in ten minutes)
    until the next start, so this is how a start that never reached `update.log` shows.
    """
    try:
        done = subprocess.run(
            ["systemctl", "show", SERVICE_UNIT, "-p", "ActiveState", "-p", "Result", "-p", "ExecMainStatus"],
            capture_output=True, text=True, timeout=SYSTEMCTL_TIMEOUT_SECONDS, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    if done.returncode != 0:
        return None
    values = dict(line.split("=", 1) for line in done.stdout.splitlines() if "=" in line)
    if "ActiveState" not in values:
        return None
    return {"active": values["ActiveState"], "result": values.get("Result", ""), "status": values.get("ExecMainStatus", "")}


def last_run(text):
    """
    `(finished, last_error)` for the last run in a log's text: whether its closing line ("=== Update done after ...")
    came, and the last error or critical line after the previous run's closing line (`None` when there is none).
    """
    finished = False
    last_error = None
    for line in text.splitlines():
        severity = logrecord.display_severity(line)
        if severity == "info" and "done after" in line and line.rstrip().endswith(logrecord.DIVIDER):
            finished = True
            last_error = None
        elif severity in ("error", "critical"):
            finished = False
            last_error = line
    return finished, last_error


def _outcome(text, status):
    """`(outcome, detail)`: `"failed"` with what is known when the last update did not finish, `"done"` when it did,
    `(None, "")` when nothing says."""
    finished, last_error = last_run(text)
    if status and status["active"] == "failed":
        if status["result"] == "start-limit-hit":
            detail = ("systemd stopped starting it after ten tries in ten minutes. Clear that with "
                      "`systemctl reset-failed planetgen-update.service`, then press Update again.")
        elif status["status"] == "1":
            detail = ("It stopped with exit status 1 before it changed anything or part way through; the usual reason is "
                      "that generation was running, which an update started this way never cancels.")
        else:
            detail = f"The service ended with {status['result'] or 'an error'} (exit status {status['status'] or '?'})."
        return "failed", detail + (f" The last error in the log: {last_error}" if last_error else "")
    if last_error:
        return "failed", f"The last error in the log: {last_error}"
    if finished:
        return "done", ""
    return None, ""


def view():
    """
    What the Update page and `GET /api/admin/update` show: `state` (one of `STATES`), `installed`, `can_request` (whether
    a request would be accepted now) and `reason` (why not), `log` (the end of `update.log`), `log_modified` (Unix
    seconds or `None`), `stale_request` and, when no update is waiting or running, `outcome` (`failed`, `done` or
    `None`) with `outcome_detail`: how the last one ended, from the service's state and the end of the log. `running` is the update guard's lock held by a live process, or the site
    still closed by an update.
    """
    on = installed()
    age = _request_age()
    updating = updateguard.is_updating()
    if not on:
        state = "not_installed"
    elif updating:
        state = "running"
    elif age is not None:
        state = "requested"
    else:
        state = "idle"
    reason = ""
    if state == "not_installed":
        reason = ("The Update button is not set up on this server: run ./update.sh once as root, which installs "
                  "planetgen-update.path.")
    elif age is not None:
        reason = "An update has been requested and is about to start." if not updating else \
            "An update has been requested and will run after the one that is running."
    else:
        reason = maintenance.write_problem() or ""
    text, modified = log_tail() if on else ("", None)
    outcome, outcome_detail = (None, "")
    if on and state == "idle":
        outcome, outcome_detail = _outcome(text, service_status())
    return {
        "outcome": outcome, "outcome_detail": outcome_detail,
        "state": state, "installed": on, "can_request": not reason, "reason": reason, "log": text,
        "log_modified": modified, "stale_request": age is not None and not updating and age > STALE_REQUEST_SECONDS,
    }


def request():
    """
    Asks root to run `update.sh`: makes the empty request file. The caller has already checked who may ask (the Owner,
    with a recent password). A press while an update runs is accepted and runs after it.

    Raises:
        UpdateRequestError: Not set up here, already requested, or the state folder cannot be written.
    """
    current = view()
    if not current["can_request"]:
        raise UpdateRequestError(current["reason"])
    try:
        fd = os.open(request_path(), os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0), 0o660)
        os.close(fd)
    except FileExistsError as exc:
        raise UpdateRequestError("An update has been requested and is about to start.") from exc
    except OSError as exc:
        raise UpdateRequestError(f"The request could not be written ({exc.strerror or exc}); is the state folder "
                                 f"{maintenance.state_dir()} writable by the web server's user?") from exc
