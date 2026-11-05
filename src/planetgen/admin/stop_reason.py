# planetgen/admin/stop_reason.py

"""
Why a command-line run was told to stop (OPS.95).

A run that ends on SIGTERM cannot tell by itself who sent it: update.sh's time limit, the update guard's cancel (an update
starting, the Stats page's cancel) and stop-all all use the same signal. A canceller leaves one short line for the pid it is
about to signal in `stop-reasons` in the state folder (`leave`), and the run reads it when the signal arrives (`take`), so
it can print one line saying what ended it instead of a traceback. A folder that cannot be written only means the run says
less. A pid's file is removed by the run when it starts, so a pid used again does not read an old reason.
"""

import contextlib
import os

from planetgen.admin import maintenance

FOLDER = "stop-reasons"
"""str: The folder of one-line reason files, named by pid, in the state folder."""

MAX_CHARS = 200
"""int: A reason is cut to this many characters."""


def folder():
    """The folder the reasons are left in."""
    return os.path.join(maintenance.state_dir(), FOLDER)


def _path(pid):
    return os.path.join(folder(), f"{int(pid)}.txt")


def leave(pid, reason):
    """Leaves `reason` (one line) for process `pid` to read; best effort, returns whether it was written."""
    text = " ".join(str(reason).split())[:MAX_CHARS]
    if not text:
        return False
    try:
        os.makedirs(folder(), exist_ok=True)
        with contextlib.suppress(OSError):
            os.chmod(folder(), 0o1777)      # every user a cancel runs as may leave one, and every run may read it
        path = _path(pid)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")
        os.chmod(path, 0o644)
    except OSError:
        return False
    return True


def take(pid):
    """The reason left for process `pid`, or `None`; the file is removed."""
    path = _path(pid)
    try:
        with open(path, "r", encoding="utf-8") as handle:
            text = handle.readline().strip()
    except OSError:
        return None
    with contextlib.suppress(OSError):
        os.unlink(path)
    return text or None


def clear(pid):
    """Removes any reason left for `pid` before it ran (the pid may be a reused one)."""
    with contextlib.suppress(OSError):
        os.unlink(_path(pid))
