# planetgen/admin/stop_reason.py

"""
Why a command-line run was told to stop (OPS.95).

A run that ends on SIGTERM cannot tell by itself who sent it: update.sh's time limit, the update guard's cancel (an update
starting, the Stats page's cancel) and stop-all all use the same signal. A canceller leaves one short line for the pid it is
about to signal in `stop-reasons` in the state folder (`leave`), and the run reads it when the signal arrives (`take`), so
it can print one line saying what ended it instead of a traceback. A folder that cannot be written only means the run says
less. A pid's file is removed by the run when it starts, so a pid used again does not read an old reason.

The folder is open to every user (mode 1777, so each user a cancel runs as may leave a reason), and the canceller is often
root, so nothing in it is trusted (SEC.38): a reason is written to a name that does not exist yet (`O_CREAT | O_EXCL |
O_NOFOLLOW`, after removing what was there, so a link a local user planted under `<pid>.txt` is replaced and never followed),
and a reason is read only from a plain file whose owner is root or the reader (a canceller can only signal a process of its
own user or, as root, any), never through a link.
"""

import contextlib
import os
import stat

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
        with open_shared_file(folder(), _path(pid)) as handle:
            handle.write(text + "\n")
    except OSError:
        return False
    return True


def open_shared_file(directory, path, mode=0o644):
    """Opens `path` in the world-open `directory` for writing, the way a file every user may leave and every user may read
    is written safely (SEC.38): the folder is made (mode 1777), whatever sits under the name is removed (a link a local
    user planted is removed itself, never followed), the file is created new (`O_CREAT | O_EXCL | O_NOFOLLOW`) and given
    `mode` whatever the umask. Raises `OSError` when any of that fails; the caller decides what that means."""
    os.makedirs(directory, exist_ok=True)
    with contextlib.suppress(OSError):
        os.chmod(directory, 0o1777)
    with contextlib.suppress(FileNotFoundError):
        os.unlink(path)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, mode)
    handle = os.fdopen(descriptor, "w", encoding="utf-8")
    os.fchmod(handle.fileno(), mode)
    return handle


def take(pid):
    """The reason left for process `pid`, or `None`; the file is removed."""
    path = _path(pid)
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        with os.fdopen(descriptor, "r", encoding="utf-8") as handle:
            status = os.fstat(handle.fileno())
            if not stat.S_ISREG(status.st_mode) or status.st_uid not in (0, os.getuid()):
                text = ""       # not a reason a canceller left: ignored, and removed below
            else:
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
