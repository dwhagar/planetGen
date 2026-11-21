"""
Whether a process is alive on this machine: the one answer the queue (workers, the run line, the stop-all, the memory gate,
the stuck-unit watch) and the admin's update guard share, so a fix to it lands everywhere at once.
"""

import os
import signal


def pid_alive(pid):
    """
    Whether a process with this id is running here. A process of another user counts as running (the signal check
    says permission denied, so it is there); a zombie (ended, not yet reaped) does not; anything that is not a
    positive int names nothing.
    """
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
            return handle.read().rpartition(b")")[2].split()[0] != b"Z"
    except (OSError, IndexError):
        return True


class Stopped(KeyboardInterrupt):
    """A command was told to stop by SIGTERM; `reason` says by what (the line a canceller left, OPS.95). A `KeyboardInterrupt`, so a
    run ends as canceled as Ctrl+C ends it."""

    def __init__(self, reason):
        super().__init__(reason)
        self.reason = reason


NO_REASON = "a SIGTERM that left no reason (a kill, or a cancel that could not write one)"
"""str: What a command says when its SIGTERM came with no reason left for it."""


def interrupt_on_sigterm():
    """Makes SIGTERM (Cancel on the Generate page, `kill`, the update guard, stop-all) raise `Stopped`, a `KeyboardInterrupt`, so a
    run ends as canceled in Redis the way Ctrl+C ends it, and carries who sent the signal when the canceller left a reason
    (`admin.stop_reason`). Returns the handler it replaced, for `signal.signal(signal.SIGTERM, previous)`."""
    from planetgen.admin import stop_reason
    stop_reason.clear(os.getpid())

    def on_term(_signum, _frame):
        raise Stopped(stop_reason.take(os.getpid()) or NO_REASON)
    return signal.signal(signal.SIGTERM, on_term)
