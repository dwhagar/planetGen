"""
Whether a process is alive on this machine: the one answer the queue (workers, the run line, the stop-all, the memory gate,
the stuck-unit watch) and the admin's update guard share, so a fix to it lands everywhere at once.
"""

import os


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
