# planetgen/queue/unitdebug.py

"""
The debug file of a restarted unit (OPS.146; Boss, 5:29 PM EST on 11-09): a unit the heartbeat stopped as stuck and queued again (`unitstuck`, OPS.140) is the one whose second attempt
has to be understood, so, with `queue.debug_restarted_units` on (the default), that attempt runs in debug mode and writes a file of its own.

* **The file** is `unit-debug-<job id>-attempt<N>.log` in the log folder (`logpaths.log_dir_path`); the job id (`units.job_id`) is the run, the stage and the unit, `N` the attempt (2, the restart). It is separate from
  `error.log` and from the run's log. It is bounded: the program rotates it at `ROTATE_BYTES` and keeps `ROTATE_COPIES` older copies (the way the activity log and `error.log` are rotated by the program where the
  system does not rotate them), its name ends in `.log` so the folder's one logrotate rule covers it as well, and the newest `KEEP_FILES` such files are kept.
* **What is in it.** At the top, the stack of every thread of the first attempt, taken when it was found stuck, and why (the follower writes this when it queues the unit again: `begin`). Then everything the
  restarted attempt logs at debug level, written by the horse (`session`): its phases, counter changes at every beat, waits for the gate, the sky lock and MySQL with how long they took, every SQL statement with its
  timing, and the memory it holds at every beat. The random draws are left out (`AREAS`): they are many times the rest.
* **The run's log** gets one line naming the file (from `begin`).

With the setting off nothing is written: no file, no line, no debug mode in the horse.
"""

import contextlib
import glob
import os
import time

from planetgen.util import log, logpaths

ROTATE_BYTES = 20 * 1024 * 1024
"""int: The size at which the program rotates a unit's debug file."""

ROTATE_COPIES = 3
"""int: The older copies of a unit's debug file that rotation keeps."""

KEEP_FILES = 40
"""int: The most debug files of restarted units kept in the log folder (the oldest go when a new one is made)."""

FILE_MODE = 0o660
"""int: The mode of a debug file when it is created (the same as `error.log`'s)."""

AREAS = "-draw"
"""str: The debug areas a restarted unit records (`log.parse_areas`): all but the random draws (and `branch`, which is off unless named)."""

NAME_PREFIX = "unit-debug-"


def enabled():
    """Whether a restarted unit is run in debug mode (`queue.debug_restarted_units`)."""
    from planetgen.util.settings import get_settings
    return bool(get_settings().queue.debug_restarted_units)


def file_name(job_id, attempt):
    """The name of the debug file of attempt `attempt` of the unit whose job is `job_id`."""
    return f"{NAME_PREFIX}{job_id}-attempt{attempt}.log"


def file_path(job_id, attempt, config=None):
    """The debug file of that attempt, in the log folder."""
    return os.path.join(os.path.abspath(logpaths.log_dir_path(config)), file_name(job_id, attempt))


def prune(folder, keep=KEEP_FILES):
    """Removes the oldest debug files (and their rotated copies) of `folder` beyond the newest `keep`. Never raises."""
    try:
        found = sorted((path for path in glob.glob(os.path.join(folder, f"{NAME_PREFIX}*.log")) if os.path.isfile(path)), key=os.path.getmtime, reverse=True)
        for path in found[keep:]:
            for old in [path, *glob.glob(path + ".*")]:
                with contextlib.suppress(OSError):
                    os.remove(old)
    except OSError:
        pass


def begin(job_id, attempt, unit_label, reason, first_stack, now=None):
    """
    The follower's side, when it queues the stopped unit again: makes the debug file of the new attempt with why the first one was stopped and its stack, and returns its path, or `None` when the setting is off
    or the file cannot be made (the folder does not admit this user), which costs the unit nothing.
    """
    if not enabled():
        return None
    try:
        path = file_path(job_id, attempt)
        os.makedirs(os.path.dirname(path), mode=0o770, exist_ok=True)
        prune(os.path.dirname(path), KEEP_FILES - 1)
        stamp = time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(time.time() if now is None else now))
        fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, FILE_MODE)
        with os.fdopen(fd, "a", encoding="utf-8") as handle:
            handle.write(f"{stamp} debug file of {unit_label}, attempt {attempt}\n{stamp} the first attempt was found stuck and stopped: {reason}\n"
                         f"{stamp} the stack of every thread of the first attempt when it was found stuck:\n{(first_stack or '(none)').rstrip()}\n"
                         f"{stamp} ---- attempt {attempt} begins; what it logs follows ----\n")
        with contextlib.suppress(OSError):
            os.chmod(path, FILE_MODE)
        return path
    except OSError:
        return None


@contextlib.contextmanager
def session(path):
    """
    The horse's side: while the block runs, the attempt's debug output goes to `path` (`log.debug_to_file`, in `AREAS`, bounded by `ROTATE_BYTES` and `ROTATE_COPIES`). Nothing is recorded, and the block still runs,
    where `path` is empty or cannot be opened.
    """
    with contextlib.ExitStack() as stack:
        if path:
            try:
                stack.enter_context(log.debug_to_file(path, AREAS, rotate=(ROTATE_BYTES, ROTATE_COPIES), create_mode=FILE_MODE))
            except OSError:
                path = None
        yield path
