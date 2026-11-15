# planetgen/queue/attemptlog.py

"""
The files an attempt leaves (OPS.223, L9 of `docs/design/operations/job-engine.md`; absorbs OPS.146).

An attempt's row in the ledger (`engine_attempts`) says how it ended and what it cost; two kinds of file hold what does not fit in a row, both in the log folder's
`unit-debug` folder, named by the unit's ledger job, its unit and its attempt:

* **The stack of an attempt that failed or was stopped** (`job<job>-attempt<n>-stack.txt`, `keep_stack`): the traceback of a failure, or the stack of every thread the follower
  had the horse write (`unitstuck`, OPS.140) before it stopped a stuck unit. The row's `stack_file` names it, so a reader of the Queue page can open it after the horse and its
  temporary stack file are gone.
* **The debug log of a restarted attempt** (`job<job>-unit<unit>-attempt<n>.log`, `capture`): the second and later attempt at a unit (a retry after a kill, a stuck stop or a failure,
  a resume after an update) run in debug mode with every debug area on, whatever the process's own settings, written to a file of their own that opens with how the attempt before
  it ended and its stack. The setting `queue.debug_restarts` (on) turns it off; with it off nothing extra is written. The run log gets one line naming the file.

Every file is bounded: a debug log is rotated at `FILE_BYTES` (one older copy), and the folder keeps the newest `KEEP_FILES` files, oldest removed as new ones are made.
"""

import contextlib
import os
import time

from planetgen.util import log, logpaths

FOLDER = "unit-debug"
"""str: The folder in the log folder that holds the files."""

FILE_BYTES = 20 * 1024 * 1024
"""int: The size at which a restart's debug log is rotated."""

FILE_COPIES = 1
"""int: How many rotated copies of a debug log are kept."""

KEEP_FILES = 400
"""int: The most files kept in the folder."""

STACK_CHARS = 200_000
"""int: The most of a stack that is kept."""

HEADER_STACK_CHARS = 12_000
"""int: How much of the earlier attempt's stack opens the debug log (the whole of it is in its file)."""


def folder():
    """The folder of the files: `unit-debug` in the log folder."""
    return os.path.join(os.path.abspath(logpaths.log_dir_path()), FOLDER)


def debug_path(row, unit, number):
    """The debug log of attempt `number` of the unit whose ledger job is `row`."""
    return os.path.join(folder(), f"job{row}-unit{unit}-attempt{number}.log")


def stack_path(row, number):
    """The stack file of attempt `number` of the ledger job `row`."""
    return os.path.join(folder(), f"job{row}-attempt{number}-stack.txt")


def enabled():
    """Whether a restarted attempt writes a debug log (`queue.debug_restarts`)."""
    from planetgen.util import settings
    return bool(settings.get_settings().queue.debug_restarts)


def _prune(keep=None):
    """Removes the oldest files of the folder past `keep` (`KEEP_FILES`)."""
    keep = KEEP_FILES if keep is None else keep
    root = folder()
    try:
        names = [os.path.join(root, name) for name in os.listdir(root)]
        names.sort(key=os.path.getmtime)
        for path in names[:max(0, len(names) - keep)]:
            os.unlink(path)
    except OSError as exc:
        log.debug(f"The unit debug folder {root} was not pruned ({exc}).", area="queue")


def keep_stack(row, number, text):
    """
    Writes `text` (the stack or traceback of an attempt that did not finish) to the attempt's stack file and returns its path, or `None` when there is nothing to write or the
    file cannot be written (said in the debug log: an attempt's history is never worth a unit).
    """
    if not text:
        return None
    path = stack_path(row, number)
    try:
        os.makedirs(os.path.dirname(path), mode=0o770, exist_ok=True)
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o660)
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text[:STACK_CHARS])
        _prune()
        return path
    except OSError as exc:
        log.debug(f"The stack of attempt {number} of job {row} was not kept ({exc}).", area="queue")
        return None


def _read(path, limit):
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            return handle.read(limit)
    except (OSError, TypeError):
        return ""


def _header(attempt, run_id, stage, unit, kind, before):
    lines = [f"=== debug log of attempt {attempt.number} of unit {unit} of stage {stage} of run {run_id} ({kind}), ledger job {attempt.row}, pid {os.getpid()}, "
             f"started {time.strftime('%Y-%m-%d %H:%M:%S %Z')}"]
    if before:
        ended = before["ended_at"].strftime("%Y-%m-%d %H:%M:%S") if before.get("ended_at") else "never (nothing wrote its end)"
        lines.append(f"=== attempt {attempt.number - 1} ended {before['outcome']} at {ended}: {before.get('error_kind') or ''} {before.get('error_text') or ''}".rstrip())
        stack = _read(before.get("stack_file"), HEADER_STACK_CHARS)
        if stack:
            lines.append(f"=== its stack ({before['stack_file']}):\n{stack}")
    lines.append("=== the debug lines of this attempt follow")
    return "\n".join(lines)


@contextlib.contextmanager
def capture(conn, attempt, run_id, stage, unit, kind):
    """
    The debug log of a restarted attempt (`attempt`: `attempts.Attempt`, or `None` for a unit with no ledger row) while the block runs, when the setting is on and the attempt is the
    second or a later one; a first attempt, a unit with no row and a setting turned off write nothing. Says in the run's log (`conn`, the run's Redis) where the file is.
    """
    if attempt is None or attempt.number < 2 or not enabled():
        yield None
        return
    path = debug_path(attempt.row, unit, attempt.number)
    try:
        os.makedirs(os.path.dirname(path), mode=0o770, exist_ok=True)
        _prune()
        before = attempt.previous()
    except OSError as exc:
        log.debug(f"The debug folder for unit {unit} was not made ({exc}).", area="queue")
        yield None
        return
    with log.captured(path, rotate=(FILE_BYTES, FILE_COPIES), header=_header(attempt, run_id, stage, unit, kind, before)) as written:
        if written:
            from planetgen.queue import runrecord
            runrecord.log(conn, run_id, f"unit {unit} is at attempt {attempt.number}: its debug log is {written}", "info", "unit-debug-file", stage=stage, unit=unit)
        yield written
