# planetgen.util.joblog

"""
The named log of a long job (OPS.96, OPS.101): `warmup.log`, `tiles.log`, `orbits.log` in the log folder
====================================================================================================

A job that runs unattended (the cache warm-up, the nightly whole-map build, the nightly orbit update) appends each run to a log of
its own in the log folder (`logpaths.job_log_path`), written by Python rather than by redirecting output, so the record does not depend on
how the run was started. A run is bracketed by a line that says what started it and a line that says how it ended and after how long;
between them are the lines the run printed, each in the logging layer's display form (`log.LineFormatter`): the local time with its zone,
the severity, the component and the message.

    2026-11-05 07:44:10 EST INFO     tiles: Whole-map build started (systemd timer)
    2026-11-05 08:01:54 EST INFO     tiles: Whole-map build finished after 17:44

A file that cannot be written (a hand run by a user the folder does not admit) is said once on the screen and never stops the job. The
folder is rotated by the rule `examples/apache/setup-debug-log.sh` installs for `*.log` there.
"""

import logging
import time

from planetgen.util import log, logpaths


def clock_text(seconds):
    """`seconds` as `1:05:12` (hours when there are any), or `4:07`."""
    seconds = int(round(max(0.0, seconds)))
    hours, rest = divmod(seconds, 3600)
    minutes, secs = divmod(rest, 60)
    return f"{hours}:{minutes:02d}:{secs:02d}" if hours else f"{minutes}:{secs:02d}"


def started_by_environment():
    """What started a command: the maintenance timer's service has systemd's `INVOCATION_ID`, a person at a terminal has not."""
    import os
    return "systemd timer" if os.environ.get("INVOCATION_ID") else "by hand"


class JobLog:
    """
    The record of one run of a job in its log file, as a context manager: it opens the file, writes the line that opens the run, appends
    every line the logging layer prints meanwhile (`log.normal`, `log.error`...) and `say(line)` adds a line that was printed to the screen
    only, and closes the run with how it ended: `finished`, `ended by <reason>` (a signal or Ctrl+C), `ended with exit status N` or
    `failed: <why>`.

    Args:
        name (str): The log's file name, one of `logpaths.JOB_LOG_NAMES` (ignored when `path` is given).
        component (str): The one component every line of the file carries (`warm`, `tiles`, `orbits`).
        title (str): What the run is, for its first and last line: `Whole-map build`.
        started_by (str): What started it, for the first line (`by hand`, `systemd timer`).
        path (str | None): The file, when it is not the log folder's `name`.
        words (str | None): The log as a person names it in the warning for a file that cannot be written (`warm-up log`); the title's, by default.
    """

    def __init__(self, name, component, title, started_by="by hand", path=None, words=None):
        self.name = name
        self.component = component
        self.title = title
        self.started_by = started_by
        self.path = path
        self.words = words or f"{title.lower()} log"
        self.detail = ""
        """str: Added to the closing line of a run that finished, in brackets (`5 families`)."""
        self.status = None
        """int | None: The exit status a job that returns one sets before it leaves the block; a non-zero one is said in the closing line."""
        self.handler = None
        self.started = time.monotonic()

    def __enter__(self):
        self.started = time.monotonic()
        try:
            path = self.path or logpaths.job_log_path(self.name)
            self.handler = log.file_handler(path, "display", self.component)
        except (OSError, TypeError, ValueError) as exc:
            log.warning(f"Cannot write the {self.words} ({exc}); the run goes on without it.")
            return self
        log.attach(self.handler)
        self._write(logging.INFO, f"{self.title} started ({self.started_by})")
        return self

    def say(self, line, level=logging.INFO):
        """Adds `line` to the file: a line that reached the screen without the logging layer."""
        self._write(level, line)

    def _write(self, level, message):
        if self.handler is None:
            return
        try:
            self.handler.handle(logging.LogRecord(f"planetgen.{self.component}", level, "", 0, message, None, None))
        except Exception:  # noqa: BLE001 -- a log that cannot be written never stops the run
            pass

    def __exit__(self, exc_type, exc, _tb):
        if self.handler is None:
            return False
        took = clock_text(time.monotonic() - self.started)
        status = exc.code if exc_type is SystemExit else self.status
        if exc_type is None or (exc_type is SystemExit and not status):
            if status:
                self._write(logging.ERROR, f"{self.title} ended with exit status {status} after {took}")
            else:
                self._write(logging.INFO, f"{self.title} finished after {took}" + (f" ({self.detail})" if self.detail else ""))
        elif exc_type is SystemExit:
            self._write(logging.ERROR, f"{self.title} ended with exit status {status} after {took}")
        elif hasattr(exc, "reason"):
            self._write(logging.ERROR, f"{self.title} ended by {exc.reason} after {took}")
        elif exc_type is KeyboardInterrupt:
            self._write(logging.ERROR, f"{self.title} ended by Ctrl+C after {took}")
        else:
            self._write(logging.ERROR, f"{self.title} failed after {took}: {exc}")
        log.detach(self.handler)
        self.handler = None
        return False
