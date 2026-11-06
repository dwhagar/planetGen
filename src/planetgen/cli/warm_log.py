# planetgen.cli.warm_log

"""
The warm-up log (OPS.96): `warmup.log` in the log folder (`logpaths.warmup_log_path`)
=====================================================================================

Every `planetgen warm` that runs a warm-up appends it to `warmup.log` itself, in Python: the update's two stages, a run by hand and
the background run that starts after a change (`warmstart`) alike, so the record does not depend on how the run was started (it
used to be a `tee` into one /tmp file that stage 2 overwrote, and a background run left nothing). Each run is bracketed by a line that
says what started it and a line that says how it ended and after how long; between them are the lines the run printed, each in the
logging layer's display form (`log.LineFormatter`): the local time with its zone, the severity, the component and the message.

    2026-11-05 07:44:10 EST INFO     warm: Cache warm-up started (update.sh, part 2 of 2, the full warm-up; time limit 4 hours)
    2026-11-05 08:01:54 EST INFO     warm: Cache warm-up finished after 17:44 (5 families)

What started the run is `PLANETGEN_WARM_STARTED_BY`, which `scripts/deploy-common.sh` sets for each of its stages; a run after a
change names the change (`warmstart` hands it as the reason); any other run is "by hand". `--follow` writes nothing (a watcher is not the
run). A file that cannot be written (a hand run by a user the folder does not admit) is said once on the screen and never stops the
warm-up. The folder is rotated by the rule `examples/apache/setup-debug-log.sh` installs for `*.log` there.
"""

import os

from planetgen.util import joblog, logpaths

STARTED_BY_ENV = "PLANETGEN_WARM_STARTED_BY"
"""str: What `scripts/deploy-common.sh` tells a stage's run about itself, for the line that opens the run in `warmup.log`."""

TIME_LIMIT_ENV = "PLANETGEN_WARM_TIME_LIMIT_SECONDS"
"""str: update.sh's time limit for the stage, in seconds (the same variable `planetgen.cli.warm` reads)."""

COMPONENT = "warm"
"""str: The component every line of the file carries, whichever module printed it."""


def _words(seconds):
    """A time limit in words for the opening line: `4 hours`, `90 minutes`, `45 seconds`."""
    if seconds >= 3600 and seconds % 3600 == 0:
        count, unit = seconds // 3600, "hour"
    elif seconds >= 60 and seconds % 60 == 0:
        count, unit = seconds // 60, "minute"
    else:
        count, unit = seconds, "second"
    return f"{count} {unit}{'' if count == 1 else 's'}"


def started_by(reason=None, environ=None):
    """What started the run, as the opening line says it: the change that started it (`reason`), else what the update told it
    (`STARTED_BY_ENV`), else "by hand"; the update's time limit follows when it set one."""
    environ = os.environ if environ is None else environ
    told = (environ.get(STARTED_BY_ENV) or "").strip()
    who = f"after a change: {reason}" if reason else told or "by hand"
    try:
        limit = int(float(environ.get(TIME_LIMIT_ENV, "")))
    except ValueError:
        limit = 0
    return f"{who}; time limit {_words(limit)}" if limit > 0 and not reason else who


class RunLog(joblog.JobLog):
    """
    The record of one warm-up in `warmup.log` (`joblog.JobLog`, titled "Cache warm-up"): the line that opens the run says what started it
    (`started_by`), `say(line)` adds a line that was printed to the screen only (the bars' console), and the closing line says how it ended
    and how many families it warmed (`families`).
    """

    def __init__(self, reason=None, path=None, environ=None):
        super().__init__(logpaths.WARMUP_LOG_NAME, COMPONENT, "Cache warm-up", path=path, words="warm-up log")
        self.reason = reason
        self.environ = environ
        self.families = 0

    def __enter__(self):
        self.started_by = started_by(self.reason, self.environ)
        return super().__enter__()

    def __exit__(self, exc_type, exc, tb):
        self.detail = f"{self.families} famil{'y' if self.families == 1 else 'ies'}"
        return super().__exit__(exc_type, exc, tb)
