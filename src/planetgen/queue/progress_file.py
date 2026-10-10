# planetgen/queue/progress_file.py

"""
Machine-readable progress for a run started from the web interface.

The admin Generate page (`planetgen/web/generate_page.py`) starts `planetgen`
as a background job and sets `PLANETGEN_PROGRESS_FILE` to a path inside
that job's directory. While it is set, `report` writes the run's current
progress there as a small JSON object, which the page reads to draw its
progress bar:

    {"description": "Sectors", "completed": 12, "total": 40, "rate": 0.8,
     "eta_s": 35.0, "detail": null, "percent": false, "updated_at": 1759236000.0}

`total` is `null` when the amount of work isn't known up front; `rate`
(units per second) and `eta_s` (seconds left, as of `updated_at`) are
`null` until the first unit is done. `detail` is a second bar under the
first, or `null`: `{"description", "completed", "total", "eta_s"}` (the
bright-star layers being drawn, while layers are slow); `percent` asks for
the bar as a share done rather than counts. Writes are
atomic (a temp file renamed over the old one) and throttled to a few per
second, so a fast loop costs nothing. Without the variable (every
terminal run) `report` does nothing.
"""

import json
import math
import os
import tempfile
import time

ENV_VAR = "PLANETGEN_PROGRESS_FILE"
"""str: The environment variable naming the progress file."""

MIN_INTERVAL_SECONDS = 0.5
"""float: Least time between two writes, unless `force` is given."""

_last_write = 0.0
_stage = None
"""dict | None: The stage the run is in (`set_stage`), written into every report."""


def report(completed, total=None, description=None, force=False, rate=None, eta_s=None, detail=None, percent=False):
    """
    Writes the current progress to `$PLANETGEN_PROGRESS_FILE`, if set.
    Never raises: progress is a nicety, never a reason to fail a run.

    Args:
        completed (float): Work done so far.
        total (float, optional): Total work, or `None` if unknown.
        description (str, optional): What is being counted ("Sectors").
        force (bool): Write even if the last write was very recent (use
            for the final count).
        rate (float, optional): Units finished per second, the decaying
            average the terminal's ETA uses (`progressRate`).
        eta_s (float, optional): Seconds left at that rate, as of
            `updated_at`.
        detail (dict, optional): A second bar under this one (PERF.4,
            the bright-star layers being drawn while layers are slow):
            `description`, `completed`, `total` and `eta_s`.
        percent (bool): Show this bar as a share done, not as counts
            (PERF.9's weighted bright-star bar).
    """
    global _last_write
    path = os.environ.get(ENV_VAR)
    if not path:
        return
    now = time.time()
    if not force and now - _last_write < MIN_INTERVAL_SECONDS:
        return
    _last_write = now
    if percent and _finite_or_none(completed) is not None and _finite_or_none(total) is not None:
        # A share is never shown past 100% (PERF.23).
        completed = min(completed, total)
    body = {
        "description": description,
        "completed": completed,
        "total": total,
        "rate": _finite_or_none(rate),
        "eta_s": _finite_or_none(eta_s),
        "detail": detail,
        "percent": bool(percent),
        "stage": _stage,
        "updated_at": now,
    }
    try:
        # Serialized up front, so a value json can't encode fails before
        # any temp file exists.
        text = json.dumps(body, allow_nan=False)
        directory = os.path.dirname(os.path.abspath(path))
    except (TypeError, ValueError):
        return
    tmp = None
    try:
        fd, tmp = tempfile.mkstemp(prefix=".progress-", dir=directory)
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(text)
        os.replace(tmp, path)
        tmp = None
    except (OSError, ValueError):
        pass
    finally:
        # Only still set if the write or the rename failed (e.g. `path` is
        # a directory) -- don't leave the half-done temp file behind.
        if tmp is not None:
            try:
                os.unlink(tmp)
            except OSError:
                pass


def set_stage(index, total, label, skipped=None):
    """Records the stage the run is in (UX.89: `index` of `total`, 1-based; `skipped` is the reason a stage that
    did not run was skipped) and writes it at once, keeping the last bar's fields."""
    global _stage
    _stage = {"index": index, "total": total, "label": label, "skipped": skipped}
    path = os.environ.get(ENV_VAR)
    if not path:
        return
    try:
        with open(path, "r", encoding="utf-8") as f:
            body = json.load(f)
    except (OSError, ValueError):
        body = {"description": None, "completed": 0, "total": None, "rate": None, "eta_s": None, "detail": None,
                "percent": False}
    body["stage"] = _stage
    body["updated_at"] = time.time()
    try:
        directory = os.path.dirname(os.path.abspath(path))
        fd, tmp = tempfile.mkstemp(prefix=".progress-", dir=directory)
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(json.dumps(body, allow_nan=False))
        os.replace(tmp, path)
    except (OSError, ValueError, TypeError):
        pass


def _finite_or_none(value):
    """`value`, or `None` when it's a NaN or infinite float -- JSON has no
    such numbers, and the Generate page's browser would refuse the file."""
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value
