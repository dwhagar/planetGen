# stellarObjects/progressFile.py

"""
Machine-readable progress for a run started from the web interface.

The admin Generate page (`html/web/generate_page.py`) starts `generate.py`
as a background job and sets `PLANETGEN_PROGRESS_FILE` to a path inside
that job's directory. While it is set, `report` writes the run's current
progress there as a small JSON object, which the page reads to draw its
progress bar:

    {"description": "Sectors", "completed": 12, "total": 40, "updated_at": 1759236000.0}

`total` is `null` when the amount of work isn't known up front (the
`plan` skeleton scan stops when it finds the galaxy's edge). Writes are
atomic (a temp file renamed over the old one) and throttled to a few per
second, so a fast loop costs nothing. Without the variable (every
terminal run) `report` does nothing.
"""

import json
import os
import tempfile
import time

ENV_VAR = "PLANETGEN_PROGRESS_FILE"
"""str: The environment variable naming the progress file."""

MIN_INTERVAL_SECONDS = 0.5
"""float: Least time between two writes, unless `force` is given."""

_last_write = 0.0


def report(completed, total=None, description=None, force=False):
    """
    Writes the current progress to `$PLANETGEN_PROGRESS_FILE`, if set.
    Never raises: progress is a nicety, never a reason to fail a run.

    Args:
        completed (float): Work done so far.
        total (float, optional): Total work, or `None` if unknown.
        description (str, optional): What is being counted ("Sectors").
        force (bool): Write even if the last write was very recent (use
            for the final count).
    """
    global _last_write
    path = os.environ.get(ENV_VAR)
    if not path:
        return
    now = time.time()
    if not force and now - _last_write < MIN_INTERVAL_SECONDS:
        return
    _last_write = now
    body = {
        "description": description,
        "completed": completed,
        "total": total,
        "updated_at": now,
    }
    try:
        directory = os.path.dirname(os.path.abspath(path))
        fd, tmp = tempfile.mkstemp(prefix=".progress-", dir=directory)
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(body, f)
        os.replace(tmp, path)
    except OSError:
        pass
