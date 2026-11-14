# planetgen/db/statementmark.py

"""
Which statement this process is waiting on, and since when (OPS.193; reports/job-engine-items-stage-1-2-2026-11-13.md, item P3).

MySQL's own progress figures are not readable by the account INSTALL.md creates (`performance_schema` is denied), so a job reports its own statement from inside: the database layer
(`db/store.py`) marks the start and the end of every statement it sends, and the beat thread of a unit (`queue/unitbeat.py`) reads the mark every 10 s. A statement that has run for
minutes then shows as "in SQL since T: label" rather than as a quiet unit.

* **The mark costs one dictionary write and one removal** per statement and nothing else: the text is only shortened into a label when the beat asks (`running`).
* **A label can be given** (`labelled`, which `ctx.sql` of OPS.195 uses), and then it names the statements inside the block instead of their text.
* **One entry per thread**, because a unit has its work thread and, at times, another that sends a statement; `running` returns the one that has waited longest.
"""

import threading
import time
from contextlib import contextmanager

LABEL_CHARS = 80
"""int: The longest label shown for a statement's own text."""

_marks = {}         # thread id -> (text, label, progress function, monotonic seconds at the start)
_given = {}         # thread id -> (the label a `labelled` block gave, the function that reads its progress or `None`)
clock = time.monotonic


class statement:
    """Marks the statement `sql` as running in this thread for the block (a plain class: it runs once for every statement the program sends)."""

    __slots__ = ("_sql", "_ident", "_before")

    def __init__(self, sql):
        self._sql = sql

    def __enter__(self):
        self._ident = ident = threading.get_ident()
        self._before = _marks.get(ident)
        _marks[ident] = (self._sql, *_given.get(ident, (None, None)), clock())

    def __exit__(self, *exc_info):
        if self._before is None:
            _marks.pop(self._ident, None)
        else:
            _marks[self._ident] = self._before
        return False


@contextmanager
def labelled(label, progress=None):
    """
    The statements this thread sends inside the block are shown as `label`. `progress`, if given, is a function the beat calls for `(done, estimated)` work of the statement now
    out (MySQL's stage events, where the account may read them: OPS.195), or `None`.
    """
    ident = threading.get_ident()
    before = _given.get(ident)
    _given[ident] = (str(label), progress)
    try:
        yield
    finally:
        if before is None:
            _given.pop(ident, None)
        else:
            _given[ident] = before


def shorten(sql):
    """One line of at most `LABEL_CHARS` characters from a statement's text."""
    text = " ".join(str(sql).split())
    return text if len(text) <= LABEL_CHARS else text[:LABEL_CHARS - 3] + "..."


def running():
    """`(label, seconds, progress)` of the statement that has run longest in this process (`progress` is a function or `None`), or `None` when no thread is in one."""
    now = clock()
    oldest = None
    for mark in list(_marks.values()):
        if oldest is None or mark[3] < oldest[3]:
            oldest = mark
    if oldest is None:
        return None
    return (oldest[1] or shorten(oldest[0]), max(now - oldest[3], 0.0), oldest[2])


class Steps:
    """
    How a clear or a drop names its statements (OPS.200): `with steps.sql("Emptying bright_stars"): conn.execute(...)` labels the statements of the block for the beat. This one only labels; a
    unit gives its own (`run_plan._UnitPhases`), which also says the step begins, so the unit's bar counts one step to a statement.
    """

    def sql(self, label):
        return labelled(label)


LABELS_ONLY = Steps()
"""Steps: What the clears and drops of `db/store.py` use where the caller gives none."""
