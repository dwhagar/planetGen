# planetgen/web/lib/slowpages.py

"""
The cost rule for slow page requests, and what a page sees of them (PERF.105; design in
docs/design/slow-page-answers.md).

A page whose answer may take a while opts in for the length of its route with `allow_pending()`. Inside it,
`apiclient` checks each request that misses the page cache against the cost rule (`cost`), once, here. A slow one
is not built in the request: `apiclient` raises `Pending`, carrying the job (`planetgen.queue.pagejobs`) that
builds the answer into the page cache, and the page shows its loading state and polls `GET /api/page-jobs/<id>`
until the answer is there. Outside `allow_pending`, and in the job itself, a miss is built in the request as it
always was.

**The cost rule.** A request is slow when, in this order:

1. its last build (`pagejobs.record_duration`, kept `tuning.SLOW_PAGE_TIMING_SECONDS`) took longer than
   `tuning.SLOW_PAGE_SECONDS`; a measured time beats any estimate, so a request that turned out quick is quick
   from then on; or
2. it has never been timed itself, and the slowest build of any request on the same route within the last
   `tuning.SLOW_PAGE_PATH_SECONDS` took longer than `tuning.SLOW_PAGE_SECONDS` (a search for a new name, just after a
   wide one); a route whose recent builds were all quick is quick for a request not seen before; or
3. nothing on the route has been timed lately and the tables its route reads in full (`PATH_TABLES`) hold more than
   `tuning.SLOW_PAGE_ROW_BUDGET` rows, by the database's own estimate (`GET /api/table-estimates`, kept a minute
   in the process).

A route that reads one page of an indexed table is not in `PATH_TABLES`, or does not call `allow_pending`.

**Redis down** (or an update running, which takes no new jobs): the answer is built in the request, as it was before
this mechanism, and the visitor is told so: `take_notice()` gives the sentence once, for the page to show.
"""

import contextlib
import contextvars
import threading
import time

from planetgen import tuning
from planetgen.util import log

PHENOMENON_TABLES = ("phenomenon_scatter", "nebulae", "asteroid_fields", "black_holes", "neutron_stars",
                     "supernova_remnants", "rogue_planets", "interstellar_comets", "quasars")
"""tuple: Every table the phenomena list reads (`db.query._PHENOMENON_TABLES` and the scatter)."""

PATH_TABLES = {
    "/phenomena": PHENOMENON_TABLES,
    "/search": ("star_systems", "bright_stars", *PHENOMENON_TABLES),
    "/uncharted-systems": ("bright_stars",),
    "/systems": ("star_systems",),
}
"""dict: Path prefix under `/api` -> the tables the route reads in full, the rows the cost rule adds up. The
longest matching prefix counts."""

_allowed = contextvars.ContextVar("planetgen_slow_pages_allowed", default=False)
_retrying = contextvars.ContextVar("planetgen_slow_pages_retrying", default=False)
_notices = threading.local()

_estimates = {}
_estimates_lock = threading.Lock()


class Pending(Exception):
    """
    The answer to a request is being built by a job; the page shows its loading state.

    Attributes:
        job (dict): `pagejobs.status` of the job (`id`, `state`, `done`, `total`, `eta_seconds` ...).
        target (str): The request that is waiting.
    """

    def __init__(self, job, target):
        super().__init__(f"{target} is being built by job {job['id']} ({job['state']})")
        self.job, self.target = job, target


@contextlib.contextmanager
def allow_pending(retry=False):
    """
    `with allow_pending():` -- inside the block, a slow request raises `Pending` instead of being built.

    Args:
        retry (bool): The visitor asked to try again: a job that failed or was cancelled is replaced by a new one at
            once, not shown again.
    """
    token, retrying = _allowed.set(True), _retrying.set(retry)
    try:
        yield
    finally:
        _retrying.reset(retrying)
        _allowed.reset(token)


def pending_allowed():
    """Whether the current route opted in with `allow_pending`."""
    return _allowed.get()


def take_notice():
    """The sentence saying a slow answer was built in the request because the queue could not take it, once; else
    `None`."""
    notice = getattr(_notices, "text", None)
    _notices.text = None
    return notice


def clear_estimates():
    """Forgets the table row estimates (a test)."""
    with _estimates_lock:
        _estimates.clear()


def _table_rows(db):
    from planetgen.web.lib import apiclient
    now = time.monotonic()
    with _estimates_lock:
        known = _estimates.get(db)
    if known is not None and now - known[0] < tuning.SLOW_PAGE_ESTIMATE_SECONDS:
        return known[1]
    try:
        rows = apiclient.get_table_estimates(db)
    except Exception as exc:  # noqa: BLE001 -- no estimate: the timed rule alone applies
        log.debug(f"Slow pages: no table estimates ({exc}).")
        return {}
    with _estimates_lock:
        _estimates[db] = (now, rows)
    return rows


def estimated_rows(db, path):
    """The rows the route `path` (under `/api`) reads in full, by the database's estimate; `None` when the route
    is not listed."""
    for prefix in sorted(PATH_TABLES, key=len, reverse=True):
        if path == prefix or path.startswith(prefix + "/") or path.startswith(prefix + "?"):
            rows = _table_rows(db)
            return sum(rows.get(table, 0) for table in PATH_TABLES[prefix])
    return None


def record_duration(db, target, seconds):
    """Remembers how long a build of `target` took, for the cost rule (`pagejobs.record_duration`)."""
    from planetgen.queue import pagejobs
    pagejobs.record_duration(db, target, seconds)


def cost(db, path, target):
    """
    Whether building `target` is slow.

    Returns:
        dict | None: `rows` (the estimate, or `None`), `expected` (seconds) and `why` (`"timed"`, `"route"` or
            `"rows"`) for a slow request; `None` for one that isn't.
    """
    from planetgen.queue import pagejobs   # not at import: a page that never opts in doesn't load the queue
    last = pagejobs.last_duration(db, target)
    rows = estimated_rows(db, path)
    if last is not None:
        return {"rows": rows, "expected": last, "why": "timed"} if last > tuning.SLOW_PAGE_SECONDS else None
    recent = pagejobs.slowest_recent_duration(db, path)
    if recent is not None:
        return {"rows": rows, "expected": recent, "why": "route"} if recent > tuning.SLOW_PAGE_SECONDS else None
    if rows is not None and rows > tuning.SLOW_PAGE_ROW_BUDGET:
        return {"rows": rows, "expected": rows / tuning.SLOW_PAGE_ROWS_PER_SECOND, "why": "rows"}
    return None


def pending_for(db, path, target):
    """
    The `Pending` to raise for a request that missed the cache, or `None` when it is built in the request: not
    slow, no database named, or the queue can't take it (the notice says so).
    """
    from planetgen.queue import pagejobs
    if db is None:
        return None
    slow = cost(db, path, target)
    if slow is None:
        return None
    try:
        job = pagejobs.ensure(db, target, slow["rows"], slow["expected"], retry=_retrying.get())
    except pagejobs.NoQueue as exc:
        _notices.text = ("This page is slow to build and the work queue is not available, so it was built while "
                         "you waited.")
        log.normal(f"Slow pages: building {path} in the request, the queue can't take it: {exc}")
        return None
    return Pending(job, target)
