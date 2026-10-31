# planetgen/admin/cachebuild.py

"""
Whether the site's caches are being rebuilt, and which requests that rules out.

The API runs as one mod_wsgi daemon with a 60-second request limit (`examples/apache/planetgen.conf.example`). A
request that has to build a big cache itself (the Galaxy Map's nodes, the shape, the lists, the search) runs far past
it on a cold cache, mod_wsgi kills the daemon, and every request in flight is cut off ("Truncated or oversized
response headers"). So while the caches are being rebuilt, such a request that the cache cannot answer is refused with
a 503 and `Retry-After` (`web/rebuild_gate.py`) and the build is left to the warm-up. Requests the cache answers, and
cheap ones (one sector, system or phenomenon by id), go on as usual.

The state needs no database and no Redis; it is read from files in the state folder, so it works while the database
is being migrated:

- maintenance mode is on, whatever its cause (the site is closed because something is being worked on, so nothing
  that could endanger its stability is started); or
- an update runs (`updateguard.is_updating`: its lock is held by a live process, or the site is still closed with the
  update's cause); or
- a generation run is live in the run record (`queue/runrecord.py`): it keeps changing the rows the galaxy's shape
  is built from, so every process that needs the shape builds it again from scratch (56 to 65 s on the test server,
  past the daemon's limit); or
- a cache warm-up asked for by hand or by `update.sh` runs (`planetgen warm` registers `warm-map` in the state
  folder's `running/` with `refuse: true`; a warm-up that only follows a change in the data does not).

A lock or marker left by a process that died is ignored (its pid is gone), and so is one older than
`MAX_REBUILD_SECONDS`, which is what keeps a reused pid from closing the site for good.
"""

import contextvars
import json
import os
import time
from datetime import datetime, timezone

from planetgen.admin import maintenance, updateguard

MAX_REBUILD_SECONDS = 12 * 3600
"""int: A warm-up marker older than this is ignored."""

RETRY_AFTER_SECONDS = 60
"""int: What `Retry-After` says, and how often the "rebuilding" page asks again."""

MESSAGE = ("The site is rebuilding its caches (an update or a generation run is under way), and this page needs one "
           "that is not ready yet. "
           "Nothing is wrong; try again in a minute or two.")
"""str: What a refused visitor is told."""

HEAVY_EXACT = frozenset({"/search", "/phenomena", "/systems", "/sectors", "/uncharted-systems"})
HEAVY_PREFIXES = ("/galaxy/",)
CHEAP = frozenset({"/galaxy/stamp", "/galaxy/changes"})
"""The API paths (without `/api`) whose answer is built from a cache that takes time to rebuild. Everything else (an
object by id, the sign-in, the health check) is cheap and always served."""


class CacheRebuilding(RuntimeError):
    """A request that needs a cache that is being rebuilt; `state` says why (`state()`)."""

    def __init__(self, state, what=None):
        super().__init__(MESSAGE)
        self.state, self.what = state, what


def is_heavy(path):
    """Whether a GET of API path `path` (with or without `/api`, without its query) needs a rebuilt cache."""
    path = path.split("?", 1)[0]
    if path.startswith("/api/") or path == "/api":
        path = path[4:]
    path = path.rstrip("/") or "/"
    if path in CHEAP:
        return False
    return path in HEAVY_EXACT or path.startswith(HEAVY_PREFIXES)


def _live_warmups():
    found = []
    folder = updateguard.running_dir()
    try:
        names = os.listdir(folder)
    except OSError:
        return found
    for name in names:
        if not name.startswith("warm-map-"):
            continue
        try:
            with open(os.path.join(folder, name), "r", encoding="utf-8") as handle:
                data = json.load(handle)
        except (OSError, ValueError):
            continue
        if not isinstance(data, dict) or not data.get("refuse") or not updateguard.pid_alive(data.get("pid")):
            continue
        try:
            since = datetime.strptime(data["since"], "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc).timestamp()
        except (KeyError, ValueError, TypeError):
            since = None
        if since is not None and time.time() - since > MAX_REBUILD_SECONDS:
            continue
        found.append(data)
    return found


_RUN_MEMO = {"at": -1e9, "since": None}
RUN_MEMO_SECONDS = 5.0
"""float: How long the answer about live generation runs is kept (it asks Redis)."""


def _generation_since():
    """When the oldest live generation run started (ISO text), `""` when it has none, `None` when it cannot be told.
    Redis that does not answer means no run is known: the guard then does not refuse on this account."""
    now = time.monotonic()
    if now - _RUN_MEMO["at"] > RUN_MEMO_SECONDS:
        since = None
        try:
            from planetgen.queue import runrecord
            for item in updateguard.RunRecordSource().items():
                if item["kind"] == runrecord.GENERATION:
                    since = min(since or item["started"] or "", item["started"] or "")
        except Exception:  # noqa: BLE001
            since = None
        _RUN_MEMO.update(at=now, since=since)
    return _RUN_MEMO["since"]


def state():
    """`{"cause": "update" | "generation" | "warm", "since": ISO time or None}` while the caches are being rebuilt, else `None`."""
    lock = updateguard.read_lock()
    closed = maintenance.read()
    if closed["enabled"] and closed["cause"] != updateguard.UPDATE_CAUSE and lock is None:
        return {"cause": "maintenance", "since": closed.get("since")}
    if lock is not None or updateguard.is_updating():
        return {"cause": "update", "since": (lock or {}).get("since")}
    running = _generation_since()
    if running is not None:
        return {"cause": "generation", "since": running or None}
    warmups = _live_warmups()
    if warmups:
        return {"cause": "warm", "since": min(item.get("since") or "" for item in warmups) or None}
    return None


visitor = contextvars.ContextVar("planetgen_cachebuild_visitor", default=False)
"""ContextVar: Whether the code running is answering a visitor's request (set by `web/rebuild_gate.gate` for every
request, `False` for the warm-up's app). A command, a queue worker or the warm-up never sets it."""


def refuse_expensive_build(what):
    """
    Called by code about to build, from scratch and in this process, something that takes a minute or more (the
    galaxy's cluster field, `db/store._component_field`). Raises `CacheRebuilding` when it would be a visitor's
    request that builds it (admins too: the daemon's limit does not care who asked) while the caches are being
    rebuilt or maintenance is on. Anywhere else (the warm-up, a worker, the command line) it does nothing.
    """
    if visitor.get():
        current = state()
        if current is not None:
            raise CacheRebuilding(current, what)


def refuse_if_rebuilding(path):
    """Raises `CacheRebuilding` when `path` (an API path) is heavy and the caches are being rebuilt."""
    if is_heavy(path):
        current = state()
        if current is not None:
            raise CacheRebuilding(current)
