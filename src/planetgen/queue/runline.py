# planetgen/queue/runline.py

"""
The one waiting list for work that cannot start yet (PERF.154)
==============================================================

Boss, 2026-10-22: "a task that arrives is never lost or failed only because another is running". Everything that has to
wait for the machine waits here, in one ordered list in Redis, up to `queue.max_waiting` entries (8 by default), first
come first served:

* a **run** (`type` `run`): a generation run (or an API job's, a one-sector regenerate or a neighbourhood) that finds
  the machine's `GENERATION` slot taken (`runrecord.start_run`), waiting inside `start_when_free`. It holds its place
  with a sign of life (`ALIVE_SECONDS`), so one that was killed drops out of the line by itself;
* a **job** (`type` `job`): a Generate page job that finds another job running, waiting for the Generate page's job
  slot (`web/jobs.py`). Its recipe is kept in a private file of the jobs folder (it holds the database environment);
  the line holds the rest. It starts by itself when the job before it succeeded or was cancelled and is held for an
  admin when that job failed (ADM.123).

The two kinds wait for two different things, so each keeps its own order: a run's turn comes when the runs ahead of it
have started, a job's when the jobs ahead of it have. (One order for both would let a waiting job hold up the stage of
the very job it waits for.) The strip and the Generate page read `entries`, which gives every entry in the order it
arrived with its place in its own order, 1 being next.

The line is refused, naming who is in it, when it holds `queue.max_waiting` entries (`LineFull`). A run that a Generate
job's own step starts is already admitted and always gets a place (`join(..., force=True)`).

Keys (`runrecord.PREFIX` first):

* `line`            sorted set: ticket -> arrival number;
* `line:seq`        the arrival counter;
* `line:<ticket>`   hash: `type`, `title`, `who`, `database`, `kind`, `queued_at`, `seq`, and for a job `root`,
                    `spec` (the file holding its steps and environment), `after_job` and `held`;
* `line:<ticket>:alive`  a run's sign of life.
"""

import secrets
import time

from planetgen.queue import runrecord
from planetgen.util.settings import get_settings

RUN = "run"
JOB = "job"
TYPES = (RUN, JOB)

ALIVE_SECONDS = 30
"""int: How long a waiting run's place is kept without a sign of life (it gives one at every poll)."""

POLL_SECONDS = 1.0
"""float: How often a waiting run looks again."""

KEEP_SECONDS = 30 * 86400
"""int: How long a waiting job's entry is kept if nobody starts or discards it."""

_JOIN = """
local cap = tonumber(ARGV[1])
if cap > 0 and redis.call('ZCARD', KEYS[1]) >= cap then return 0 end
local seq = redis.call('INCR', KEYS[2])
redis.call('ZADD', KEYS[1], seq, ARGV[2])
return seq
"""


class LineFull(runrecord.Busy):
    """The line already holds `queue.max_waiting` entries. A `runrecord.Busy` so that the commands that refuse to
    wait answer it like any run that may not start; the message names the line."""


def _line_key():
    return runrecord._k("line")


def _seq_key():
    return runrecord._k("line", "seq")


def _entry_key(ticket):
    return runrecord._k("line", ticket)


def _alive_key(ticket):
    return runrecord._k("line", ticket, "alive")


def max_waiting():
    """How many entries the line holds, `queue.max_waiting`."""
    return get_settings().queue.max_waiting


def _fields(conn, ticket):
    return {runrecord._s(key): runrecord._s(value) for key, value in conn.hgetall(_entry_key(ticket)).items()}


def _drop(conn, ticket):
    conn.zrem(_line_key(), ticket)
    conn.delete(_entry_key(ticket), _alive_key(ticket))


def _prune(conn):
    """Drops from the line each run whose sign of life has run out and each ticket with no entry."""
    for member in conn.zrange(_line_key(), 0, -1):
        ticket = runrecord._s(member)
        kind = runrecord._s(conn.hget(_entry_key(ticket), "type"))
        if kind is None or kind == "" or (kind == RUN and not conn.exists(_alive_key(ticket))):
            _drop(conn, ticket)


def entries(conn):
    """
    Every entry in the order it arrived, as dicts: `ticket`, `type`, `title`, `who`, `database`, `kind`, `queued_at`,
    `seq`, `position` (its place in its own order, 1 being next) and, for a job, `root`, `after_job` and `held`
    (`None` when not held).
    """
    _prune(conn)
    found, places = [], {RUN: 0, JOB: 0}
    for member, _score in conn.zrange(_line_key(), 0, -1, withscores=True):
        ticket = runrecord._s(member)
        fields = _fields(conn, ticket)
        if not fields.get("type"):
            continue
        places[fields["type"]] = places.get(fields["type"], 0) + 1
        found.append({
            "ticket": ticket, "type": fields["type"], "title": fields.get("title") or "", "who": fields.get("who") or None,
            "database": fields.get("database") or None, "kind": fields.get("kind") or None,
            "queued_at": float(fields.get("queued_at") or 0), "seq": int(fields.get("seq") or 0),
            "position": places[fields["type"]], "root": fields.get("root") or None, "spec": fields.get("spec") or None,
            "after_job": fields.get("after_job") or None, "held": fields.get("held") or None,
        })
    return found


def of_type(conn, type_, root=None):
    """The entries of one type in their order, a job's only for the jobs folder `root` when it is given."""
    return [entry for entry in entries(conn) if entry["type"] == type_ and (root is None or entry["root"] == root)]


def describe(found):
    """The line in words, for a refusal: `'"Redo scatters" (boss), run "Regenerate sector"'`."""
    return ", ".join(f'"{entry["title"]}"' + (f" ({entry['who']})" if entry["who"] else "") for entry in found)


def join(conn, type_, title, who=None, database=None, kind=None, fields=None, force=False, now=None, ticket=None):
    """
    Puts an entry at the end of the line.

    Args:
        type_ (str): `RUN` or `JOB`.
        title (str): What the pages call it.
        who (str, optional): Who asked for it.
        fields (dict, optional): More to keep with it (a job's `root`, `spec`, `after_job`).
        force (bool): Take a place even in a full line (a run that a running job's own step starts).
        ticket (str, optional): The ticket to use (a job writes its recipe under it first).

    Returns:
        str: The entry's ticket.

    Raises:
        LineFull: The line holds `queue.max_waiting` entries; the message names them.
    """
    if type_ not in TYPES:
        raise ValueError(f"an entry is one of {TYPES}")
    _prune(conn)
    ticket = ticket or secrets.token_hex(8)
    cap = 0 if force else max_waiting()
    seq = int(conn.register_script(_JOIN)(keys=[_line_key(), _seq_key()], args=[cap, ticket]))
    if not seq:
        waiting = entries(conn)
        raise LineFull(f"The line is full: {len(waiting)} {'is' if len(waiting) == 1 else 'are'} waiting "
                       f"({describe(waiting)}); the line holds {cap} (queue.max_waiting). Try again when one has started.")
    now = time.time() if now is None else now
    conn.hset(_entry_key(ticket), mapping={
        "type": type_, "title": title, "who": who or "", "database": database or "", "kind": kind or "",
        "queued_at": now, "seq": seq, **{key: "" if value is None else str(value) for key, value in (fields or {}).items()}})
    conn.expire(_entry_key(ticket), KEEP_SECONDS)
    if type_ == RUN:
        conn.set(_alive_key(ticket), "1", ex=ALIVE_SECONDS)
    return ticket


def leave(conn, ticket):
    """Takes the entry out of the line; whether it was in it."""
    present = conn.zscore(_line_key(), ticket) is not None
    _drop(conn, ticket)
    return present


def claim(conn, ticket):
    """Takes the entry off the line so that one process alone starts it; whether this call took it. Its fields stay
    until `leave` (or `restore`), so a start that fails can put it back."""
    return bool(conn.zrem(_line_key(), ticket))


def restore(conn, ticket):
    """Puts a claimed entry back at the place it had."""
    seq = conn.hget(_entry_key(ticket), "seq")
    if seq is not None:
        conn.zadd(_line_key(), {ticket: int(runrecord._s(seq))})


def set_field(conn, ticket, name, value):
    conn.hset(_entry_key(ticket), name, "" if value is None else str(value))


def sign_of_life(conn, ticket):
    conn.set(_alive_key(ticket), "1", ex=ALIVE_SECONDS)


def start_when_free(conn, take, limit_seconds, on_wait=None, kind=runrecord.GENERATION, title=None, who=None,
                    database=None, admitted=False, clock=time.time, sleep=time.sleep):
    """
    Takes the machine with `take()` (it returns the run id, or raises `runrecord.Busy`), waiting in the line for it
    when another run holds it.

    Only the first run in the line tries the machine, so a late comer never jumps the line. It gives up when
    `limit_seconds` have passed (0 or less never waits) and fails naming the run it waited for.

    Args:
        conn: The Redis client.
        take (callable): Starts the run; raises `runrecord.Busy` when the machine is taken.
        limit_seconds (float): The longest to wait; 0 or less does not wait (the first `Busy` is raised).
        on_wait (callable, optional): `on_wait(holder, position)` each time the wait changes (the run that holds the
            machine, and the place among the waiting runs, 1 being next), and `on_wait(None, None)` once it is over.
        kind (str): The kind of run, to name the holder.
        title, who, database (str, optional): What the line shows for this run.
        admitted (bool): The run belongs to the job that is running (its step started it): it gets a place in a full
            line.

    Returns:
        str: What `take` returned.

    Raises:
        LineFull: The line is full.
        runrecord.Busy: The machine stayed taken for `limit_seconds`; the message names the run.
    """
    if limit_seconds <= 0 or not of_type(conn, RUN):
        try:
            return take()
        except runrecord.Busy:
            if limit_seconds <= 0:
                raise
    ticket = join(conn, RUN, title or "A run", who=who, database=database, kind=kind, force=admitted)
    deadline = clock() + limit_seconds
    shown = None
    try:
        while True:
            sign_of_life(conn, ticket)
            waiting = of_type(conn, RUN)
            places = [entry["ticket"] for entry in waiting]
            position = places.index(ticket) + 1 if ticket in places else 1
            holder = runrecord.active_run(conn, kind)
            if position == 1:
                try:
                    return take()
                except runrecord.Busy as busy:
                    holder = runrecord.active_run(conn, kind) or holder
                    last = busy
            else:
                last = runrecord.Busy(runrecord.busy_message(conn, holder))
            if clock() >= deadline:
                raise runrecord.Busy(f"{last}; waited {limit_seconds / 60:g} minutes for it")
            if on_wait is not None and (holder, position) != shown:
                shown = (holder, position)
                on_wait(holder, position)
            sleep(POLL_SECONDS)
    finally:
        leave(conn, ticket)
        if shown is not None and on_wait is not None:
            on_wait(None, None)
