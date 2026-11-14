# planetgen/queue/runline.py

"""
The one waiting list for work that cannot start yet (PERF.154)
==============================================================

Boss, 2026-10-22: "a task that arrives is never lost or failed only because another is running". Everything that has to
wait for the machine waits here, in one ordered list in Redis, up to `queue.max_waiting` entries (8 by default), first
come first served:

* a **run** (`type` `run`): a generation run (or an API job's, a one-sector regenerate or a neighbourhood) that finds
  the machine's `GENERATION` slot taken (`runrecord.start_run`), waiting inside `start_when_free`. The supervisor starts it (OPS.211: `supervisor.RunStarts`), the first of the line when the
  machine is free, and the waiting process follows the run it started. It holds its place
  with a sign of life (`ALIVE_SECONDS`), so one that was killed drops out of the line by itself;
* a **job** (`type` `job`): a Generate page job that finds another job running, waiting for the Generate page's job
  slot (`web/jobs.py`). Its recipe is kept in a private file of the jobs folder (it holds the database environment);
  the line holds the rest. It starts by itself when the job before it succeeded or was canceled and is held for an
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
* `line:<ticket>:alive`  a run's sign of life;
* `line:<ticket>:verdict`  the supervisor's answer to a run that asked to start (`run_id`, `busy` or `failed`);
* `line:starts`     hash: ticket -> the start of a run that does not wait (no place in the line), JSON; a run in the line keeps its own in the entry's `start` field.
"""

import json
import secrets
import time
import traceback

from planetgen.queue import runrecord
from planetgen.util import log
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
            "start": fields.get("start") or None,
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
        sign_of_life(conn, ticket)
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


STARTS_KEY_NAME = "line:starts"
START_WAIT_SECONDS = 60.0
"""float: How long a run that nothing holds up waits for the supervisor to start it before it says the supervisor did not."""


class NotStarted(runrecord.Busy):
    """Nothing holds the machine and the supervisor did not start the run in time."""


def _starts_key():
    return runrecord.PREFIX + STARTS_KEY_NAME


def _verdict_key(ticket):
    return runrecord.PREFIX + f"line:{ticket}:verdict"


def verdict(conn, ticket):
    """What the supervisor decided about the run `ticket` asked for: `{"run_id": ...}`, `{"busy": message}` or `{"failed": message}`; `None` while it has not decided."""
    fields = {runrecord._s(key): runrecord._s(value) for key, value in conn.hgetall(_verdict_key(ticket)).items()}
    return fields or None


def decide(conn, ticket, **fields):
    """The supervisor's answer to the run `ticket` asked for (`verdict`)."""
    conn.hset(_verdict_key(ticket), mapping=fields)
    conn.expire(_verdict_key(ticket), KEEP_VERDICT_SECONDS)


KEEP_VERDICT_SECONDS = 600


def direct_starts(conn):
    """`{ticket: start}` of the runs that asked for the machine without a place in the line (a run that does not wait), `start` the arguments of `runrecord.start_run`."""
    return {runrecord._s(ticket): json.loads(runrecord._s(raw)) for ticket, raw in conn.hgetall(_starts_key()).items()}


def take_direct(conn, ticket):
    """Takes one direct ask so that one process alone acts on it; whether this call took it."""
    return bool(conn.hdel(_starts_key(), ticket))


def waiting_start(conn):
    """The run first in the line and the start it asks for, as `(ticket, start)`; `None` when the line holds no run or its first run starts itself (it has no start to carry out)."""
    waiting = of_type(conn, RUN)
    if waiting and waiting[0]["start"]:
        return waiting[0]["ticket"], json.loads(waiting[0]["start"])
    return None


def start_when_free(conn, start, limit_seconds, ask, on_wait=None, kind=runrecord.GENERATION, title=None, who=None,
                    database=None, admitted=False, clock=time.time, sleep=time.sleep,
                    nested_guard=False):
    """
    Asks the supervisor to start a run (OPS.211) and waits for its answer, in the line for the machine when another run holds it.

    The run does not start here: the ask is left on Redis (a place in the line for a run that waits, a request of its own for one that does not) and the supervisor's `RunStarts` starts it
    when the machine is free and the run is first in the line, so a run is never started by a process that may be gone a moment later. Only the first run in the line is started, so a late
    comer never jumps it. This gives up when `limit_seconds` have passed (0 or less never waits) and fails naming the run it waited for.

    Args:
        start (dict): The arguments of `runrecord.start_run`, as JSON: the stages as dicts and `owner` the process that follows the run.
        limit_seconds (float): The longest to wait for the machine; 0 or less does not wait (the first `Busy` is raised).
        ask (callable): `ask()` makes sure a supervisor is there to answer, called at every poll (the caller's `supervisor.ensure_started`; a test passes one turn of `RunStarts`).
        on_wait (callable, optional): `on_wait(holder, position)` each time the wait changes (the run that holds the
            machine, and the place among the waiting runs, 1 being next), and `on_wait(None, None)` once it is over.
        kind (str): The kind of run, to name the holder.
        title, who, database (str, optional): What the line shows for this run.
        nested_guard (bool): Fail at once when the run holding the machine is followed by this very process (a run started
            inside a run would wait for itself for ever; ADM.136).
        admitted (bool): The run belongs to the job that is running (its step started it): it gets a place in a full
            line.

    Returns:
        str: The id of the run that was started.

    Raises:
        LineFull: The line is full.
        runrecord.Busy: The machine stayed taken for `limit_seconds`; the message names the run.
        NotStarted: Nothing held the machine and no supervisor started the run within `START_WAIT_SECONDS`.
    """
    waits = limit_seconds > 0
    if waits:
        ticket = join(conn, RUN, title or "A run", who=who, database=database, kind=kind,
                      force=admitted or not of_type(conn, RUN), fields={"start": json.dumps(start)})
    else:
        ticket = secrets.token_hex(8)
        conn.hset(_starts_key(), ticket, json.dumps(start))
        conn.expire(_starts_key(), KEEP_VERDICT_SECONDS)
    now = clock()
    deadline = now + limit_seconds
    free_since = now
    shown = None
    try:
        while True:
            if waits:
                sign_of_life(conn, ticket)
            ask()
            found = verdict(conn, ticket)
            if found:
                if found.get("run_id"):
                    return found["run_id"]
                if found.get("failed"):
                    raise RuntimeError(found["failed"])
                raise runrecord.Busy(found.get("busy") or "another run is running")
            holder = runrecord.active_run(conn, kind)
            places = [entry["ticket"] for entry in of_type(conn, RUN)] if waits else []
            position = places.index(ticket) + 1 if ticket in places else 1
            if holder is None:
                if clock() - free_since >= START_WAIT_SECONDS:
                    raise NotStarted("nothing holds the machine, but no supervisor started the run in time; the supervisor may not be running (planetgen.cli.supervisor)")
            else:
                free_since = clock()
                if not waits and clock() - now >= START_WAIT_SECONDS:
                    raise NotStarted("no supervisor answered the ask for a run in time; the supervisor may not be running (planetgen.cli.supervisor)")
            if nested_guard and holder and runrecord.followed_here(conn, holder):
                log.error("A run was started inside a run of the same process; where from:\n"
                                                     + "".join(traceback.format_stack(limit=25)))
                raise runrecord.Busy(f"{runrecord.busy_message(conn, holder)}; this run is started by the same process "
                                     f"that follows it, so it would wait for itself for ever (a run started inside a run)")
            if waits and clock() >= deadline and holder is not None:
                raise runrecord.Busy(f"{runrecord.busy_message(conn, holder)}; waited {limit_seconds / 60:g} minutes for it")
            if on_wait is not None and waits and holder is not None and (holder, position) != shown:
                shown = (holder, position)
                on_wait(holder, position)
            sleep(POLL_SECONDS)
    finally:
        leave(conn, ticket)
        conn.hdel(_starts_key(), ticket)
        if shown is not None and on_wait is not None:
            on_wait(None, None)
