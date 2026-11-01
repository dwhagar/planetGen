# planetgen/queue/runrecord.py

"""
A run's record in Redis (PERF.96, step 2 of the work queue rebuild, docs/design/operations/redis-and-the-work-queue.md 5.1
and 5.2): the one place a run lives. A library with no pages and no workers of its own; PERF.97 to 99 use it.

Keys (`PREFIX` is `planetgen`; the tests give it a prefix of their own):

* `<PREFIX>:run:<id>`            hash: title, who, database, argv, seed, state, the stage list (JSON), the current
                                  stage, times, and what a run needs to be finished later (PERF.182): the galaxy's seed
                                  and the build the run was started on, and the run it finishes, if it does. No expiry
                                  while live; `KEEP_SECONDS` after the run ends.
* `<PREFIX>:run:<id>:s<k>`       hash: stage k's counters (units, units done, predicted cost and cost done, seconds
                                  the units took, first and last finish) and the decayed sums the ETA reads.
* `<PREFIX>:run:<id>:s<k>:done`  set: the indexes of the units that finished; also the barrier.
* `<PREFIX>:run:<id>:s<k>:advanced`  string: set once, by the unit that completes the stage.
* `<PREFIX>:run:<id>:s<k>:results`   hash: unit -> the JSON result it returned (PERF.97), written before the unit
                                  joins the done set, so a unit found done has its result.
* `<PREFIX>:run:<id>:s<k>:progress`  hash: unit -> `"done,estimate"`, a unit's work so far while it runs (PERF.97).
* `<PREFIX>:run:<id>:s<k>:plan`      string: the JSON summary the stage's command gave when it was expanded into
                                  units (PERF.118): what any process following the run needs to read the units'
                                  results. Written before the stage's units are queued.
* `<PREFIX>:run:<id>:banner`     string: the whole hours the site banner shows for the run (UX.3), so every visitor
                                  sees one figure; expires a day after it is written.
* `<PREFIX>:run:<id>:alive`      string: set with a time to live of `ALIVE_SECONDS` by the process that follows the run
                                  (`touch`, at every poll), so a run whose process was killed can be told from a live one.
* `<PREFIX>:run:<id>:events`     stream: the log and timeline, `XADD` with `MAXLEN ~ EVENTS_MAXLEN`.
* `<PREFIX>:job:<job id>:runs`   list: the ids of the runs a Generate page job started, in order (PERF.99), so the
                                  page that shows the job finds its runs; expires with them.
* `<PREFIX>:run:active`          string: the run that owns the machine, `SET NX`, so a second run gets "busy".
* `<PREFIX>:run:active:<kind>`   string: the same for a run of `MAINTENANCE` kind (the cache warm-up, PERF.127) or
                                  `CHECK` kind (the database check, PERF.126): it does not own the machine, so it
                                  neither blocks a generation run nor shows the galaxy-changing banner, and only one
                                  of each kind runs at a time.

A unit that finishes calls `unit_done`, one script run by Redis: it adds the index to the done set (a unit that
runs twice counts once), updates the counters and the decayed sums, and, for the unit that makes the set
complete, sets the `advanced` flag with `NX`, so exactly one unit moves the run on, even if a retried unit runs
twice. The decayed sums are `progress_rate.DecayingRate`'s (a ratio of two sums, each decayed with the time
constant `max(60 s, 20 x mean unit seconds / workers)`), so the live rate here and the terminal bar agree.

Every function takes `now` (Unix seconds, default the wall clock) so tests drive a fake clock. The ETA is a pure
function of what `read_run` returns: `planetgen.queue.runeta`.
"""

import json
import logging
import os
import re
import secrets
import socket
import time
from dataclasses import asdict, dataclass

from planetgen.util import log as debug_log
from planetgen.util import logrecord

PREFIX = "planetgen"
"""str: Every key starts with this."""

KEEP_SECONDS = 30 * 86400
"""int: How long a finished run's keys are kept."""

EVENTS_MAXLEN = 100_000
"""int: The event stream's approximate length limit."""

GENERATION = "generation"
MAINTENANCE = "maintenance"
CHECK = "check"
KINDS = (GENERATION, MAINTENANCE, CHECK)
"""A run changes the galaxy (`GENERATION`, owns the machine), only keeps the site in order (`MAINTENANCE`, PERF.127:
building caches changes nothing a visitor reads) or only reads it (`CHECK`, PERF.126: the database check). Each kind has
one live run at a time, in a place of its own, so a check is never refused for a warm-up and neither holds up a
generation run."""

LIVE_STATES = ("running", "cancelling")
FINAL_STATES = ("done", "failed", "cancelled")

ALIVE_SECONDS = 120
"""int: How long a run's `alive` key lives after the last `touch`: the run's process is taken as gone when it has
run out (or, on the same host, when its process id is no longer running)."""

STALE_FLOOR_SECONDS = 120.0
"""float: A live run with no unit finished and no event for this long (and 3 times its longest unit) is stale."""


class Busy(RuntimeError):
    """Another run owns the machine (`active_run`)."""


@dataclass
class Stage:
    """One stage of a run, as declared when it starts.

    Attributes:
        key (str): The stage's name (`generation/stages.py`).
        label (str): What the pages call it.
        units (int): How many units it queues (1 for a one-shot stage).
        cost (float, optional): The predicted cost of all its units together, in the units of `prior_rate`
            (default `units`: every unit costs 1).
        expected_seconds (float, optional): How long earlier runs took this stage (PERF.32's recorded timing),
            the ETA's fallback while the stage reports nothing.
        prior_rate (float, optional): Cost units per second recorded for this kind of work (`R0`).
        worker_rate (float, optional): Cost units per worker second this server recorded for the stage's last runs
            (`stagetimes`), which `runeta` turns into a rate once the units are known (PERF.184).
        workers (int): Workers the stage runs on, for the ETA's time constant.
        row (dict, optional): The `generation_stage_runs` row the unit that completes the stage writes (PERF.106:
            `command`, `n`, `total`, `label` and `settings`; `generation.stages.hand_over_open_stage`); none for a
            stage that stores no row.
    """

    key: str
    label: str
    units: int = 1
    cost: float = None
    expected_seconds: float = None
    prior_rate: float = None
    worker_rate: float = None
    workers: int = 1
    row: dict = None


RUN_ID_RE = re.compile(r"^\d{8}-\d{6}-[0-9a-f]{8}$")
"""re.Pattern: A run's id (`new_run_id`), checked before an id from a URL is used."""


def _k(*parts):
    return ":".join([PREFIX, *map(str, parts)])


def run_key(run_id):
    return _k("run", run_id)


def stage_key(run_id, index):
    return _k("run", run_id, f"s{index}")


def done_key(run_id, index):
    return _k("run", run_id, f"s{index}", "done")


def advanced_key(run_id, index):
    return _k("run", run_id, f"s{index}", "advanced")


def events_key(run_id):
    return _k("run", run_id, "events")


def results_key(run_id, index):
    return _k("run", run_id, f"s{index}", "results")


def progress_key(run_id, index):
    return _k("run", run_id, f"s{index}", "progress")


def plan_key(run_id, index):
    return _k("run", run_id, f"s{index}", "plan")


def debug_key(run_id):
    """The stream of debug lines the run's workers relay to the process that started it (`planetgen.util.log.relayed`)."""
    return _k("run", run_id, "debug")


def banner_key(run_id):
    return _k("run", run_id, "banner")


def job_runs_key(job_id):
    """The list of the runs the Generate page job `job_id` started (PERF.99)."""
    return _k("job", job_id, "runs")


def alive_key(run_id):
    return _k("run", run_id, "alive")


def active_key(kind=GENERATION):
    """The key that names the live run of `kind`."""
    return _k("run", "active") if kind == GENERATION else _k("run", "active", kind)


def new_run_id(now=None):
    """A run's id: the UTC date and time and eight hex digits, so ids sort by start."""
    return time.strftime("%Y%m%d-%H%M%S", time.gmtime(time.time() if now is None else now)) + "-" + secrets.token_hex(4)


def _now(now):
    return time.time() if now is None else float(now)


def _s(value):
    return value.decode("utf-8") if isinstance(value, bytes) else value


_UNIT_DONE = """
local added = redis.call('SADD', KEYS[1], ARGV[1])
if added == 0 then return {0, 0, redis.call('SCARD', KEYS[1])} end
local cost, secs, now, workers = tonumber(ARGV[2]), tonumber(ARGV[3]), tonumber(ARGV[4]), tonumber(ARGV[5])
local d = redis.call('HMGET', KEYS[2], 'units', 'done', 'cost_done', 'seconds', 'last', 'u', 't', 'share',
                     'pending', 'longest_gap', 'longest_unit', 'peak_mb')
local units = tonumber(d[1])
local done = tonumber(d[2] or 0) + 1
local seconds = tonumber(d[4] or 0) + secs
local last = tonumber(d[5] or now)
local u, t = tonumber(d[6] or 0), tonumber(d[7] or 0)
local share, pending, gap = tonumber(d[8] or 0), tonumber(d[9] or 0), tonumber(d[10] or 0)
local longest = math.max(tonumber(d[11] or 0), secs)
local peak = math.max(tonumber(d[12] or 0), tonumber(ARGV[6] or 0))
local tau = math.max(60, 20 * (seconds / done) / math.max(workers, 1))
if now > last then
  local interval = now - last
  if interval > gap then gap = interval end
  local decay = math.exp(-interval / tau)
  local weight = tau * (1 - decay)
  share = weight / interval
  t = t * decay + weight
  u = u * decay + (pending + cost) * share
  pending = 0
  last = now
elseif share == 0 then
  pending = pending + cost
else
  u = u + cost * share
end
redis.call('HSET', KEYS[2], 'done', done, 'cost_done', tonumber(d[3] or 0) + cost, 'seconds', seconds,
           'last', last, 'last_end', now, 'u', u, 't', t, 'share', share, 'pending', pending, 'longest_gap', gap,
           'longest_unit', longest, 'peak_mb', peak)
redis.call('HSETNX', KEYS[2], 'first_end', now)
redis.call('HSET', KEYS[4], 'updated', now)
local advance = 0
if redis.call('SCARD', KEYS[1]) >= units then
  if redis.call('SET', KEYS[3], ARGV[4], 'NX') then
    advance = 1
    redis.call('HSET', KEYS[2], 'ended', now)
  end
end
return {1, advance, redis.call('SCARD', KEYS[1])}
"""

_RELEASE = """
if redis.call('GET', KEYS[1]) == ARGV[1] then return redis.call('DEL', KEYS[1]) end
return 0
"""


@dataclass
class UnitDone:
    """What `unit_done` found: `counted` (false for a unit that had already finished), `advance` (this unit
    completed the stage and is the one to start the next), `done` (units finished so far)."""

    counted: bool
    advance: bool
    done: int


def active_run(conn, kind=GENERATION):
    """The id of the live run of `kind` (by default the one that owns the machine), or `None`."""
    value = conn.get(active_key(kind))
    return _s(value) if value else None


def start_run(conn, title, stages, run_id=None, who=None, database=None, argv=None, seed=None, now=None,
              steps=None, kind=GENERATION, disregard_missing_data=False, debug="", job=None, galaxy_seed=None,
              build=None, finishes=None):
    """
    Records a new run and takes the machine for it (a `MAINTENANCE` or `CHECK` run takes only its own place, so a
    generation run and one of those may be live together).

    Args:
        kind (str): One of `KINDS`.
        stages (list[Stage]): The run's stages, in order (at least one).
        steps (list[dict], optional): The wrapper's step list (`queue.commands.Step.as_dict`), one per stage; the
            unit that completes a stage reads the next one from here (PERF.118).
        disregard_missing_data (bool): Whether the run was told to take a labelled fallback for a missing downloaded data
            file (GEN.599); every unit's worker is set the same (`queue.units.run_unit`).
        debug (str): `""` for an ordinary run; for a debug run, the areas it shows (`planetgen.util.log.parse_areas`,
            `"all"` for every one). Every unit's worker then relays the debug lines it writes to the stream
            `debug_key(run_id)`, and the process that started the run prints them (`planetgen.util.log.RelayReader`).
        job (str, optional): The Generate page job that started the run (PERF.99): kept in the record and appended to
            `job_runs_key(job)`.
        galaxy_seed (str, optional): The galaxy's seed as hex and `build` the version key of the code that started the
            run (PERF.182): what a later `Finish run` checks are still the same before it carries on.
        finishes (str, optional): The id of the run this one carries on from (`generation/finish_run.py`).

    Returns:
        str: The run id.

    Raises:
        Busy: Another run of that kind is running.
    """
    if not stages:
        raise ValueError("a run has at least one stage")
    if kind not in KINDS:
        raise ValueError(f"a run is one of {KINDS}")
    now = _now(now)
    run_id = run_id or new_run_id(now)
    touch(conn, run_id)         # before the machine is taken: a run that holds it always has a sign of life
    if not conn.set(active_key(kind), run_id, nx=True):
        if not (_free_dead_holder(conn, kind) and conn.set(active_key(kind), run_id, nx=True)):
            conn.delete(alive_key(run_id))
            raise Busy(busy_message(conn, active_run(conn, kind)))
    declared = [asdict(stage) for stage in stages]
    for stage in declared:
        stage["cost"] = float(stage["units"] if stage["cost"] is None else stage["cost"])
    pipe = conn.pipeline()
    pipe.hset(run_key(run_id), mapping={
        "title": title, "kind": kind, "who": who or "", "database": database or "", "argv": json.dumps(argv or []),
        "seed": "" if seed is None else str(seed), "state": "running", "stages": json.dumps(declared),
        "steps": json.dumps(steps or []), "disregard_missing_data": "1" if disregard_missing_data else "",
        "debug": debug or "", "current": 0, "created": now, "started": now, "updated": now,
        "finished": "", "error": "", "owner": _owner(), "job": job or "", "galaxy_seed": galaxy_seed or "",
        "build": build or "", "finishes": finishes or ""})
    if job:
        pipe.rpush(job_runs_key(job), run_id)
        pipe.expire(job_runs_key(job), KEEP_SECONDS)
    for index, stage in enumerate(declared):
        pipe.hset(stage_key(run_id, index), mapping={"units": stage["units"], "cost_total": stage["cost"], "done": 0,
                                                     "cost_done": 0, "seconds": 0})
    pipe.execute()
    log(conn, run_id, f"run started: {title}", kind="run.start", now=now)
    return run_id


def busy_message(conn, run_id):
    """
    What a refusal says about the run that holds the machine: `run <id> is running` and, from its record, what it is
    (title, kind, when it started) and who follows it, so the real blocker is named instead of an id alone.
    """
    if not run_id:
        return "another run is running"
    record = {_s(k): _s(v) for k, v in (conn.hgetall(run_key(run_id)) or {}).items()}
    if not record:
        return f"run {run_id} is running"
    started = record.get("started") or ""
    try:
        started = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime(float(started)))
    except ValueError:
        started = "at an unknown time"
    return (f"run {run_id} is running ({record.get('title') or 'untitled'}, {record.get('kind') or GENERATION} kind, "
            f"started {started}, followed by {record.get('owner') or 'an unknown process'})")


def _owner():
    """Who follows a run started here: `host:process id`."""
    return f"{socket.gethostname()}:{os.getpid()}"


def touch(conn, run_id):
    """The process that follows the run says it is alive (`ALIVE_SECONDS`); called at every poll."""
    conn.set(alive_key(run_id), "1", ex=ALIVE_SECONDS)


def _pid_alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except OSError:
        return True         # of another user: running
    return True


def holder_is_dead(conn, run_id):
    """
    Whether the process that follows the run `run_id` is gone, so the run can never end by itself: it started on
    this host and its process id is not running, or (another host, or a run recorded before runs said who follows
    them) its `alive` key has run out. A run whose record is gone or already over is dead too: nothing follows it.
    """
    owner = _s(conn.hget(run_key(run_id), "owner")) or ""
    state = _s(conn.hget(run_key(run_id), "state")) or ""
    if state and state not in LIVE_STATES:
        return True
    host, _, pid = owner.rpartition(":")
    if state and host == socket.gethostname() and pid.isdigit():
        return not _pid_alive(int(pid))
    return not conn.exists(alive_key(run_id))


def _free_dead_holder(conn, kind):
    """
    If the run that holds the machine for `kind` is dead (`holder_is_dead`), ends it as failed (its waiting units
    then do nothing) and frees the machine. A run killed with `kill -9`, by the out-of-memory killer or by a restart
    never ended itself, and every new run was refused as busy for ever.

    Returns:
        bool: Whether the machine was held by a dead run and is free now.
    """
    holder = active_run(conn, kind)
    if holder is None or not holder_is_dead(conn, holder):
        return False
    finish_run(conn, holder, "failed", error="its process was gone (killed or restarted) before the run ended")
    conn.register_script(_RELEASE)(keys=[active_key(kind)], args=[holder])
    return True


def read_steps(conn, run_id):
    """The step list `start_run` stored (`[]` for a run started without one)."""
    return json.loads(_s(conn.hget(run_key(run_id), "steps")) or "[]")


def set_stage_plan(conn, run_id, index, units, cost, plan, now=None):
    """
    Records what stage `index` turned out to be when its command was expanded (PERF.118): `units` units of total
    predicted `cost`, and `plan`, the JSON summary a process following the run reads its results by. Called once,
    before the units are queued, so no unit can finish against the old count.
    """
    now = _now(now)
    declared = json.loads(_s(conn.hget(run_key(run_id), "stages")))
    declared[index]["units"] = int(units)
    declared[index]["cost"] = float(cost)
    pipe = conn.pipeline()
    pipe.hset(run_key(run_id), mapping={"stages": json.dumps(declared), "updated": now})
    pipe.hset(stage_key(run_id, index), mapping={"units": int(units), "cost_total": float(cost)})
    pipe.set(plan_key(run_id, index), json.dumps(plan))
    pipe.execute()


def read_plan(conn, run_id, index):
    """Stage `index`'s plan (`set_stage_plan`), or `None` while the stage has not been expanded."""
    raw = conn.get(plan_key(run_id, index))
    return None if raw is None else json.loads(_s(raw))


def set_conclusion(conn, run_id, index, conclusion, now=None):
    """Stores what the end of stage `index` found (`CommandSpec.conclude`) and marks the stage concluded, for the
    process following the run (`conclusion`)."""
    conn.hset(stage_key(run_id, index), mapping={"conclusion": json.dumps(conclusion), "concluded": _now(now)})


def conclusion(conn, run_id, index):
    """What `set_conclusion` stored for stage `index`, or `None` while the stage is not concluded."""
    raw = conn.hget(stage_key(run_id, index), "conclusion")
    return None if raw is None else json.loads(_s(raw))


def stored_conclusion(conn, run_id, key):
    """What the end of the stage `key` of run `run_id` found (`conclusion`), `None` when the run has no such stage or the stage
    was not concluded (a run carries on from a run's finished stages with theirs, PERF.182)."""
    declared = json.loads(_s(conn.hget(run_key(run_id), "stages")) or "[]")
    for index, stage in enumerate(declared):
        if _s(stage["key"]) == key:
            return conclusion(conn, run_id, index)
    return None


def begin_stage(conn, run_id, index, now=None):
    """Marks stage `index` as started now (the clock its ETA counts from)."""
    now = _now(now)
    pipe = conn.pipeline()
    pipe.hset(stage_key(run_id, index), mapping={"started": now, "last": now})
    pipe.hset(run_key(run_id), mapping={"current": index, "updated": now})
    pipe.execute()
    declared = json.loads(conn.hget(run_key(run_id), "stages"))
    log(conn, run_id, logrecord.stage_start(index + 1, len(declared), _s(declared[index]["label"])),
        kind=logrecord.STAGE_START, now=now, stage=index, stages=len(declared))


def set_opened_clock(conn, run_id, index, clock):
    """Records the database's clock (as text) when stage `index` opened, for a stage a later step's parameter refers to
    (`wrapper.opened_by`, PERF.213)."""
    conn.hset(stage_key(run_id, index), "opened_clock", clock)


def opened_clock(conn, run_id, index):
    """The database's clock when stage `index` opened (`set_opened_clock`), `None` for a stage that recorded none."""
    value = conn.hget(stage_key(run_id, index), "opened_clock")
    return _s(value) if value else None


def end_stage(conn, run_id, index, now=None):
    """Marks stage `index` as ended now: the end of a stage that had no unit to end it."""
    conn.hset(stage_key(run_id, index), mapping={"ended": _now(now)})


def unit_done(conn, run_id, index, unit, cost=None, seconds=0.0, now=None, peak_mb=0.0):
    """
    Records that unit `unit` of stage `index` finished (after `seconds` of work, at cost `cost`: default the
    stage's mean unit cost, and `peak_mb` the most resident memory the process that did it took, PERF.226: the stage keeps the
    largest). One atomic script; see the module docstring.

    Returns:
        UnitDone: `advance` is true for the one unit that completed the stage.
    """
    now = _now(now)
    stage = json.loads(conn.hget(run_key(run_id), "stages"))[index]
    if cost is None:
        cost = stage["cost"] / max(stage["units"], 1)
    script = conn.register_script(_UNIT_DONE)
    counted, advance, done = script(
        keys=[done_key(run_id, index), stage_key(run_id, index), advanced_key(run_id, index), run_key(run_id)],
        args=[unit, cost, seconds, now, max(int(stage.get("workers") or 1), 1), peak_mb])
    if counted:
        log(conn, run_id, f"unit {unit} done", kind="unit.done", now=now, stage=index, unit=unit)
    if advance:
        started = conn.hget(stage_key(run_id, index), "started")
        seconds_taken = max(now - float(_s(started)), 0.0) if started else 0.0
        declared = json.loads(conn.hget(run_key(run_id), "stages"))
        peak = float(_s(conn.hget(stage_key(run_id, index), "peak_mb")) or 0.0)
        log(conn, run_id, logrecord.stage_finish(index + 1, len(declared), seconds_taken, peak_mb=peak),
            kind=logrecord.STAGE_FINISH, now=now, stage=index, stages=len(declared), seconds=round(seconds_taken, 1),
            **({"peak_mb": round(peak)} if peak else {}))
        _record_stage_time(conn, run_id, index, stage["key"])
    return UnitDone(bool(counted), bool(advance), int(done))


def _record_stage_time(conn, run_id, index, key):
    """Folds the stage that has just ended into what the server knows of its kind (`stagetimes`, PERF.184). Best
    effort: an estimate's memory is not worth failing a unit for."""
    from planetgen.queue import stagetimes
    try:
        counters = {_s(k): _s(v) for k, v in conn.hgetall(stage_key(run_id, index)).items()}
        started, ended = counters.get("started"), counters.get("ended")
        if started is None or ended is None:
            return
        stagetimes.record(conn, key, float(ended) - float(started), float(counters.get("cost_done") or 0.0),
                          float(counters.get("seconds") or 0.0))
    except Exception:  # noqa: BLE001
        pass


def request_cancel(conn, run_id, now=None):
    """Marks a live run as cancelling (the unit that completes a stage then queues no next stage). Returns
    whether it was live."""
    now = _now(now)
    if _s(conn.hget(run_key(run_id), "state")) not in LIVE_STATES:
        return False
    conn.hset(run_key(run_id), mapping={"state": "cancelling", "updated": now})
    log(conn, run_id, "cancel requested", kind="run.cancelling", now=now)
    return True


def is_cancelling(conn, run_id):
    return _s(conn.hget(run_key(run_id), "state")) == "cancelling"


def finish_run(conn, run_id, state="done", error=None, now=None):
    """
    Ends the run in `state` (`done`, `failed` or `cancelled`): records it, frees the machine if this run owns
    it, and gives the run's keys `KEEP_SECONDS` to live. A run already over is left as it is.

    Returns:
        bool: Whether this call ended it.
    """
    if state not in FINAL_STATES:
        raise ValueError(f"a run ends as one of {FINAL_STATES}")
    now = _now(now)
    if _s(conn.hget(run_key(run_id), "state")) not in LIVE_STATES:
        return False
    conn.hset(run_key(run_id), mapping={"state": state, "finished": now, "updated": now, "error": error or ""})
    kind = _s(conn.hget(run_key(run_id), "kind")) or GENERATION
    conn.register_script(_RELEASE)(keys=[active_key(kind)], args=[run_id])
    log(conn, run_id, error or f"run {state}", kind=f"run.{state}", level="error" if state == "failed" else "info",
        now=now)
    pipe = conn.pipeline()
    for key in _keys_of(conn, run_id):
        pipe.expire(key, KEEP_SECONDS)
    pipe.execute()
    return True


def _keys_of(conn, run_id, every=False):
    """The keys the run keeps for `KEEP_SECONDS` once it is over (the stage count comes from the record), and with
    `every` the ones that expire by themselves too (the banner figure, the sign of life, the debug stream)."""
    stages = len(json.loads(conn.hget(run_key(run_id), "stages") or "[]"))
    keys = [run_key(run_id), events_key(run_id)]
    if every:
        keys += [banner_key(run_id), alive_key(run_id), debug_key(run_id)]
    for index in range(stages):
        keys += [stage_key(run_id, index), done_key(run_id, index), advanced_key(run_id, index),
                 results_key(run_id, index), progress_key(run_id, index), plan_key(run_id, index)]
    return keys


def delete_run(conn, run_id):
    """
    Removes a finished run's record, counters and events (the Queue page's delete). A run still live is left alone.

    Returns:
        bool: Whether it was deleted.
    """
    state = _s(conn.hget(run_key(run_id), "state"))
    if state is None or state in LIVE_STATES:
        return False
    conn.delete(*_keys_of(conn, run_id, every=True))
    return True


def log(conn, run_id, message, level="info", kind="log", now=None, **fields):
    """Appends a line to the run's event stream; returns its id."""
    entry = {"t": _now(now), "level": level, "kind": kind, "message": message,
             **{name: value for name, value in fields.items() if value is not None}}
    if debug_log.debug_log_active():
        debug_log.trace(f"run {run_id}: {message} ({kind}" + "".join(f", {name}={value}" for name, value in fields.items()
                                                                  if value is not None) + ")", area="queue")
    return _s(conn.xadd(events_key(run_id), entry, maxlen=EVENTS_MAXLEN, approximate=True))


def read_events(conn, run_id, after="0", count=500):
    """The events after stream id `after` (`"0"` for all), oldest first: `[(id, {field: value})]` with text
    values. Follow a run by passing the last id seen."""
    rows = conn.xrange(events_key(run_id), min="(" + after if after != "0" else "-", max="+", count=count)
    return [(_s(entry_id), {_s(k): _s(v) for k, v in fields.items()}) for entry_id, fields in rows]


EVENT_FIELD_NAMES = ("stage", "stages", "unit", "seconds")
"""tuple: The fields of an event, beside its time, level, kind and message, that `event_record` carries into the record."""


def event_record(fields):
    """
    The logging layer's record (`logrecord.Record`) of an event `read_events` gave (UX.122): its time, its level as the
    severity, the component `run`, its message, and the fields `event` (its kind), `stage` (the stage's number from 1),
    `stages` (how many the run has), `unit` and `seconds` where it has them.
    """
    level = fields.get("level", "info")
    named = {"event": fields.get("kind", "log")}
    if "stage" in fields:
        named["stage"] = int(fields["stage"]) + 1
    for name in ("stages", "unit"):
        if name in fields:
            named[name] = int(fields[name])
    if "seconds" in fields:
        named["seconds"] = float(fields["seconds"])
    return logrecord.Record.at(float(fields.get("t", 0)), level if level in logrecord.SYSLOG_CODES else "info", "run",
                               fields.get("message", ""), **named)


def log_records(conn, run_id):
    """Every event of the run, oldest first, as records (`event_record`)."""
    records, after = [], "0"
    while True:
        rows = read_events(conn, run_id, after)
        if not rows:
            return records
        for entry_id, fields in rows:
            records.append(event_record(fields))
            after = entry_id


def log_text(conn, run_id):
    """The run's log as text: one file line (`logrecord.file_line`) to an event. What `write_log_file` writes and the
    admin run log shows, so the file and the page are the same lines."""
    return "".join(logrecord.file_line(record) + "\n" for record in log_records(conn, run_id))


def write_log_file(conn, run_id, directory, keep=20):
    """
    Writes the run's events to `<directory>/<run id>.log` (the live log, kept at the run's end), then removes
    the oldest `.log` files beyond the newest `keep`. Returns the path.
    """
    import os
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, f"{run_id}.log")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(log_text(conn, run_id))
    logs = sorted(name for name in os.listdir(directory) if name.endswith(".log"))
    for name in logs[:-keep] if keep else logs:
        os.remove(os.path.join(directory, name))
    return path


_STAGE_NUMBERS = ("units", "done", "cost_total", "cost_done", "seconds", "started", "ended", "last", "last_end",
                  "first_end", "u", "t", "longest_gap", "longest_unit", "peak_mb")


def read_run(conn, run_id):
    """
    The run as one dict for pages and the ETA (`runeta.estimate`): the record (`id`, `title`, `kind`, `who`, `database`,
    `job`, `argv`, `state`, `current`, `created`, `started`, `updated`, `finished`, `error`) and `stages`, each the declared
    stage with its counters merged in (numbers, `None` where not set yet), and `galaxy_seed`, `build` and `finishes`
    (PERF.182: `None` for a run that recorded none). `None` for an unknown run.
    """
    pipe = conn.pipeline()
    pipe.hgetall(run_key(run_id))
    record = {_s(k): _s(v) for k, v in pipe.execute()[0].items()}
    if not record:
        return None
    declared = json.loads(record["stages"])
    pipe = conn.pipeline()
    for index in range(len(declared)):
        pipe.hgetall(stage_key(run_id, index))
    run = {"id": run_id, "title": record["title"], "kind": record.get("kind") or GENERATION, "who": record["who"],
           "database": record["database"], "job": record.get("job") or None, "argv": json.loads(record.get("argv") or "[]"),
           "galaxy_seed": record.get("galaxy_seed") or None, "build": record.get("build") or None,
           "finishes": record.get("finishes") or None,
           "state": record["state"], "error": record["error"] or None, "current": int(float(record["current"])),
           "stages": []}
    for name in ("created", "started", "updated", "finished"):
        run[name] = float(record[name]) if record.get(name) else None
    for stage, counters in zip(declared, pipe.execute()):
        counters = {_s(k): _s(v) for k, v in counters.items()}
        stage = dict(stage)
        for name in _STAGE_NUMBERS:
            stage[name] = float(counters[name]) if name in counters else None
        for name in ("units", "done"):
            stage[name] = int(stage[name] or 0)
        stage["cost_total"] = stage["cost_total"] if stage["cost_total"] is not None else stage["cost"]
        stage["cost_done"] = stage["cost_done"] or 0.0
        run["stages"].append(stage)
    return run


def is_stale(run, now=None):
    """Whether a live run has shown no sign of life -- no unit finished, no event -- for `max(120 s, 3 x its
    longest unit so far)`: the Queue page offers to clear it."""
    if run is None or run["state"] not in LIVE_STATES:
        return False
    now = _now(now)
    longest = max([stage.get("longest_gap") or 0.0 for stage in run["stages"]] or [0.0])
    return now - (run["updated"] or run["started"] or now) > max(STALE_FLOOR_SECONDS, 3.0 * longest)


class EventHandler(logging.Handler):
    """
    Writes what the process that follows a run logs at the normal level or above to the run's event stream, so the
    stream is the run's whole live log (PERF.99): the lines the console shows and the stage and unit events beside
    them. A Redis that does not answer costs the line, never the run.
    """

    def __init__(self, conn, run_id):
        super().__init__(level=logging.INFO)
        self.conn, self.run_id = conn, run_id

    def emit(self, record):
        try:
            log(self.conn, self.run_id, record.getMessage(), level=record.levelname.lower())
        except Exception:  # noqa: BLE001 -- the line is lost, the run goes on
            pass


def job_runs(conn, job_id):
    """The ids of the runs the Generate page job `job_id` started, oldest first."""
    return [_s(run_id) for run_id in conn.lrange(job_runs_key(job_id), 0, -1)]


def recent_runs(conn, limit=20):
    """Ids of the runs whose records exist, newest first (the ids sort by start time)."""
    prefix = _k("run") + ":"
    found = set()
    for key in conn.scan_iter(match=prefix + "*", count=500):
        rest = _s(key)[len(prefix):]
        if rest and rest[:2].isdigit() and ":" not in rest:
            found.add(rest)
    return sorted(found, reverse=True)[:limit]


__all__ = ["CHECK", "GENERATION", "KINDS", "MAINTENANCE", "Busy", "Stage", "UnitDone", "active_run", "begin_stage", "finish_run", "is_cancelling", "is_stale",
           "EventHandler", "job_runs", "log", "new_run_id", "read_events", "read_plan", "read_run", "read_steps", "recent_runs", "request_cancel",
           "set_stage_plan", "start_run", "stored_conclusion", "unit_done", "write_log_file", "log_text", "log_records", "event_record"]
