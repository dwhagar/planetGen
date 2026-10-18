# planetgen/queue/runrecord.py

"""
A run's record in Redis (PERF.96, step 2 of the work queue rebuild, docs/design/redis-and-the-work-queue.md 5.1
and 5.2): the one place a run lives. A library with no pages and no workers of its own; PERF.97 to 99 use it.

Keys (`PREFIX` is `planetgen`; the tests give it a prefix of their own):

* `<PREFIX>:run:<id>`            hash: title, who, database, argv, seed, state, the stage list (JSON), the current
                                  stage, times. No expiry while live; `KEEP_SECONDS` after the run ends.
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
* `<PREFIX>:run:<id>:events`     stream: the log and timeline, `XADD` with `MAXLEN ~ EVENTS_MAXLEN`.
* `<PREFIX>:run:active`          string: the run that owns the machine, `SET NX`, so a second run gets "busy".

A unit that finishes calls `unit_done`, one script run by Redis: it adds the index to the done set (a unit that
runs twice counts once), updates the counters and the decayed sums, and, for the unit that makes the set
complete, sets the `advanced` flag with `NX`, so exactly one unit moves the run on, even if a retried unit runs
twice. The decayed sums are `progress_rate.DecayingRate`'s (a ratio of two sums, each decayed with the time
constant `max(60 s, 20 x mean unit seconds / workers)`), so the live rate here and the terminal bar agree.

Every function takes `now` (Unix seconds, default the wall clock) so tests drive a fake clock. The ETA is a pure
function of what `read_run` returns: `planetgen.queue.runeta`.
"""

import json
import secrets
import time
from dataclasses import asdict, dataclass

PREFIX = "planetgen"
"""str: Every key starts with this."""

KEEP_SECONDS = 30 * 86400
"""int: How long a finished run's keys are kept."""

EVENTS_MAXLEN = 100_000
"""int: The event stream's approximate length limit."""

LIVE_STATES = ("running", "cancelling")
FINAL_STATES = ("done", "failed", "cancelled")

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
        workers (int): Workers the stage runs on, for the ETA's time constant.
    """

    key: str
    label: str
    units: int = 1
    cost: float = None
    expected_seconds: float = None
    prior_rate: float = None
    workers: int = 1


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


def banner_key(run_id):
    return _k("run", run_id, "banner")


def active_key():
    return _k("run", "active")


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
                     'pending', 'longest_gap')
local units = tonumber(d[1])
local done = tonumber(d[2] or 0) + 1
local seconds = tonumber(d[4] or 0) + secs
local last = tonumber(d[5] or now)
local u, t = tonumber(d[6] or 0), tonumber(d[7] or 0)
local share, pending, gap = tonumber(d[8] or 0), tonumber(d[9] or 0), tonumber(d[10] or 0)
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
           'last', last, 'last_end', now, 'u', u, 't', t, 'share', share, 'pending', pending, 'longest_gap', gap)
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


def active_run(conn):
    """The id of the run that owns the machine, or `None`."""
    value = conn.get(active_key())
    return _s(value) if value else None


def start_run(conn, title, stages, run_id=None, who=None, database=None, argv=None, seed=None, now=None,
              steps=None):
    """
    Records a new run and takes the machine for it.

    Args:
        stages (list[Stage]): The run's stages, in order (at least one).
        steps (list[dict], optional): The wrapper's step list (`queue.commands.Step.as_dict`), one per stage; the
            unit that completes a stage reads the next one from here (PERF.118).

    Returns:
        str: The run id.

    Raises:
        Busy: Another run owns the machine.
    """
    if not stages:
        raise ValueError("a run has at least one stage")
    now = _now(now)
    run_id = run_id or new_run_id(now)
    if not conn.set(active_key(), run_id, nx=True):
        raise Busy(f"run {active_run(conn)} is running")
    declared = [asdict(stage) for stage in stages]
    for stage in declared:
        stage["cost"] = float(stage["units"] if stage["cost"] is None else stage["cost"])
    pipe = conn.pipeline()
    pipe.hset(run_key(run_id), mapping={
        "title": title, "who": who or "", "database": database or "", "argv": json.dumps(argv or []),
        "seed": "" if seed is None else str(seed), "state": "running", "stages": json.dumps(declared),
        "steps": json.dumps(steps or []), "current": 0, "created": now, "started": now, "updated": now,
        "finished": "", "error": ""})
    for index, stage in enumerate(declared):
        pipe.hset(stage_key(run_id, index), mapping={"units": stage["units"], "cost_total": stage["cost"], "done": 0,
                                                     "cost_done": 0, "seconds": 0})
    pipe.execute()
    log(conn, run_id, f"run started: {title}", kind="run.start", now=now)
    return run_id


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


def begin_stage(conn, run_id, index, now=None):
    """Marks stage `index` as started now (the clock its ETA counts from)."""
    now = _now(now)
    pipe = conn.pipeline()
    pipe.hset(stage_key(run_id, index), mapping={"started": now, "last": now})
    pipe.hset(run_key(run_id), mapping={"current": index, "updated": now})
    pipe.execute()
    label = _s(json.loads(conn.hget(run_key(run_id), "stages"))[index]["label"])
    log(conn, run_id, f"stage started: {label}", kind="stage.start", now=now, stage=index)


def unit_done(conn, run_id, index, unit, cost=None, seconds=0.0, now=None):
    """
    Records that unit `unit` of stage `index` finished (after `seconds` of work, at cost `cost`: default the
    stage's mean unit cost). One atomic script; see the module docstring.

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
        args=[unit, cost, seconds, now, max(int(stage.get("workers") or 1), 1)])
    if counted:
        log(conn, run_id, f"unit {unit} done", kind="unit.done", now=now, stage=index, unit=unit)
    if advance:
        log(conn, run_id, f"stage done: {stage['label']}", kind="stage.done", now=now, stage=index)
    return UnitDone(bool(counted), bool(advance), int(done))


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
    conn.register_script(_RELEASE)(keys=[active_key()], args=[run_id])
    log(conn, run_id, error or f"run {state}", kind=f"run.{state}", level="error" if state == "failed" else "info",
        now=now)
    stages = len(json.loads(conn.hget(run_key(run_id), "stages")))
    keys = [run_key(run_id), events_key(run_id)]
    for index in range(stages):
        keys += [stage_key(run_id, index), done_key(run_id, index), advanced_key(run_id, index),
                 results_key(run_id, index), progress_key(run_id, index), plan_key(run_id, index)]
    pipe = conn.pipeline()
    for key in keys:
        pipe.expire(key, KEEP_SECONDS)
    pipe.execute()
    return True


def log(conn, run_id, message, level="info", kind="log", now=None, **fields):
    """Appends a line to the run's event stream; returns its id."""
    entry = {"t": _now(now), "level": level, "kind": kind, "message": message,
             **{name: value for name, value in fields.items() if value is not None}}
    return _s(conn.xadd(events_key(run_id), entry, maxlen=EVENTS_MAXLEN, approximate=True))


def read_events(conn, run_id, after="0", count=500):
    """The events after stream id `after` (`"0"` for all), oldest first: `[(id, {field: value})]` with text
    values. Follow a run by passing the last id seen."""
    rows = conn.xrange(events_key(run_id), min="(" + after if after != "0" else "-", max="+", count=count)
    return [(_s(entry_id), {_s(k): _s(v) for k, v in fields.items()}) for entry_id, fields in rows]


def write_log_file(conn, run_id, directory, keep=20):
    """
    Writes the run's events to `<directory>/<run id>.log` (the live log, kept at the run's end), then removes
    the oldest `.log` files beyond the newest `keep`. Returns the path.
    """
    import os
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, f"{run_id}.log")
    after = "0"
    with open(path, "w", encoding="utf-8") as handle:
        while True:
            rows = read_events(conn, run_id, after)
            if not rows:
                break
            for entry_id, fields in rows:
                stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(float(fields.get("t", 0))))
                handle.write(f"{stamp} {fields.get('level', 'info').upper():<5} {fields.get('message', '')}\n")
                after = entry_id
    logs = sorted(name for name in os.listdir(directory) if name.endswith(".log"))
    for name in logs[:-keep] if keep else logs:
        os.remove(os.path.join(directory, name))
    return path


_STAGE_NUMBERS = ("units", "done", "cost_total", "cost_done", "seconds", "started", "ended", "last", "last_end",
                  "first_end", "u", "t", "longest_gap")


def read_run(conn, run_id):
    """
    The run as one dict for pages and the ETA (`runeta.estimate`): the record (`id`, `title`, `who`, `database`,
    `state`, `current`, `created`, `started`, `updated`, `finished`, `error`) and `stages`, each the declared
    stage with its counters merged in (numbers, `None` where not set yet). `None` for an unknown run.
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
    run = {"id": run_id, "title": record["title"], "who": record["who"], "database": record["database"],
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


def recent_runs(conn, limit=20):
    """Ids of the runs whose records exist, newest first (the ids sort by start time)."""
    prefix = _k("run") + ":"
    found = set()
    for key in conn.scan_iter(match=prefix + "*", count=500):
        rest = _s(key)[len(prefix):]
        if rest and rest[:2].isdigit() and ":" not in rest:
            found.add(rest)
    return sorted(found, reverse=True)[:limit]


__all__ = ["Busy", "Stage", "UnitDone", "active_run", "begin_stage", "finish_run", "is_cancelling", "is_stale",
           "log", "new_run_id", "read_events", "read_plan", "read_run", "read_steps", "recent_runs", "request_cancel",
           "set_stage_plan", "start_run", "unit_done", "write_log_file"]
