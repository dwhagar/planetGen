# planetgen/queue/units.py

"""
A stage of a run on the new work queue (PERF.97, step 3 of the rebuild, docs/design/operations/redis-and-the-work-queue.md
5.1): the units of one stage queued on the bulk queue, run by the idle-timeout workers (`workers.ensure_workers`),
and followed by the process that started the run until the stage is done. The wrapper (`queue/wrapper.py`, PERF.118)
cuts a stage into its units and starts a run; this file queues, runs and follows them.

* **A unit is data.** It is an RQ job for `run_unit` with the run id, the stage, the unit's index, a `kind` and a
  JSON dict (`planetgen.queue.redisqueue.JsonSerializer`: the queue carries data, not code). The `kind` names
  a command in `commands.REGISTRY`, a table in code, so what Redis holds can pick a command but can't name a
  function. A unit carries no password: the worker reads the database's credentials from its own settings.
* **The barrier is the done set.** A unit that finishes stores its result (a JSON value), then adds its index to
  the stage's done set (`runrecord.unit_done`); the unit that completes the set is the one that moves the run on.
  A unit that runs twice (RQ's retry after a killed worker) counts once. The result is written before the index,
  so a unit found done has its result.
* **Nothing depends on the worker count.** The units are cut from the work (`static_chunks`), each draws from a
  stream of its own, and the order the workers take them in changes nothing.
* **The starter follows** (`follow_stage`): it reads the done set, the stored results and the units' progress
  every poll, hands them to the caller, and replaces a dead worker as it goes. It is a client: if it dies, the run
  goes on as far as its workers take it, and the Queue page (PERF.99) can show and clear it.
* **A failed unit ends the run** (`unit_failed`, RQ's `on_failure`, which RQ calls on every failure, before it decides
  on the retry); the units still waiting see the run is no longer live and do nothing. A unit that is already in
  the done set did its work: a failure after that is RQ's own bookkeeping, and neither ends the run nor runs the
  unit again.

The unit that completes the last stage finishes the run; the unit that completes an earlier one has the wrapper queue
the next (`wrapper.advance`).
"""

import json
import os
import resource
import time

import rq
from rq.job import Callback

from planetgen.datasets import fetch
from planetgen.db import store
from planetgen.queue import commands, redisqueue, runrecord
from planetgen.util import log

UNIT_FUNCTION = "planetgen.queue.units.run_unit"
FAILURE_FUNCTION = "planetgen.queue.units.unit_failed"

MAX_UNITS = 4000
"""int: The most units one stage queues; Redis memory and queuing cost follow the unit count (about 7 MB at this
many), so a stage with more work uses bigger units."""

STATIC_UNITS = 512
"""int: About how many units a stage is cut into when the work allows. The count comes from the work and never
from the number of workers, so the same work is the same units at any worker count; 512 keeps the largest
sensible pool (the 2 million-object ceiling's 95% of 64 cores) busy to the end."""

FAILED_KEEP_SECONDS = 86400
"""int: How long RQ keeps a failed unit's job."""

POLL_SECONDS = 0.25
"""float: How often `follow_stage` reads the run."""

ENSURE_SECONDS = 5.0
"""float: How often `follow_stage` has the workers topped up, so a dead one is replaced."""

REPORT_SECONDS = 0.25
"""float: How often a unit writes its progress."""

ENQUEUE_BATCH = 500
"""int: How many units go to Redis in one pipeline."""

class RunStopped(RuntimeError):
    """The run ended (failed, was cancelled or went stale) before its stage was done."""


def job_id(run_id, stage, unit):
    """The RQ job id of a unit, so a unit can be found from the run (cancel, PERF.99)."""
    return f"{run_id}-s{stage}-u{unit}"


def unit_timeout(expected_seconds=None):
    """
    How long a unit of a stage may run before RQ ends it (PERF.106): `redisqueue.job_timeout` of the stage's
    expected seconds (`runrecord.Stage.expected_seconds`); the ceiling (6 hours) when the stage has no recorded
    time. Never -1 (no limit).
    """
    return redisqueue.job_timeout(expected_seconds)


def coarsen(count, size):
    """
    The size of a piece, `size` or a whole multiple of it, that cuts `count` items into at most `MAX_UNITS` pieces: a
    stage whose cut has a fixed piece size (a chunk of tiles, a range of ids) uses this so a big galaxy makes bigger
    units rather than more than the cap, which `enqueue_units` refuses. The same `count` and `size` always give the same
    piece, whatever the number of workers.
    """
    pieces = -(-count // size) if count > 0 else 0
    return size * max(1, -(-pieces // MAX_UNITS))


def static_chunks(rows, units=STATIC_UNITS, max_rows=None):
    """
    Cuts `rows` (consecutive items, each `row[-1]` its weight) into about `units` units of about equal total
    weight, in order: the cut depends on the rows alone, so it is the same whatever number of workers takes
    them. A unit holds at most `max_rows` rows (default: any number), and never more than `MAX_UNITS` units are
    made: a cut that would make more is made coarser (about half the cap's worth by weight, and by rows).

    Returns:
        list[list]: The rows of each unit, in the order given.
    """
    rows = list(rows)
    if not rows:
        return []
    half = max(MAX_UNITS // 2 - 1, 1)
    target = sum(row[-1] for row in rows) / max(min(units, half), 1)
    limit = max(max_rows or len(rows), -(-len(rows) // half))
    chunks, current, weight = [], [], 0.0
    for row in rows:
        current.append(row)
        weight += row[-1]
        if weight >= target or len(current) >= limit:
            chunks.append(current)
            current, weight = [], 0.0
    if current:
        chunks.append(current)
    return chunks


def enqueue_units(conn, run_id, stage, kind, units, queue_name=None, expected_seconds=None):
    """
    Queues the units of stage `stage` of the run.

    Args:
        units (list[tuple]): `(args, cost)` per unit, in unit order (the index is the position): a JSON-able dict
            and the unit's predicted cost for the stage's bar and ETA.
        expected_seconds (float, optional): How long the stage is expected to take, for each unit's time limit
            (`unit_timeout`).

    Raises:
        ValueError: More than `MAX_UNITS` units, or a kind that is not a command with units.
    """
    if kind not in commands.REGISTRY or commands.spec(kind).inline:
        raise ValueError(f"no unit kind {kind!r}")
    if len(units) > MAX_UNITS:
        raise ValueError(f"a stage queues at most {MAX_UNITS} units, not {len(units)}")
    queue = redisqueue.queue(queue_name or redisqueue.BULK_QUEUE, conn)
    timeout = unit_timeout(expected_seconds)
    for start in range(0, len(units), ENQUEUE_BATCH):
        queue.enqueue_many([
            queue.prepare_data(
                UNIT_FUNCTION, args=(run_id, stage, index, kind, args, cost), job_id=job_id(run_id, stage, index),
                timeout=timeout, result_ttl=0, failure_ttl=FAILED_KEEP_SECONDS,
                retry=rq.Retry(max=1), on_failure=Callback(FAILURE_FUNCTION), description=f"{kind} {index}")
            for index, (args, cost) in enumerate(units[start:start + ENQUEUE_BATCH], start)])


class _Reporter:
    """What a unit calls with its work so far: at most every `REPORT_SECONDS` it writes `done,estimate` to the
    stage's progress hash for the follower."""

    def __init__(self, conn, run_id, stage, unit):
        self.conn, self.key, self.unit = conn, runrecord.progress_key(run_id, stage), unit
        self.last = 0.0

    def __call__(self, done, estimate):
        now = time.monotonic()
        if now - self.last >= REPORT_SECONDS:
            self.last = now
            self.conn.hset(self.key, self.unit, f"{float(done)},{float(estimate)}")


def run_unit(run_id, stage, unit, kind, args, cost=None):
    """
    The RQ job of one unit: does it, stores its result, joins the done set, and, when it completes its stage, has the
    wrapper queue the next one or, after the last, ends the run. A unit that finds the run no longer live (cancelled,
    failed) does nothing.

    Returns:
        dict: `seconds` the unit took, or `skipped`.
    """
    conn = rq.get_current_job().connection
    state, disregard, debug = conn.hmget(runrecord.run_key(run_id), "state", "disregard_missing_data", "debug")
    if runrecord._s(state) != "running" or conn.sismember(runrecord.done_key(run_id, stage), unit):
        return {"skipped": True}
    fetch.set_disregard_missing_data(bool(runrecord._s(disregard)))      # GEN.599: the run's own switch, in every worker
    wanted = runrecord._s(debug)
    if not wanted:
        return _do_unit(conn, run_id, stage, unit, kind, args, cost)
    # A debug run (`--debug`, or debug mode in config.json): this worker's stdout goes nowhere, so what the unit
    # logs is relayed to the process that started the run, which prints it as it arrives.
    with log.relayed(conn, runrecord.debug_key(run_id), "" if wanted == "all" else wanted):
        return _do_unit(conn, run_id, stage, unit, kind, args, cost)


def peak_resident_mb():
    """The most resident memory this process has held, MB (PERF.226): a unit's horse is a fork of its worker, so it counts the worker's
    pages too, which is what one worker of the stage takes of the host."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0     # kilobytes on Linux


def _do_unit(conn, run_id, stage, unit, kind, args, cost):
    started = time.monotonic()
    log.trace(f"unit {unit} of stage {stage} of run {run_id} starts in worker {os.getpid()}: {kind}, cost {cost}, "
              f"arguments {_short(args)}", area="queue")
    result = commands.unit_function(kind)(args, _Reporter(conn, run_id, stage, unit))
    seconds = time.monotonic() - started
    log.trace(f"unit {unit} of stage {stage} of run {run_id} finished its work in {seconds:.3f} s", area="queue")
    log.flush_relayed()         # before the unit counts as done, so the starter has every line when it sees it
    conn.hset(runrecord.results_key(run_id, stage), unit, json.dumps(result))
    outcome = runrecord.unit_done(conn, run_id, stage, unit, cost=cost, seconds=seconds, peak_mb=peak_resident_mb())
    conn.hdel(runrecord.progress_key(run_id, stage), unit)
    log.trace(f"unit {unit} of stage {stage} of run {run_id} recorded; "
              f"{'it completes the stage: the next is queued' if outcome.advance else f'{outcome.done} unit(s) done'}",
              area="queue")
    if outcome.advance:
        from planetgen.queue import wrapper         # the wrapper queues units, so it imports this file
        _end_stage(conn, run_id, stage, wrapper)
    return {"seconds": seconds}


def _end_stage(conn, run_id, stage, wrapper):
    """
    What the unit that completes `stage` does for the run (PERF.172): the command's conclusion, the stage's
    statistics, and the opening of the next stage. A failure in the conclusion or in opening the next stage (its
    `prepare`, or reading the galaxy to cut it) is the run's failure, said at once with the stage named; left to
    RQ it would read as a unit that failed after it was done, and the run would wait for a stage nobody queues.
    """
    try:
        took = wrapper.conclude_stage(conn, run_id, stage)
    except Exception as exc:
        runrecord.finish_run(conn, run_id, "failed", error=f"the end of stage {wrapper.stage_label(conn, run_id, stage)} "
                                                           f"failed: {type(exc).__name__}: {exc}")
        raise
    wrapper.finish_stage(conn, run_id, stage, took)
    try:
        wrapper.advance(conn, run_id, stage + 1)
    except Exception as exc:
        runrecord.finish_run(conn, run_id, "failed", error=f"stage {wrapper.stage_label(conn, run_id, stage + 1)} could not "
                                                           f"start: {type(exc).__name__}: {exc}")
        raise


def _short(value, limit=300):
    text = json.dumps(value, default=str)
    return text if len(text) <= limit else text[:limit] + f"... ({len(text)} chars)"


def warm_unit(job):
    """
    PERF.146: runs the command's `warm` function, if it has one, in the worker before RQ forks the horse for `job`.
    A horse is a new process, so everything a unit works out from the galaxy alone (its shape, which holds the component
    field built from tens of thousands of stored clusters, and the sampling tables of the scatters) was worked out
    again in every unit: seconds of arithmetic before the first star, a unit after a unit. Done here, in the
    worker, the caches it fills (`store._FIELD_BY_SIGNATURE`, `star_population._bright_table`) are in the worker's
    memory and each horse starts with them; a table that changed meanwhile is read again, by its signature. Best
    effort: a unit that cannot warm meets the same error itself.
    """
    if job.func_name != UNIT_FUNCTION:
        return
    kind, args = job.args[3], job.args[4]
    spec = commands.spec(kind)
    if spec.warm is None:
        return
    try:
        commands.warm_function(kind)(args)
    except Exception:
        pass


def warm_shape(unit_args):
    """A `warm` function for the units that read the galaxy's shape first thing (`unit_args["database"]`): reads it, so
    the worker holds its component field (`warm_unit`). The connection is closed and its pool dropped before the fork,
    so no socket is shared with a horse."""
    named = unit_args["database"]
    config = store.MySQLConfig(host=named["host"], port=named["port"], user=named["user"], database=named["database"])
    try:
        conn = store.get_connection(config)
        try:
            store.get_galaxy_shape(conn)
        finally:
            conn.close()
    finally:
        store.close_pool(config)


def unit_failed(job, connection, exc_type, exc_value, _traceback):
    """RQ's `on_failure` for a unit: the run fails, with the unit and the error recorded, unless the unit's result is
    already stored and counted (RQ failed it after it was done, say `ValueError: Execution ... not found in Redis`
    when it cleaned up): its sector or shells are saved, so the run goes on."""
    run_id, stage, unit = job.args[:3]
    if connection.sismember(runrecord.done_key(run_id, stage), unit):
        log.debug(f"Unit {unit} of stage {stage} of run {run_id} was done when RQ failed it ({exc_type.__name__}: "
                  f"{exc_value}); the run goes on.")
        return
    from planetgen.queue import wrapper             # the wrapper queues units, so it imports this file
    runrecord.finish_run(connection, run_id, "failed",
                         error=f"unit {unit} of stage {wrapper.stage_name(connection, run_id, stage)} failed: "
                               f"{exc_type.__name__}: {exc_value}")


def horse_killed(job, retpid, ret_val, _rusage):
    """
    RQ's `work_horse_killed_handler` (PERF.106): the process that was running a unit was killed (out of memory, or
    `kill -9`). The run fails now, with the unit named, instead of when the worker's heartbeat runs out; the retry RQ
    makes finds the run no longer live and does nothing. A job that is not a unit of a run is left to RQ.
    """
    if job.func_name != UNIT_FUNCTION:
        return
    run_id, stage, unit = job.args[:3]
    from planetgen.queue import wrapper             # the wrapper queues units, so it imports this file
    runrecord.finish_run(job.connection, run_id, "failed",
                         error=f"unit {unit} of stage {wrapper.stage_name(job.connection, run_id, stage)} was killed: "
                               f"its worker process {retpid} ended with "
                               f"status {ret_val} (the machine may be out of memory)")


def _progress(conn, run_id, stage):
    """`{unit: (done, estimate)}` of the units reporting now."""
    found = {}
    for unit, raw in conn.hgetall(runrecord.progress_key(run_id, stage)).items():
        try:
            done, estimate = (float(part) for part in runrecord._s(raw).split(","))
            found[int(unit)] = (done, estimate)
        except ValueError:
            continue
    return found


def follow_stage(conn, run_id, stage, units, on_result, on_progress=None, ensure=None, clock=time.monotonic,
                 sleep=time.sleep, on_debug=None):
    """
    Follows stage `stage` of the run until its `units` units are done: calls `on_result(unit, result)` once for
    each as it finishes and `on_progress(unit, done, estimate)` for those still working, and `ensure()` every
    `ENSURE_SECONDS` (starting workers, replacing dead ones), and `on_debug()` at every poll and when the stage ends
    (a debug run prints the lines its workers relayed, `planetgen.util.log.RelayReader`).

    Raises:
        RunStopped: The run is no longer live and the stage is not complete (a unit failed, or it was cancelled).
    """
    seen = set()
    last_ensure = float("-inf")

    def drain():
        done = {int(unit) for unit in conn.smembers(runrecord.done_key(run_id, stage))} - seen
        if done:
            ordered = sorted(done)
            for unit, raw in zip(ordered, conn.hmget(runrecord.results_key(run_id, stage), ordered)):
                on_result(unit, json.loads(runrecord._s(raw)))
            seen.update(ordered)

    while True:
        runrecord.touch(conn, run_id)
        drain()
        if on_debug is not None:
            on_debug()
        if len(seen) >= units:
            return
        state = runrecord._s(conn.hget(runrecord.run_key(run_id), "state"))
        if state not in runrecord.LIVE_STATES:
            drain()         # the last units may have finished between the two reads
            if on_debug is not None:
                on_debug()
            if len(seen) >= units:
                return
            error = runrecord._s(conn.hget(runrecord.run_key(run_id), "error"))
            raise RunStopped(error or f"the run {state} before the stage was done")
        if on_progress is not None:
            for unit, (done, estimate) in _progress(conn, run_id, stage).items():
                if unit not in seen:
                    on_progress(unit, done, estimate)
        if ensure is not None and clock() - last_ensure >= ENSURE_SECONDS:
            last_ensure = clock()
            ensure()
        sleep(POLL_SECONDS)
