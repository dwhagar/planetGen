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

The unit that completes the last stage finishes the run; the unit that completes an earlier one queues the opening of the next
(`open_stage`, PERF.254): the cut of a stage is a job of its own, with its own claim on the memory gate (`claims.CUT_MB`,
`planetgen.queue.memorygate`), so it runs when the host can give it the bulk lane to itself, in a horse that is named when it fails or is
killed, and not inside the horse of the last unit of the stage before.
"""

import json
import os
import resource
import signal
import time

import rq
from rq.job import Callback

from planetgen.datasets import fetch
from planetgen.db import store
from planetgen.queue import attemptlog, attempts, commands, commandgaps, engine, failures, held, jobprofiles, planrows, redisqueue, requires, runrecord, unitbeat, unitrates, usage
from planetgen.util import log
from planetgen.util.format import format_duration_seconds, format_megabytes

UNIT_FUNCTION = "planetgen.queue.units.run_unit"
FAILURE_FUNCTION = "planetgen.queue.units.unit_failed"
OPEN_FUNCTION = "planetgen.queue.units.open_stage"
OPEN_FAILURE_FUNCTION = "planetgen.queue.units.open_failed"

MAX_UNITS = 4000
"""int: The most units one stage queues; Redis memory and queuing cost follow the unit count (about 7 MB at this
many), so a stage with more work uses bigger units."""

STATIC_UNITS = 512
"""int: About how many units a stage is cut into when the work allows. The count comes from the work and never
from the number of workers, so the same work is the same units at any worker count; 512 keeps the largest
sensible pool (the 2 million-object ceiling's 95% of 64 cores) busy to the end."""

FAILED_KEEP_SECONDS = 86400
"""int: How long RQ keeps a failed unit's job."""


ENDING_LEASE_SECONDS = 1800
"""int: How long the process that ends a stage holds the lease (`runrecord.ending_key`) against a second attempt; a kill with a retry left frees it (`horse_killed`)."""

OPEN_TRIES = 4
"""int: How many times the queuing of an opening is tried before the run fails (a Redis hiccup is not the run's end)."""

POLL_SECONDS = 0.25
"""float: How often `follow_stage` reads the run."""

ENSURE_SECONDS = 5.0
"""float: How often `follow_stage` has the workers topped up, so a dead one is replaced."""

QUIET_LINE_SECONDS = 1800.0
"""float: How long a stage may finish no unit before `follow_stage` says so in the log, and how often it says so after (OPS.253): a stage whose units are long, or wait on one another, was silent for hours."""

ENQUEUE_BATCH = 500
"""int: How many units go to Redis in one pipeline."""

class RunStopped(RuntimeError):
    """The run ended (failed, was canceled or went stale) before its stage was done."""


KILLED_EXCEPTIONS = failures.KILLED_EXCEPTIONS
"""tuple[str]: The names of the exceptions with which RQ ends a unit from outside (its time limit, its worker lost): a run a unit ended with one of them was killed (OPS.144), not failed by its own code."""


def job_id(run_id, stage, unit):
    """The RQ job id of a unit, so a unit can be found from the run (cancel, PERF.99)."""
    return f"{run_id}-s{stage}-u{unit}"


def unit_timeout(expected_seconds=None):
    """
    The time limit RQ puts on a unit: the ceiling (`redisqueue.MAX_TIMEOUT_SECONDS`, 6 hours), whatever the stage's expected time (OPS.140). It is a backstop only: a unit
    that is dead or stuck is found by its beat (`unitwatch`, `unitstuck`) in a minute or a few, and a unit whose counter moves is never ended for taking longer than
    estimated (Boss, 11-09). `expected_seconds` is kept in the signature for the callers that name it. Never -1 (no limit).
    """
    return redisqueue.MAX_TIMEOUT_SECONDS


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


EXPECTED_SECONDS_ARG = "expected_seconds"
"""str: The key of a unit's arguments that gives the seconds that unit alone is predicted to take, for the heartbeat judge (its expected rate, OPS.135); the time limit does not read it (`unit_timeout`)."""


def enqueue_units(conn, run_id, stage, kind, units, queue_name=None, expected_seconds=None, rows=None, skip=frozenset()):
    """
    Queues the units of stage `stage` of the run.

    Args:
        units (list[tuple]): `(args, cost)` per unit, in unit order (the index is the position): a JSON-able dict
            and the unit's predicted cost for the stage's bar and ETA. The dict may name `expected_seconds`, the unit's own
            predicted time (`EXPECTED_SECONDS_ARG`).
        expected_seconds (float, optional): How long the stage is expected to take, for each unit's time limit
            (`unit_timeout`).
        rows (list[int], optional): The id of each unit's row in the job ledger (`planrows.record_cut`, OPS.216), in unit order, handed to the unit (`run_unit`'s `row`).
        skip (frozenset[int]): The numbers of units that are not queued: the ones the run this one carries on from did (`resumed.carry`, OPS.222). The others keep their numbers.

    Raises:
        ValueError: More than `MAX_UNITS` units, or a kind that is not a command with units.
    """
    if kind not in commands.REGISTRY or commands.spec(kind).inline:
        raise ValueError(f"no unit kind {kind!r}")
    if len(units) > MAX_UNITS:
        raise ValueError(f"a stage queues at most {MAX_UNITS} units, not {len(units)}")
    if rows is not None and len(rows) != len(units):
        raise ValueError(f"the plan has {len(rows)} rows for {len(units)} units")
    queue = redisqueue.queue(queue_name or redisqueue.BULK_QUEUE, conn)
    timeout = unit_timeout(expected_seconds)
    for start in range(0, len(units), ENQUEUE_BATCH):
        batch = [_prepared(queue, run_id, stage, index, kind, args, cost, None if rows is None else rows[index], timeout)
                 for index, (args, cost) in enumerate(units[start:start + ENQUEUE_BATCH], start) if index not in skip]
        if batch:
            queue.enqueue_many(batch)


def _prepared(queue, run_id, stage, index, kind, args, cost, row, timeout):
    """The RQ job of one unit, ready to queue."""
    return queue.prepare_data(
        UNIT_FUNCTION, args=(run_id, stage, index, kind, args, cost, row), job_id=job_id(run_id, stage, index),
        timeout=timeout, result_ttl=0, failure_ttl=FAILED_KEEP_SECONDS,
        retry=rq.Retry(max=1), on_failure=Callback(FAILURE_FUNCTION), description=f"{kind} {index}")


def requeue_unit(conn, run_id, stage, index, kind, args, cost, row, queue_name=None, expected_seconds=None):
    """
    Queues one unit of the run again with the index and the row it was cut with (OPS.221: a held unit after "Retry held"). What RQ kept of the unit's failed job (its record and its place in the failed
    registry, which share the job id) is removed first.
    """
    queue = redisqueue.queue(queue_name or redisqueue.BULK_QUEUE, conn)
    try:
        rq.job.Job.fetch(job_id(run_id, stage, index), connection=conn).delete()
    except rq.exceptions.NoSuchJobError:
        pass
    queue.enqueue_many([_prepared(queue, run_id, stage, index, kind, args, cost, row, unit_timeout(expected_seconds))])


def row_of(job):
    """The ledger row id (OPS.216) a unit's RQ job was queued with, or `None` for a job that has none (an opening, a run cut before the plan was rows)."""
    args = getattr(job, "args", None) or ()
    return args[6] if getattr(job, "func_name", None) == UNIT_FUNCTION and len(args) > 6 else None


def run_unit(run_id, stage, unit, kind, args, cost=None, row=None):
    """
    The RQ job of one unit: does it, stores its result, joins the done set, and, when it completes its stage, has the
    wrapper queue the next one or, after the last, ends the run. A unit that finds the run no longer live (canceled,
    failed) does nothing.

    `row` is the unit's row in the job ledger (OPS.216), which the fence of OPS.217 will hold.

    Returns:
        dict: `seconds` the unit took, or `skipped`.
    """
    return run_unit_on(rq.get_current_job().connection, run_id, stage, unit, kind, args, cost, row)


def run_unit_on(conn, run_id, stage, unit, kind, args, cost=None, row=None, halt=None):
    """
    A unit run on the Redis connection `conn`, whoever calls it: `run_unit` for the RQ job of a unit, `forkedjob.run` for the job process the template forked (OPS.226).

    `halt`, when given, is a callable the engine asks between two pieces: it returns the reason the unit is to stop there (a request from the supervisor, the supervisor gone) or `None`.
    """
    live, wanted = _run_settings(conn, run_id)
    if not live:
        return {"skipped": True}
    if conn.sismember(runrecord.done_key(run_id, stage), unit):
        if stage_unended(conn, run_id, stage):          # OPS.186: the attempt that finished the stage was killed before it queued the next opening
            repair_unended_stages(conn, run_id)
            return {"resumed": True}
        return {"skipped": True}
    if held.is_held(conn, run_id, stage, unit):         # OPS.221: RQ's retry of a unit that failed for good; "Retry held" queues it again by hand
        return {"skipped": True}
    result = _relayed(conn, run_id, wanted, lambda: _do_unit(conn, run_id, stage, unit, kind, args, cost, row, halt))
    held.settle(conn, run_id)                           # a unit that ends beside a held one may have been the last thing that could run
    return result


def _run_settings(conn, run_id):
    """`(live, debug areas)` of the run for a job of it: whether it is running, and, having set the worker's GEN.599 switch to the run's
    (`fetch.set_disregard_missing_data`), the debug areas its workers relay (`""` for an ordinary run)."""
    state, disregard, debug = conn.hmget(runrecord.run_key(run_id), "state", "disregard_missing_data", "debug")
    if runrecord._s(state) != "running":
        return False, ""
    fetch.set_disregard_missing_data(bool(runrecord._s(disregard)))      # GEN.599: the run's own switch, in every worker
    return True, runrecord._s(debug) or ""


def _relayed(conn, run_id, wanted, work):
    """`work()`, with what it logs relayed to the process that started the run for a debug run (`wanted`, the areas: this worker's stdout
    goes nowhere, so the starter prints the lines as they arrive)."""
    if not wanted:
        return work()
    with log.relayed(conn, runrecord.debug_key(run_id), "" if wanted == "all" else wanted):
        return work()


def _proportional_mb(pid):
    """The memory a process holds, MB, with the pages it shares split among the processes that share them (the `Pss` of
    /proc/<pid>/smaps_rollup); 0.0 where the system does not say."""
    try:
        with open(f"/proc/{pid}/smaps_rollup", encoding="ascii") as handle:
            for line in handle:
                if line.startswith("Pss:"):
                    return int(line.split()[1]) / 1024.0                   # kilobytes
    except (OSError, ValueError):
        pass
    return 0.0


def peak_resident_mb():
    """
    What the worker that ran this unit holds of the host's memory, MB (PERF.226): the horse's own and its worker's together, shared
    pages counted once (`Pss` of the two), or, where the system does not say, the most resident memory this process has held. A horse is
    a fork of its worker, so its resident size counts the worker's tables but not the worker's own pages beside it: a worker and its
    horse measured 1,487 MB resident in the horse and 2,180 MB together.
    """
    resident = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0     # kilobytes on Linux
    return max(resident, resident_mb())


def resident_mb():
    """What this process and its parent (the worker beside a work horse) hold now, MB, shared pages counted once (`_proportional_mb`); 0.0 where the system does not say."""
    return _proportional_mb(os.getpid()) + _proportional_mb(os.getppid())


def _profile_unit(conn, run_id, stage, unit, kind, finished, worked, peak, figures):
    """
    OPS.155: the unit's row in `job_profiles` (items, worked seconds, memory and processor, peak and average) and its figures in the stage's usage hash, which the stage's own row is
    worked out from when the stage ends (`runrecord._record_stage_time`). Best effort: a profile is not worth a unit.
    """
    if not jobprofiles.enabled():
        return
    try:
        jobprofiles.record(f"unit:{kind}", items=finished, seconds=worked, peak_mb=peak, avg_mb=figures["mem_avg_mb"], cpu_peak=figures["cpu_peak"], cpu_avg=figures["cpu_avg"],
                           seconds_peak=worked, seconds_avg=worked)
        conn.hset(runrecord.usage_key(run_id, stage), unit, json.dumps({"seconds": worked, "wall": figures["seconds"], "peak_mb": peak, "avg_mb": figures["mem_avg_mb"],
                                                                      "cpu_peak": figures["cpu_peak"], "cpu_seconds": figures["cpu_seconds"]}))
        conn.expire(runrecord.usage_key(run_id, stage), runrecord.KEEP_SECONDS)
    except Exception as exc:        # noqa: BLE001
        log.debug(f"Job profiles: unit {unit} of stage {stage} of run {run_id} was not profiled ({exc}).", area="queue")


def _do_unit(conn, run_id, stage, unit, kind, args, cost, row=None, halt=None):
    """
    Takes the attempt in the ledger (`attempts.begin`, OPS.217: the unit's attempt number is raised and written to the galaxy), does the unit (`_work_of_unit`) and writes how the attempt ended.
    An attempt that finds another has taken the unit over (`attempts.Replaced`) has rolled back, counts nothing and ends `replaced`; the other attempt's result stands. A unit run again after a
    passing failure (OPS.220) has an attempt for each run; the one that ends here is the last.
    """
    attempt = attempts.begin(conn, run_id, row, args, pid=os.getpid())
    latest = [attempt]            # `_work_of_unit` puts the attempt of each run again here
    with attemptlog.capture(conn, attempt, run_id, stage, unit, kind):        # OPS.223: a restarted attempt writes a debug log of its own
        try:
            result = _work_of_unit(conn, run_id, stage, unit, kind, args, cost, attempt, halt, latest)
        except attempts.Replaced as replaced:
            log.debug(f"Unit {unit} of stage {stage} of run {run_id} ({kind}) was replaced by a later attempt and stops: {replaced}.", area="queue")
            conn.hdel(runrecord.progress_key(run_id, stage), unit)
            if latest[0] is not None:
                latest[0].end("replaced", peak_mb=peak_resident_mb())
            return {"replaced": True}
        except BaseException as exc:
            if latest[0] is not None:
                kind_of_error, text, outcome = attempts.failure(exc)
                latest[0].end(outcome, error_kind=kind_of_error, error_text=f"[{failures.classify(exc)}] {text}", stack_text=text, peak_mb=peak_resident_mb())
            raise
        if latest[0] is not None:
            latest[0].end("done" if "seconds" in result else "lost", peak_mb=peak_resident_mb())
    return result


SLEEP = time.sleep
"""callable: How the engine waits before it runs a unit again (a test puts its own clock here)."""


def _wait_to_retry(conn, run_id, seconds):
    """Waits `seconds`, a second at a time; returns whether the run is still live at the end (a cancel or a failure of the run in the meantime ends the wait)."""
    left = seconds
    while left > 0:
        SLEEP(min(left, 1.0))
        left -= 1.0
        if not _run_settings(conn, run_id)[0]:
            return False
    return _run_settings(conn, run_id)[0]


def _retry_attempt(conn, run_id, attempt, args, latest, exc):
    """The attempt of the next run of a unit that failed with `exc` (OPS.220): this one ends `failed`, with its kind, and the next is begun (a unit with no row has none to end or begin)."""
    if attempt is None:
        return None
    kind_of_error, text, outcome = attempts.failure(exc)
    attempt.end(outcome, error_kind=kind_of_error, error_text=f"[{failures.classify(exc)}] {text}", stack_text=text, peak_mb=peak_resident_mb())
    latest[0] = attempts.begin(conn, run_id, attempt.row, args, pid=os.getpid())
    return latest[0]


def _work_of_unit(conn, run_id, stage, unit, kind, args, cost, attempt, halt=None, latest=None):
    started = time.monotonic()
    latest = [attempt] if latest is None else latest
    log.trace(f"unit {unit} of stage {stage} of run {run_id} starts in worker {os.getpid()}: {kind}, cost {cost}, "
              f"arguments {_short(args)}", area="queue")
    used = usage.Usage(resident_mb)         # OPS.155: the memory and processor of this unit, peak and average
    with unitbeat.UnitBeat(conn, run_id, stage, unit, resident=resident_mb, stack_file=unitbeat.stack_path(run_id, stage, unit), usage=used) as report:       # OPS.136: a thread of this horse writes the beat while the work loop only sets numbers
        retries = 0
        while True:
            try:
                result = engine.run(kind, args, report, redis=conn, run_id=run_id, stage=stage, unit=unit, job=job_id(run_id, stage, unit), attempt=attempt, halt=halt)
                break
            except engine.UnitStopped as stopped:         # OPS.191: the engine found the run no longer live between two pieces; the unit stays undone, as one a cancel stops does
                log.debug(f"Unit {unit} of stage {stage} of run {run_id} ({kind}) ends between two pieces: {stopped}.", area="queue")
                conn.hdel(runrecord.progress_key(run_id, stage), unit)
                return {"stopped": True}
            except Exception as exc:        # noqa: BLE001 -- what to do with a failure is `failures`' to say (OPS.220)
                wait = failures.retry_wait(conn, run_id, exc, retries)
                if wait is None:
                    raise
                retries += 1
                log.warning(f"Unit {unit} of stage {stage} of run {run_id} ({kind}) failed with a passing error ({type(exc).__name__}: {_short(str(exc))}); "
                            f"it is run again in {format_duration_seconds(wait)} (run again {retries} time(s) in a row).")
                runrecord.log(conn, run_id, f"unit {unit} of stage {stage} run again in {format_duration_seconds(wait)} after {type(exc).__name__}", level="warning", kind="unit.retry",
                              stage=stage)
                attempt = _retry_attempt(conn, run_id, attempt, args, latest, exc)
                with unitbeat.waiting("retry"):
                    live = _wait_to_retry(conn, run_id, wait)
                if not live:
                    conn.hdel(runrecord.progress_key(run_id, stage), unit)
                    return {"stopped": True}
        gap = report.longest_gap()
        finished, worked = report.totals()
    commandgaps.note(conn, run_id, kind, gap)       # OPS.135
    unitrates.note(conn, run_id, kind, finished, worked)        # PERF.147: how fast this command's units work, for the ETA of a unit's bar
    peak = peak_resident_mb()
    _profile_unit(conn, run_id, stage, unit, kind, finished, worked, peak, used.summary())       # OPS.155: what this kind of unit costs, kept for the self-tuning jobs (nothing reads it yet)
    seconds = time.monotonic() - started
    log.trace(f"unit {unit} of stage {stage} of run {run_id} finished its work in {format_duration_seconds(seconds)}", area="queue")
    log.flush_relayed()         # before the unit counts as done, so the starter has every line when it sees it
    if attempt is not None:
        attempt.finish()        # OPS.217: the done mark is true in the galaxy first, and only where this attempt is still the unit's (`attempts.Replaced` otherwise)
    conn.hset(runrecord.results_key(run_id, stage), unit, json.dumps(result))
    outcome = runrecord.unit_done(conn, run_id, stage, unit, cost=cost, seconds=seconds, peak_mb=peak)
    conn.hdel(runrecord.progress_key(run_id, stage), unit)
    log.trace(f"unit {unit} of stage {stage} of run {run_id} recorded; "
              f"{'it completes the stage: the next is queued' if outcome.advance else f'{outcome.done} unit(s) done'}",
              area="queue")
    if outcome.advance:
        from planetgen.queue import wrapper         # the wrapper queues units, so it imports this file
        _end_stage(conn, run_id, stage, wrapper)
    return {"seconds": seconds}


def stage_unended(conn, run_id, stage):
    """OPS.186: whether `stage` was completed (its `advanced` flag is set) and its end has not queued the next opening: the state a kill or a failure in the end leaves."""
    return bool(conn.exists(runrecord.advanced_key(run_id, stage))) and planrows.RELEASE not in planrows.ended(conn, run_id, stage)


def reopen_pending(conn, run_id, queue_name, queued=()):
    """
    OPS.222: queues the opening of the first stage that was never opened when the run was held (the opening job found the run not live and skipped itself, so nothing else would queue it), once the stage before
    it has ended and unless the queue still holds the opening (`queued`: the job ids waiting or running). Returns the index of the stage whose opening it queued, or `None`.
    """
    count = len(runrecord.declared_stages(conn, run_id))
    for index in range(count):
        if conn.exists(runrecord.opened_key(run_id, index)):
            continue
        if (index == 0 or runrecord._s(conn.hget(runrecord.stage_key(run_id, index - 1), "ended"))) and open_job_id(run_id, index) not in queued:
            enqueue_open(conn, run_id, index, queue_name)
            return index
        return None
    return None


def unended_stages(conn, run_id):
    """OPS.186: the stages of a live run that were completed and whose end has not queued the next opening, in order (what the reconciler of OPS.210 looks for)."""
    if runrecord._s(conn.hget(runrecord.run_key(run_id), "state")) != "running":
        return []
    count = len(json.loads(runrecord._s(conn.hget(runrecord.run_key(run_id), "stages")) or "[]"))
    return [stage for stage in range(count) if stage_unended(conn, run_id, stage)]


def repair_unended_stages(conn, run_id, wrapper=None):
    """
    OPS.186: finds the stages of the run that finished and did not advance and ends each (`_end_stage`: concludes what was not concluded, queues the opening
    that was not queued, and only that). The retry of a killed unit calls it for its stage (`run_unit`), and the supervisor's reconciler (OPS.210) calls it for
    every live run, so a stage whose end was lost is never waited on for ever.

    Returns:
        list: The indexes of the stages it ended.
    """
    if wrapper is None:
        from planetgen.queue import wrapper             # the wrapper queues units, so it imports this file
    repaired = unended_stages(conn, run_id)
    for stage in repaired:
        _end_stage(conn, run_id, stage, wrapper)
    return repaired


def _end_stage(conn, run_id, stage, wrapper, sleep=time.sleep):
    """
    What the unit that completes `stage` does for the run (PERF.172): the command's conclusion, the stage's
    statistics, and the queuing of the next stage's opening (`enqueue_open`). A failure in the conclusion is the run's failure, said at
    once with the stage named; left to RQ it would read as a unit that failed after it was done, and the run would wait for a stage
    nobody queues.

    OPS.186: idempotent and resumable. `runrecord.ended_key` says how far it got (`ENDED_CONCLUDED`, `ENDED_QUEUED`), so the retry of a unit
    whose horse was killed here (`run_unit`) goes on from there and never concludes twice or queues the opening twice, and the queuing is tried
    `OPEN_TRIES` times before it fails the run (`unit_failed` fails it loudly, with the stage named, when the opening is still not queued). A lease
    keeps a second attempt from ending the stage beside the first.
    """
    lease = runrecord.ending_key(run_id, stage)
    if not conn.set(lease, os.getpid(), nx=True, ex=ENDING_LEASE_SECONDS):
        return                      # another attempt is ending this stage
    try:
        got = planrows.ended(conn, run_id, stage)
        if planrows.RELEASE in got:
            return
        took = 0.0 if planrows.CONCLUDE in got else wrapper.conclude_stage_job(conn, run_id, stage)
        if planrows.FINISH not in got:
            wrapper.finish_stage_job(conn, run_id, stage, took)
        queue_name = wrapper.queue_of(runrecord._s(conn.hget(runrecord.run_key(run_id), "kind")))
        for attempt in range(OPEN_TRIES):
            try:
                enqueue_open(conn, run_id, stage + 1, queue_name)
                break
            except Exception as exc:        # noqa: BLE001 -- a Redis hiccup is tried again; the last one is raised
                if attempt == OPEN_TRIES - 1:
                    raise
                log.warning(f"Run {run_id}: queuing the opening of stage {stage + 1} failed ({type(exc).__name__}: {exc}); trying again.")
                sleep(0.5 * 2 ** attempt)
        planrows.mark_ended(conn, run_id, stage, planrows.RELEASE)
    finally:
        conn.delete(lease)


def open_job_id(run_id, index):
    """The RQ job id of the opening of stage `index` (cancel, PERF.99)."""
    return f"{run_id}-open{index}"


def enqueue_open(conn, run_id, index, queue_name):
    """
    Queues the opening of stage `index` of the run (PERF.254): the cut of the stage into its units, and of each stage after it that has no
    units, or the end of the run, as a job of its own (`open_stage`). The job claims the lane to itself on the memory gate (`claims.CUT_MB`,
    learned per command), has the time limit of a job with no estimate, and is retried once if its horse is killed for memory, with a larger
    claim (`horse_killed`).
    """
    from planetgen.queue import wrapper         # the wrapper queues units, so it imports this file
    if wrapper.forked(conn, run_id):            # OPS.228: the supervisor's launcher hands the opening out as a job
        planrows.queue_opening(conn, run_id, index)
        return
    queue = redisqueue.queue(queue_name, conn)
    queue.enqueue(OPEN_FUNCTION, args=(run_id, index), job_id=open_job_id(run_id, index), job_timeout=redisqueue.job_timeout(None), result_ttl=0,
                  failure_ttl=FAILED_KEEP_SECONDS, retry=rq.Retry(max=1), on_failure=Callback(OPEN_FAILURE_FUNCTION),
                  description=f"open stage {index}")


def open_stage(run_id, index):
    """
    The RQ job of an opening (PERF.254): `wrapper.advance` for stage `index`, which cuts it and queues its units, or ends the run when no stage is
    left. A failure is the run's failure, said at once with the stage named. A run that is no longer live is left alone.

    Returns:
        dict: `opened`, or `skipped`.
    """
    conn = rq.get_current_job().connection
    live, wanted = _run_settings(conn, run_id)
    if not live or conn.exists(runrecord.opened_key(run_id, index)):
        return {"skipped": True}
    return _relayed(conn, run_id, wanted, lambda: _open_stage(conn, run_id, index))


def run_opening_on(conn, run_id, index):
    """
    The opening of stage `index` run on the Redis connection `conn` by a job process the launcher forked (OPS.228): what `open_stage` does for an RQ job. A run that is no longer live, or a stage already opened, is
    left alone. Returns `{"opened": True}` or `{"skipped": True}`.
    """
    live, wanted = _run_settings(conn, run_id)
    if not live or conn.exists(runrecord.opened_key(run_id, index)):
        return {"skipped": True}
    return _relayed(conn, run_id, wanted, lambda: _open_stage(conn, run_id, index))


def _open_stage(conn, run_id, index):
    from planetgen.queue import wrapper         # the wrapper queues units, so it imports this file
    try:
        wrapper.advance(conn, run_id, index)
    except Exception as exc:
        runrecord.finish_run(conn, run_id, "failed", error=f"stage {wrapper.stage_label(conn, run_id, index)} could not "
                                                           f"start: {type(exc).__name__}: {exc}")
        raise
    conn.set(runrecord.opened_key(run_id, index), 1, ex=runrecord.KEEP_SECONDS)     # OPS.187: from here a failure RQ reports for this job is not the run's
    log.flush_relayed()
    return {"opened": True}


def open_failed(job, connection, exc_type, exc_value, _traceback):
    """RQ's `on_failure` for an opening: a failure that was not the cut's own (the job could not start, say) fails the run with the stage named,
    unless the run has ended already."""
    from planetgen.queue import wrapper             # the wrapper queues units, so it imports this file
    run_id, index = job.args[:2]
    if connection.exists(runrecord.opened_key(run_id, index)):       # OPS.187: the cut was finished; RQ's own housekeeping failed after it
        log.debug(f"The opening of stage {index} of run {run_id} was done when RQ failed it ({exc_type.__name__}: {exc_value}); the run goes on.")
        return
    runrecord.finish_run(connection, run_id, "failed",
                         error=f"stage {wrapper.stage_name(connection, run_id, index)} could not start: {exc_type.__name__}: {exc_value}")
    requires.job_failed(job, connection, exc_type, exc_value)


def _short(value, limit=300):
    text = json.dumps(value, default=str)
    return text if len(text) <= limit else text[:limit] + f"... ({len(text)} chars)"


TABLE_KINDS = frozenset({"scatter-mass", "scatter-luminosity", "scatter-backfill", "scatter-band", "generate-start-sector",
                         "generate-neighborhood", "generate-sectors", "generate-unplaced-sector"})
"""frozenset: The commands whose units draw stars from the bright-star sampling tables (PERF.226), so a worker builds the tables
before it forks a horse for one (`warm_unit`) and keeps them; any other unit runs in a worker without them. `generate-system` is not one
(PERF.266): a one-off system draws its star from the generator's own laws and never reads a sampling table, and building them all at the
worker's start took 56 s on a cold table cache (a fresh checkout) for nothing."""


def tables_held():
    """Whether this process holds bright-star sampling tables (`star_population._bright_table`)."""
    import sys
    module = sys.modules.get("planetgen.generation.star_population")
    return module is not None and module._bright_table.cache_info().currsize > 0


def sheds_tables_for(job):
    """Whether `job` is a unit that needs no star tables and this worker holds them: the worker starts itself over for it (PERF.226)."""
    if job.func_name != UNIT_FUNCTION:
        return False
    return job.args[3] not in TABLE_KINDS and tables_held()


def warm_star_tables():
    """Builds in the worker the sampling table of every population at the sectors' floor (PERF.42, PERF.226), which each horse it forks
    starts with, for the units that draw stars (`TABLE_KINDS`, `warm_unit`). Opens no connection."""
    from planetgen import tuning
    from planetgen.generation import star_population
    for population in star_population.POPULATIONS:
        star_population.warm_bins(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL, population)


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
    try:
        if kind in TABLE_KINDS:
            warm_star_tables()
        if spec.warm is not None:
            commands.warm_function(kind)(args)
    except Exception as exc:
        log.swallowed("queue/units", "warm_unit", exc)


def warm_shape(unit_args):
    """A `warm` function for the units that read the galaxy's shape first thing (`unit_args["database"]`): reads it, so
    the worker holds its component field (`warm_unit`). The connection is closed and its pool dropped before the fork,
    so no socket is shared with a horse."""
    config = store.MySQLConfig.from_args(unit_args["database"])
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
    kind = job.args[3] if len(job.args) > 3 else "?"
    if fail_unit(connection, run_id, stage, unit, kind, getattr(job, "retries_left", None), exc_type, exc_value, _traceback, row=row_of(job)):
        requires.job_failed(job, connection, exc_type, exc_value)      # PERF.247: a job held for this unit does not wait for it for ever


def fail_unit(connection, run_id, stage, unit, kind, retries_left, exc_type, exc_value, _traceback, row=None):
    """
    What a unit's failure does to its run, whoever ran the unit (`unit_failed` for RQ, `forkedjob.run` for a forked job): the run fails with the unit and the error recorded, unless the
    unit's result is stored and counted already. `row` is the unit's row in the job ledger: a unit that has one is held for a person instead (OPS.221, `held.hold_unit`) and the run goes on.
    Returns whether the unit failed for good (the run was failed or the unit held).
    """
    unqueued = ""
    if connection.sismember(runrecord.done_key(run_id, stage), unit):
        unended = stage_unended(connection, run_id, stage) and runrecord._s(connection.hget(runrecord.run_key(run_id), "state")) == "running"
        if not unended or (retries_left or 0) > 0:
            # OPS.186: with the stage's end unfinished, a retry left is what finishes it (`run_unit` resumes it); with none left the run fails below
            log.debug(f"Unit {unit} of stage {stage} of run {run_id} was done when RQ failed it ({exc_type.__name__}: "
                      f"{exc_value}); the run goes on.")
            return False
        unqueued = " (the stage was finished but the opening of the next could not be queued)"
    from planetgen.queue import wrapper             # the wrapper queues units, so it imports this file
    stage_name = wrapper.stage_name(connection, run_id, stage)
    # OPS.162: the failure reaches error.log with the run, stage, unit and kind beside the traceback (RQ's own copy has none of them)
    log.record_error("queue", f"Unit {unit} of stage {stage_name} failed: {exc_type.__name__}: {exc_value}",
                     exc_info=(exc_type, exc_value, _traceback), run_id=run_id, stage=stage, stage_name=stage_name, unit=unit, kind=kind)
    if row is not None and not unqueued:
        # OPS.221: a unit with a row in the ledger is held for a person; the rest of the run goes on, and the run is held when nothing more can run
        held.hold_unit(connection, run_id, stage, unit, row, f"{exc_type.__name__}: {exc_value}", killed=exc_type.__name__ in KILLED_EXCEPTIONS)
        return True
    runrecord.finish_run(connection, run_id, "failed",
                         error=f"unit {unit} of stage {stage_name} failed{failures.describe(connection, run_id, exc_value)}: "
                               f"{exc_type.__name__}: {exc_value}{unqueued}", killed=exc_type.__name__ in KILLED_EXCEPTIONS)
    return True


def horse_killed(job, retpid, ret_val, _rusage):
    """
    RQ's `work_horse_killed_handler` (PERF.106): the process that was running a job was killed (out of memory, or `kill -9`).

    A job killed with SIGKILL (what the out-of-memory killer sends) that has a retry left is not the end of the run (PERF.254): the memory gate
    notes what the worker claimed and the most it held, learns a larger size for the kind and gives the retry a claim `memorygate.GROW_FACTOR`
    times the old one (`memorygate.grow`), so the retry waits for room that fits; RQ queues it again. A kill with no retry left fails the run
    now, with the unit (or the opening) named and what it claimed, instead of when the worker's heartbeat runs out; the retry RQ makes of a job
    of a run that has ended finds the run no longer live and does nothing. A job that is not part of a run is only learned from and left to RQ.
    """
    from planetgen.queue import memorygate
    conn = job.connection
    status = log.decode_wait_status(ret_val) if isinstance(ret_val, int) else -1
    name = getattr(job, "id", None) or "?"
    note = memorygate.kill_note(conn, getattr(job, "worker_name", None) or "")
    memory_kill = _killed_by_sigkill(ret_val)
    if row_of(job) is not None:         # OPS.223: the attempt could not write its own end; its row says how the horse went
        attempts.close_lost(conn, job.args[0], row_of(job), "killed_memory" if memory_kill else "lost", error_kind="killed",
                            error_text=f"the work horse ({retpid}) ended with {_ended_with(ret_val)}")
    detail = ""
    if note is not None and memory_kill:
        detail = f"; it claimed {format_megabytes(note['claimed_mb'])} on the memory gate and the most seen was {format_megabytes(note['peak_mb'])}"
        bigger = memorygate.grow(conn, note["kind"], job.id, note["claimed_mb"], note["peak_mb"])
        jobprofiles.record(note["kind"], peak_mb=bigger / memorygate.LEARN_MARGIN, outcome="killed_memory")     # OPS.155: a kill is data; its peak is the size the retry was given less the margin (OPS.156 floors the kind there)
        retry = (getattr(job, "retries_left", None) or 0) > 0
        log.warning(f"Job {name} ({job.func_name}) was killed with SIGKILL while it claimed {format_megabytes(note['claimed_mb'])} (most seen {format_megabytes(note['peak_mb'])}); "
                    + (f"it is queued again with a claim of {format_megabytes(bigger)}." if retry else "it has no retry left."))
        if retry:
            if job.func_name == UNIT_FUNCTION:
                conn.delete(runrecord.ending_key(job.args[0], job.args[1]))     # OPS.186: the lease of the end this horse was in, so the retry can finish it
            return
    if job.func_name == UNIT_FUNCTION and _was_ending_its_stage(conn, job):
        # OPS.186: the horse that completed the stage was killed in the end of it, whatever killed it (`kill -9`, the memory killer): the unit is done, the retry RQ makes of it finishes the end
        # (`run_unit` resumes it), and a run is not failed for what the next attempt carries on
        log.warning(f"Job {name} ({job.func_name}) was killed in the end of its stage ({_ended_with(ret_val)}); it is queued again to finish it.")
        conn.delete(runrecord.ending_key(job.args[0], job.args[1]))
        return
    if job.func_name in (UNIT_FUNCTION, OPEN_FUNCTION):
        from planetgen.queue import wrapper             # the wrapper queues units, so it imports this file
        run_id, index = job.args[:2]
        if job.func_name == UNIT_FUNCTION:
            unit = job.args[2]
            doing, what = f"running unit {unit} of run {run_id}", f"unit {unit} of stage {wrapper.stage_name(conn, run_id, index)}"
        else:
            doing, what = f"opening stage {index} of run {run_id}", f"the opening of stage {wrapper.stage_name(conn, run_id, index)}"
        log.process_ended("queue", f"The process {doing}, pid {retpid},", status, deliberate=(), pid=retpid)
        error = f"{what} was killed: its worker process {retpid} ended with {_ended_with(ret_val)}{detail}"
        if job.func_name == UNIT_FUNCTION and row_of(job) is not None and runrecord.run_state(conn, run_id) in runrecord.LIVE_STATES:
            held.hold_unit(conn, run_id, index, unit, row_of(job), error, killed=True)         # OPS.221: held for a person, as a unit that failed
        else:
            runrecord.finish_run(conn, run_id, "failed", error=error, killed=True)
    else:
        log.process_ended("queue", f"The process running job {name} ({job.func_name}), pid {retpid},", status, deliberate=(), pid=retpid)


def _was_ending_its_stage(conn, job):
    """Whether the unit of `job` is done and its stage's end is not (the state a kill in the end leaves), the run is live and the job has a retry left to finish it."""
    run_id, stage, unit = job.args[:3]
    return (bool(conn.sismember(runrecord.done_key(run_id, stage), unit)) and (getattr(job, "retries_left", None) or 0) > 0
            and runrecord._s(conn.hget(runrecord.run_key(run_id), "state")) == "running" and stage_unended(conn, run_id, stage))


def _killed_by_sigkill(status):
    """Whether a work horse's wait status says SIGKILL, what the kernel's out-of-memory killer sends."""
    return _signal_number(status) == signal.SIGKILL


def _signal_number(status):
    """The signal a work horse's wait status says ended it (a negative status is the signal already; a raw wait status keeps it in the low bits),
    0 for a process that exited, `None` for a status that is no number."""
    try:
        number = int(status)
    except (TypeError, ValueError):
        return None
    return -number if number < 0 else number & 0x7f


def _ended_with(status):
    """
    What a work horse's wait status says: the signal that ended it, by name, and what that usually means. SIGKILL is what the
    out-of-memory killer sends; SIGTERM is a request to end, sent by something that stopped it on purpose (a stop-all, a cancel, a
    memory guard such as earlyoom), which is not the kernel's out-of-memory killer. A status that is no signal is the exit code.
    """
    number = _signal_number(status)
    if not number:
        return f"status {status}"
    try:
        name = signal.Signals(number).name
    except ValueError:
        return f"status {status}"
    if number == signal.SIGKILL:
        return f"status {status} ({name}: the machine may be out of memory)"
    if number == signal.SIGTERM:
        return (f"status {status} ({name}: something asked it to end, a stop-all, a cancel or a memory guard on the machine; "
                f"the kernel's out-of-memory killer sends SIGKILL)")
    return f"status {status} ({name})"


def _progress(conn, run_id, stage):
    """`{unit: (done, estimate)}` of the units reporting now (OPS.136: read from the beat record of each)."""
    found = {}
    for unit, raw in conn.hgetall(runrecord.progress_key(run_id, stage)).items():
        record = unitbeat.parse(raw)
        if record is None:
            continue
        try:
            done, estimate = float(record["done"]), float(record["estimate"])
            if estimate > 0.0:        # a unit that has not said how much work it has is on the beat record but has no bar
                found[int(unit)] = (done, estimate)
        except (KeyError, ValueError, TypeError):
            continue
    return found


def _step_now(conn, run_id, stage, seen):
    """
    The second-level count of the stage as one (OPS.192): `{"name", "done", "total"}` for the units still working, from the beat records. The count is the first working unit's
    step; the units that are in a step of the same name add their counts to it (four units each removing their share of a table are one count of the table). `None` where no working
    unit is in a step.
    """
    found = {}
    for unit, raw in conn.hgetall(runrecord.progress_key(run_id, stage)).items():
        record = unitbeat.parse(raw)
        if record is None or not record.get("step"):
            continue
        try:
            number = int(unit)
            if number not in seen:
                found[number] = (str(record["step"]), float(record.get("step_done") or 0.0), float(record.get("step_total") or 0.0))
        except (ValueError, TypeError):
            continue
    if not found:
        return None
    name = found[min(found)][0]
    same = [counts for counts in found.values() if counts[0] == name]
    return {"name": name, "done": sum(counts[1] for counts in same), "total": sum(counts[2] for counts in same)}


def follow_stage(conn, run_id, stage, units, on_result, on_progress=None, ensure=None, clock=time.monotonic,
                 sleep=time.sleep, on_debug=None, touch=True, on_step=None):
    """
    Follows stage `stage` of the run until its `units` units are done: calls `on_result(unit, result)` once for
    each as it finishes and `on_progress(unit, done, estimate)` for those still working, and `ensure()` every
    `ENSURE_SECONDS` (starting workers, replacing dead ones), and `on_debug()` at every poll and when the stage ends
    (a debug run prints the lines its workers relayed, `planetgen.util.log.RelayReader`).

    Args:
        touch (bool): Whether this process says the run is alive (`runrecord.touch`) at every poll. The process that
            started the run does; one that only watches it (`wrapper.follow`, PERF.235) does not, so a run whose own
            process is gone is not taken for a live one on its account.

        on_step (callable, optional): `on_step(step)` at every poll with the stage's second-level count (`_step_now`: `{"name", "done", "total"}`, or `None` when no unit is in one), and once with
            `None` when the stage ends, so a display that drew it can take it down (OPS.192).

    Raises:
        RunStopped: The run is no longer live and the stage is not complete (a unit failed, or it was canceled).
    """
    try:
        _follow(conn, run_id, stage, units, on_result, on_progress, ensure, clock, sleep, on_debug, touch, on_step)
    finally:
        if on_step is not None:
            on_step(None)


def _gate_note(conn, seen, units):
    """` ` and the gate's figures (`memorygate.gate_words`) when units are still to start, so a stage that runs fewer units than it has workers says what limits it (OPS.270); `""` when the gate cannot be read."""
    if seen >= units:
        return ""
    from planetgen.queue import memorygate
    try:
        return f" ({memorygate.gate_words(conn)})"
    except Exception as exc:  # noqa: BLE001 -- the note is an addition to a line that is said anyway
        log.swallowed("queue/units", "_gate_note", exc)
        return ""


def _quiet_line(conn, run_id, stage, units, seen, quiet_for):
    """The line the log gets when a stage has finished no unit for `quiet_for` seconds (OPS.253): how long, how many units are done, and what the running ones are doing."""
    from planetgen.queue import runview, unitwatch, wrapper
    cause = runview.stall_cause(unitwatch.read(conn, run_id, stage))
    return f"Stage '{wrapper.stage_name(conn, run_id, stage)}': no unit has finished for {format_duration_seconds(quiet_for)}; {seen} of {units} units are done" + (f", the others: {cause}." if cause else ".") + _gate_note(conn, seen, units)


def _follow(conn, run_id, stage, units, on_result, on_progress, ensure, clock, sleep, on_debug, touch, on_step=None):
    seen = set()
    last_ensure = float("-inf")
    last_finish = last_said = clock()

    def drain():
        nonlocal last_finish
        done = {int(unit) for unit in conn.smembers(runrecord.done_key(run_id, stage))} - seen
        if done:
            ordered = sorted(done)
            for unit, raw in zip(ordered, conn.hmget(runrecord.results_key(run_id, stage), ordered)):
                on_result(unit, json.loads(runrecord._s(raw)))
            seen.update(ordered)
            last_finish = clock()

    while True:
        if touch:
            runrecord.touch(conn, run_id)
        drain()
        if on_debug is not None:
            on_debug()
        if len(seen) >= units:
            return
        state = runrecord.run_state(conn, run_id)
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
        if on_step is not None:
            on_step(_step_now(conn, run_id, stage, seen))
        now = clock()
        if now - max(last_finish, last_said) >= QUIET_LINE_SECONDS:
            last_said = now
            log.normal(_quiet_line(conn, run_id, stage, units, len(seen), now - last_finish))
        if ensure is not None and clock() - last_ensure >= ENSURE_SECONDS:
            last_ensure = clock()
            ensure()
        sleep(POLL_SECONDS)
