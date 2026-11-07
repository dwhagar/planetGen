# planetgen/queue/requires.py

"""
Job dependencies (PERF.247, docs/design/operations/redis-and-the-work-queue.md 5.13): a job can declare that it requires other jobs,
so no worker is ever handed one that cannot run yet.

A queue entry made with `requires` (`TimedQueue.create_job(..., requires=[job ids])`) is not put in any queue while a required job
is not done: RQ holds it in its deferred registry, so the round robin of a worker that serves several queues (PERF.245) cannot pick
it, and the worker that finishes the last required job queues it, with no polling. Only a declared requirement holds a job: a job
that declares none is queued at once and runs beside every other, as before.

* **Held.** A job whose requirement is not done is logged as held (`held`), with the ids it waits for.
* **Released.** The job that a held one waits for logs the release as it ends (`job_succeeded`), naming the jobs it frees.
* **A failed or stopped requirement.** RQ would leave the dependant deferred for good. `job_failed` (the on-failure callback of every
  job the queue makes, and called by a unit's own, `units.unit_failed`) fails each dependant that is still held, with the line
  `required job X failed`, and so on down the chain; a unit that is a dependant fails its run with that line. `drop_held` removes the
  held jobs of a queue (Stop All, a canceled run) the same way. A requirement declared `optional` is not a barrier to failure: its
  dependant runs after the requirement ends, whichever way it ended.
* **A requirement that is gone** (finished and deleted, which is what `result_ttl=0` does) counts as done: it cannot be told from a
  job that never existed.

Both jobs must use one serializer (the queues of a lane do, `redisqueue.serializer_for`): RQ reads the dependants with the serializer
of the queue the required job ran on.
"""

import importlib

import rq
from rq.job import Dependency, Job, JobStatus
from rq.registry import DeferredJobRegistry, FailedJobRegistry

from planetgen.util import log

FAILED_FUNCTION = "planetgen.queue.requires.job_failed"
"""str: The on-failure callback the queue gives a job that has none of its own."""

SUCCEEDED_FUNCTION = "planetgen.queue.requires.job_succeeded"
"""str: The on-success callback the queue gives a job that has none of its own."""

FAILED_KEEP_SECONDS = 600
"""int: How long a dependant that was failed for its requirement is kept in the failed registry."""


class RequiredJobFailed(RuntimeError):
    """The error a held job is failed with: a job it required failed or was stopped."""


def job_ids(requires):
    """
    The ids of `requires` (RQ jobs or ids), in order, once each.

    Raises:
        ValueError: A requirement that is empty.
    """
    ids = []
    for item in requires or ():
        found = item.id if isinstance(item, Job) else str(item or "")
        if not found:
            raise ValueError("a requirement needs a job id")
        if found not in ids:
            ids.append(found)
    return ids


def dependency(requires, optional=False):
    """
    The RQ dependency of a job that requires `requires` (RQ jobs or ids), or `None` for no requirement.

    Args:
        optional (bool): The job runs after the requirements end whichever way they end (`allow_failure`); by default a requirement
            that failed or was stopped fails it.
    """
    ids = job_ids(requires)
    return Dependency(jobs=ids, allow_failure=bool(optional)) if ids else None


def enqueue(queue, function, *args, requires=(), optional=False, **options):
    """
    Queues `function(*args)` on `queue` (a `redisqueue.queue`) held until the jobs `requires` (RQ jobs or ids) are done, or at once
    when there is none; `options` are RQ's (`job_timeout`, `job_id`, `result_ttl`, ...).

    Returns:
        rq.job.Job: The job; its status is `deferred` while it is held.

    Raises:
        ValueError: A requirement that is empty, is the job itself or was queued with another serializer.
    """
    wanted = dependency(requires, optional)
    if wanted is not None:
        options["depends_on"] = wanted
    return queue.enqueue(function, *args, **options)


def check_serializers(connection, ids, serializer, own_id=None):
    """
    Refuses a requirement that RQ could not read the dependants of: a required job queued with another serializer than the dependant's.

    Raises:
        ValueError: A job requires itself, or a required job uses another serializer.
    """
    from planetgen.queue import redisqueue
    for found in ids:
        if own_id is not None and found == own_id:
            raise ValueError(f"job {found} cannot require itself")
        origin = connection.hget(f"rq:job:{found}", "origin")
        if origin is None:
            continue            # done and gone, or not queued yet: nothing to read
        if redisqueue.serializer_for(origin.decode() if isinstance(origin, bytes) else origin) is not serializer:
            raise ValueError(f"job {found} was queued with another serializer than the job that requires it; "
                             f"a requirement and its dependant must share one (the queues of a lane do)")


def log_held(job, ids):
    """Writes that `job` is held until the jobs `ids` are done: to the log, and to its run's when it is a unit."""
    text = f"held: job {job.id} waits for {', '.join(ids)}"
    _say(job.connection, job, text, "job.held")


def _say(connection, job, text, kind):
    log.normal(text)
    run = _run_of(job)
    if run is not None:
        from planetgen.queue import runrecord
        try:
            runrecord.log(connection, run[0], text, level="info", kind=kind, stage=run[1])
        except Exception:  # noqa: BLE001 -- a line in the run's log must never fail the job
            pass


def _run_of(job):
    """`(run id, stage)` of a unit of a run, `None` for any other job."""
    from planetgen.queue import units
    if getattr(job, "func_name", None) != units.UNIT_FUNCTION or len(job.args) < 2:
        return None
    return job.args[0], job.args[1]


def _dependents(connection, job_id):
    return sorted(item.decode() if isinstance(item, bytes) else item for item in connection.smembers(f"rq:job:{job_id}:dependents"))


def _held_job(connection, dependent_id):
    """`(job, origin)` of a dependant that is still held, `None` when it is not (queued, running, done, or gone)."""
    from planetgen.queue import redisqueue
    origin = connection.hget(f"rq:job:{dependent_id}", "origin")
    if origin is None:
        return None
    origin = origin.decode() if isinstance(origin, bytes) else origin
    job = redisqueue.fetch_job(dependent_id, connection, redisqueue.serializer_for(origin))
    if job is None or job.get_status(refresh=True) != JobStatus.DEFERRED:
        return None
    return job, origin


def job_succeeded(job, connection, result, *_ignored):
    """RQ's `on_success` for a job that has none of its own: writes which held jobs this one's end frees (all their other requirements
    done too) or leaves still held."""
    try:
        for dependent_id in _dependents(connection, job.id):
            held = _held_job(connection, dependent_id)
            if held is None:
                continue
            others = [other for other in held[0].dependency_ids if other != job.id and connection.exists(f"rq:job:{other}")]
            if others:
                _say(connection, held[0], f"held: job {dependent_id} still waits for {', '.join(others)} after job {job.id} finished",
                     "job.held")
            else:
                _say(connection, held[0], f"released: job {dependent_id} was waiting for job {job.id}, which finished", "job.released")
    except Exception as exc:  # noqa: BLE001 -- the log must never fail the job that finished
        log.debug(f"Release of the dependants of job {job.id} not logged ({exc}).")


def job_failed(job, connection, exc_type=None, exc_value=None, _traceback=None):
    """RQ's `on_failure` for a job that has none of its own (and called by a unit's own): fails the jobs held for this one."""
    job_id = getattr(job, "id", None)
    if not job_id:
        return []
    reason = f"{getattr(exc_type, '__name__', 'an error')}: {exc_value}" if exc_type is not None else "it failed"
    try:
        return fail_dependents(connection, job_id, f"failed ({reason})")
    except Exception as exc:  # noqa: BLE001 -- failing what waits must never hide the failure being reported
        log.warning(f"The jobs held for job {job_id} were not failed: {exc}")
        return []


def fail_dependents(connection, job_id, why, _seen=None):
    """
    Fails every job held for the required job `job_id` (it `why`: failed, or was stopped) and, down the chain, those held for them:
    each is taken out of the deferred registry, put in the failed registry with `required job X <why>`, and, when it is a unit or
    has an on-failure callback of its own, that callback is called, so a run that a held unit belongs to fails with the line. A
    dependant that is `optional` is left to RQ, which queues it.

    Returns:
        list[str]: The ids of the jobs failed.
    """
    seen = _seen if _seen is not None else set()
    failed = []
    for dependent_id in _dependents(connection, job_id):
        if dependent_id in seen:
            continue
        seen.add(dependent_id)
        held = _held_job(connection, dependent_id)
        if held is None:
            continue
        job, origin = held
        if getattr(job, "allow_dependency_failures", False):
            continue
        line = f"required job {job_id} {why}"
        _fail_held(connection, job, origin, line)
        failed.append(dependent_id)
        failed.extend(fail_dependents(connection, dependent_id, "did not run: " + line, seen))
    return failed


def _fail_held(connection, job, origin, line):
    from planetgen.queue import redisqueue
    serializer = redisqueue.serializer_for(origin)
    DeferredJobRegistry(origin, connection=connection, serializer=serializer).remove(job)
    job.set_status(JobStatus.FAILED)
    job._exc_info = line
    job.save()
    FailedJobRegistry(origin, connection=connection, serializer=serializer).add(job, ttl=job.failure_ttl or FAILED_KEEP_SECONDS,
                                                                                 exc_string=line)
    _say(connection, job, f"dropped: job {job.id} will not run: {line}", "job.dropped")
    callback = getattr(job, "failure_callback_name", None) or getattr(job, "_failure_callback_name", None)
    if callback and callback != FAILED_FUNCTION:
        module, _, name = callback.rpartition(".")
        getattr(importlib.import_module(module), name)(job, connection, RequiredJobFailed, RequiredJobFailed(line), None)


def drop_held(connection, queue_name):
    """
    Removes the held jobs of queue `queue_name` (Stop All, a canceled run): each is deleted with the entries that held for it, so
    nothing waits for a requirement that will not run.

    Returns:
        int: How many held jobs were dropped.
    """
    from planetgen.queue import redisqueue
    serializer = redisqueue.serializer_for(queue_name)
    registry = DeferredJobRegistry(queue_name, connection=connection, serializer=serializer)
    dropped = 0
    for job_id in registry.get_job_ids():
        job = redisqueue.fetch_job(job_id, connection, serializer)
        registry.remove(job_id)
        if job is not None:
            job.delete()
        dropped += 1
    if dropped:
        log.normal(f"dropped: {dropped} job(s) held in {queue_name} for a requirement")
    return dropped


def held_jobs(connection, queue_name):
    """The jobs held in queue `queue_name`: `[(job id, [the job ids it waits for])]`, for the Queue page and the log."""
    from planetgen.queue import redisqueue
    registry = DeferredJobRegistry(queue_name, connection=connection, serializer=redisqueue.serializer_for(queue_name))
    held = []
    for job_id in registry.get_job_ids():
        job = redisqueue.fetch_job(job_id, connection, redisqueue.serializer_for(queue_name))
        if job is not None:
            held.append((job_id, list(job.dependency_ids)))
    return held


HELD_SHOWN = 100
"""int: Most held jobs the Queue page lists (the count it gives is the whole count)."""


def held_overview(connection, limit=HELD_SHOWN):
    """
    The jobs held in any queue, for the Queue page (PERF.248): what each waits for, and which run and stage it belongs to when it is a unit.

    Returns:
        dict: `{"total": how many are held, "items": [{"id", "queue", "label", "waits_for": [job ids still not done], "run", "stage"}]}`,
        the first `limit` ones, queue by queue.
    """
    from planetgen.queue import redisqueue
    items = []
    total = 0
    for name in redisqueue.QUEUES:
        registry = DeferredJobRegistry(name, connection=connection, serializer=redisqueue.serializer_for(name))
        ids = registry.get_job_ids()
        total += len(ids)
        for job_id in ids:
            if len(items) >= limit:
                break
            job = redisqueue.fetch_job(job_id, connection, redisqueue.serializer_for(name))
            if job is None:
                continue
            run = _run_of(job)
            items.append({"id": job_id, "queue": name, "label": getattr(job, "func_name", None) or "",
                          "waits_for": [other for other in job.dependency_ids if connection.exists(f"rq:job:{other}")],
                          "run": run[0] if run else None, "stage": run[1] if run else None})
    return {"total": total, "items": items}
