# planetgen/queue/jobbeat.py

"""
The beat of the jobs that are not units of a run (OPS.142): an API job (`queue/api_jobs.py`), a web job (`web/job_runner.py`, the Generate page's) and a tile prefetch chunk
(`web/lib/tileprefetch.py`). They relied on RQ's time limit alone; now each writes a beat while it runs, the way a unit does (`unitbeat.UnitBeat`: a thread of the process, the
Redis server clock), and the Queue page lists them with the state `beatjudge` gives.

* **One hash, `planetgen:jobbeats`**, a field to a running job (`<kind>:<job id>`), removed when the job ends. A job that is killed leaves its record, which then reads `dead` for
  `FORGET_SECONDS` (so the page says a job died) and is dropped by the next listing after that.
* **No counter.** These jobs have no unit of work to count (a web job is steps of unknown length, an API job one call), so the judge can say only `working` or `dead`; what the job is
  doing is its `phase` (the step of a web job). A job is never ended on its state: RQ's limit stays what ends one, the page is how a person sees a dead one.
* **Fails open.** A Redis that does not answer means no beat, never a job that does not run.

    with jobbeat.beating("api", job_id, label="regenerate_sector") as beat:
        beat.phase("step 2 of 4")
"""

import time
from contextlib import contextmanager

from planetgen.queue import beatjudge, redisqueue, runrecord, unitbeat

FORGET_SECONDS = 3600.0
"""float: How long a dead job's record is listed before it is dropped."""


def hash_key():
    return runrecord._k("jobbeats")


class _NoBeat:
    """What a job gets where no Redis answers: every call does nothing."""

    def __call__(self, done, estimate):
        pass

    def phase(self, text):
        pass

    def waiting(self, tag):
        return _null()


@contextmanager
def _null():
    yield


@contextmanager
def beating(kind, job_id, label="", conn=None, beat_seconds=None):
    """
    The job `job_id` of `kind` (`api`, `web`, `tile-prefetch`) is running inside the block: a thread writes its beat, and the record goes when the block ends. Yields the `UnitBeat`
    (`phase(text)`), or a stand-in that does nothing where Redis does not answer.
    """
    try:
        conn = conn or redisqueue.connect(redisqueue.redis_url())
        field = f"{kind}:{job_id}"
        beat = unitbeat.UnitBeat(conn, "", 0, field, hash_key=hash_key(), beat_seconds=beat_seconds, extra={"kind": kind, "job": str(job_id), "label": str(label)})
    except Exception:       # noqa: BLE001 -- a job runs without a beat before it runs with a failing Redis
        yield _NoBeat()
        return
    try:
        with beat:
            yield beat
    finally:
        try:
            conn.hdel(hash_key(), field)
        except Exception:       # noqa: BLE001 -- the record expires into `dead` and is dropped
            pass


def states(conn, now=None):
    """
    The running jobs and what each is doing: `[{kind, job, label, state, reason, phase, running_seconds, beat_age_seconds}]`, oldest first. A dead job older than `FORGET_SECONDS`
    since its last beat is dropped from the hash.
    """
    now = unitbeat.server_time(conn) if now is None else now
    rows, forget = [], []
    for field, raw in conn.hgetall(hash_key()).items():
        found = unitbeat.parse(raw)
        if found is None or "beat" not in found:
            continue
        try:
            beat = unitbeat.beat_of(found)
        except (ValueError, TypeError, KeyError):
            continue
        if now - beat.beat > FORGET_SECONDS:
            forget.append(field)
            continue
        verdict = beatjudge.judge(beat, beatjudge.History(), now)
        rows.append({"kind": str(found.get("kind", "")), "job": str(found.get("job", "")), "label": str(found.get("label", "")), "state": verdict.state,
                     "reason": verdict.reason, "phase": beat.phase, "running_seconds": max(now - beat.started, 0.0), "beat_age_seconds": max(now - beat.beat, 0.0)})
    if forget:
        conn.hdel(hash_key(), *forget)
    return sorted(rows, key=lambda row: -row["running_seconds"])
