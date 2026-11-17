# planetgen/queue/lightjobs.py

"""
The site's jobs as light jobs (OPS.230, D5 of `/mnt/project-files/reports/job-engine-items-stage-3-6-2026-11-13.md`; design section 3.2 of `/mnt/project-files/reports/job-engine-design-2026-11-13.md`).

An API job (`api_jobs.submit`), a slow page's answer (`pagejobs.ensure`) and a visitor's node prefetch are not units of a run: they are re-derivable (a page asks again, `singleflight` joins
the duplicates) and must add no MySQL write to the site's fast path. So they are not in the ledger, and no queue holds them. This module is their whole life:

* **A record in Redis** (`job_key`): what to call (`module:name`, the API's own functions, `allowlist.API_FUNCTIONS`), its arguments, the lane (`interactive`, or `long` when it is expected to run
  longer than `queue.long_job_seconds`), when it was queued, started and ended, and, once it ends, its result or its error. Kept `RESULT_TTL_SECONDS` after the end. The list of those not yet
  started is `pending_key` (oldest first), and of those running `running_key`.
* **The supervisor's duty `LightJobs`** takes the pending jobs oldest first, asks the dispatch rule (`dispatch.dispatch`, which serves the site's lanes before every other and keeps their places
  clear) which go now, and forks each from the `plain` template (`jobhost.TemplateSet`) running `run`. It ends what ran past its time limit (OPS.251: `LEARNED_MULTIPLE` times the longest of the kind's own recorded runs once it has `LEARNED_SAMPLES` of them; until then ten times what it was expected to take, as RQ's limit was,
  `time_limit`), settles a job whose process ended without saying how it ended, and drops a job nobody is waiting for.
* **A deadline** (W12 of the audit): a job that says how long it stays wanted (`wanted_seconds`: a page's answer) is dropped, if it has not started, once nobody has asked about it for that long.
  Every `status` is such an ask. A job that began is never dropped: it ends, and its answer is kept for the next visitor. A job that says nothing (an edit, a regeneration) is never dropped.
* **The job's own process writes its end.** `run` writes the result (or the failure, or a refusal the work raised on purpose) to the record before it exits, so a supervisor that ends meanwhile
  loses nothing; one that finds a job whose process is gone with no end written (`recover`) fails it, naming how it died.

A caller sees what it saw with RQ: `submit`, `status`, `wait` and `cancel` of `api_jobs`, with the same states (`queued`, `running`, `succeeded`, `failed`).
"""

import json
import math
import os
import secrets
import signal
import re
import time
import traceback
import types

from planetgen.queue import dispatch, jobhost, memorygate, redisqueue, runrecord, template, work
from planetgen.util import log

RESULT_TTL_SECONDS = 86400
"""int: How long a finished job's result or error is kept."""

UNESTIMATED_SECONDS = 3600
"""int: The time limit of a job with no estimate (PERF.144): a one-system edit, a regenerate of one sector. A fixed ceiling, not an estimate."""

ESTIMATED_CEILING_SECONDS = 7 * 86400
"""int: The longest limit of a job with an estimate (a neighbourhood of a great many sectors can run for days)."""

LEARNED_SAMPLES = 8
"""int: How many finished runs of a kind (`jobprofiles.seconds_of`) its time limit needs before the limit is the kind's own figure and not the typed one (OPS.251). Stand-in."""

LEARNED_MULTIPLE = 3.0
"""float: The learned time limit is this many times the longest of the kind's recorded runs. Stand-in; OPS.249 (H4) learns the edges of the controllers from recorded runs, this one is for it to take too."""

LEARNED_FLOOR_SECONDS = 30.0
"""float: The least a learned limit is, so a kind that always finishes in a blink is not ended by one slow moment of the host. A clamp, not an estimate."""

LEARNED_KEEP_SECONDS = 300.0
"""float: How long the supervisor keeps what it read of a kind, so a start reads the database once in this time and not at every job."""

JOB_ID_RE = re.compile(r"^[0-9a-f]{16}$")
"""re.Pattern: A job id. Checked before any id from a URL reaches Redis."""

CANCELED_TEXT = "the job ended without finishing: it was canceled before it started"
"""str: The error of a job taken off the pending list by `cancel`."""

DROPPED_TEXT = "the job was dropped: nobody was waiting for its answer any more"
"""str: The error of a job `LightJobs` dropped for want of anyone waiting."""

LANE_INTERACTIVE, LANE_LONG = memorygate.LANE_INTERACTIVE, memorygate.LANE_LONG

OWNER = "light"
"""str: What the supervisor's template set calls the jobs this duty forked, so the launcher's poll does not take their ends (`jobhost.TemplateSet.poll`)."""

RUNNER = "planetgen.queue.lightjobs.run"
"""str: What a forked light job runs in its process."""

KIND = "light"
"""str: The kind of command the template set is asked to fork (`jobhost.template_key`): every light job is forked from the `plain` template."""

WINDOW = 64
"""int: How many pending jobs a turn reads (the oldest; `dispatch` judges as many)."""

RECOVER_SECONDS = 30.0
"""float: How often the jobs left `running` by a supervisor that ended are looked for."""

NO_PID_SECONDS = 20.0
"""float: How long a job that is `running` with no process named is waited for before it is called lost: the supervisor marks it and then forks it."""

LANE_MAINTENANCE = memorygate.LANE_MAINTENANCE

PROCESS_NICE = {LANE_INTERACTIVE: 5, LANE_LONG: work.WORKER_NICENESS, LANE_MAINTENANCE: work.WORKER_NICENESS}
"""dict: How much `os.nice` lowers the process of a light job of each lane: the visitor's job a little (the short jobs get the CPU before a generation run's units), a long job and the upkeep a visit asks for (the
node prefetch, `LANE_MAINTENANCE`) as a bulk worker is."""


class NoQueue(Exception):
    """No Redis server answers, or an update is running, so nothing can be queued."""


# -- the keys ------------------------------------------------------------------------------------------------------------------

def job_key(job_id):
    """The hash of a light job's record."""
    return runrecord._k("light", "job", job_id)


def pending_key():
    """The sorted set of the jobs not started yet, scored by when they were queued."""
    return runrecord._k("light", "pending")


def running_key():
    """The set of the ids of the jobs a supervisor has started and not seen end."""
    return runrecord._k("light", "running")


def _text(value):
    return value.decode("utf-8") if isinstance(value, bytes) else value


def _record(connection, job_id):
    """The record of `job_id` as text, or `None` for an unknown or expired id."""
    if not isinstance(job_id, str) or not JOB_ID_RE.match(job_id):
        return None
    found = connection.hgetall(job_key(job_id))
    return {_text(name): _text(value) for name, value in found.items()} if found else None


def time_limit(expected_seconds=None):
    """The seconds a job may run before the supervisor ends it (PERF.144): `redisqueue.job_timeout` of `expected_seconds`, at most `ESTIMATED_CEILING_SECONDS`, and `UNESTIMATED_SECONDS` for a job
    with no estimate."""
    return redisqueue.job_timeout(expected_seconds, ceiling=ESTIMATED_CEILING_SECONDS, unestimated=UNESTIMATED_SECONDS)


def learned_limit(seconds):
    """The time limit a kind's own recorded run times give (OPS.251), or `None` while it has fewer than `LEARNED_SAMPLES` of them: `LEARNED_MULTIPLE` times the longest, at least `LEARNED_FLOOR_SECONDS`.
    The typed limit (`time_limit`) stays what a kind with no history is held to."""
    if len(seconds) < LEARNED_SAMPLES:
        return None
    return max(LEARNED_FLOOR_SECONDS, LEARNED_MULTIPLE * max(seconds))


def is_long(expected_seconds):
    """Whether a job expected to take `expected_seconds` goes to the long lane (PERF.151): more than `queue.long_job_seconds`. A job with no estimate counts as short."""
    from planetgen.util.settings import get_settings
    if expected_seconds is None or not math.isfinite(expected_seconds):
        return False
    return expected_seconds > get_settings().queue.long_job_seconds


def lane_of(expected_seconds=None):
    """The lane of a job expected to take `expected_seconds`: `long` when `is_long`, else `interactive`."""
    return LANE_LONG if is_long(expected_seconds) else LANE_INTERACTIVE


# -- the submitter's side ------------------------------------------------------------------------------------------------------

def add_job(connection, reference, args, expected_seconds=None, wanted_seconds=None, drop_keys=(), url=None, lane=None, limit_seconds=None):
    """
    Records a job and puts it on the pending list. The supervisor is asked for (`ask_supervisor`) by the caller (`api_jobs.submit`).

    Args:
        connection (redis.Redis): The server the supervisor reads.
        reference (str): `module:name` of the API function (`allowlist.API_FUNCTIONS`).
        args (list): The JSON arguments (`api_jobs.encode_argument`).
        expected_seconds (float, optional): How long it is expected to take: the lane and the time limit follow.
        wanted_seconds (float, optional): How long the job stays wanted once nobody asks about it (the deadline); a job without it is never dropped.
        drop_keys (iterable): Keys the supervisor deletes when it drops the job (the page job's claim on its request).
        url (str, optional): The Redis server the job process writes its end to; the default one.
        lane (str, optional): The lane the job goes by when it is not the one its expected time gives (`LANE_MAINTENANCE`: the node prefetch is upkeep that serves a visitor).
        limit_seconds (int, optional): The time limit when it is not `time_limit(expected_seconds)` (a chunk of nodes has a floor of its own).

    Returns:
        str: The job's id.
    """
    job_id = secrets.token_hex(8)
    now = round(time.time(), 3)
    lane = lane or lane_of(expected_seconds)
    limit = limit_seconds or time_limit(expected_seconds)
    fields = {"function": reference, "args": json.dumps(args), "lane": lane, "state": "queued", "created": f"{now:.3f}", "limit": f"{limit:.0f}",
              "expected": "" if expected_seconds is None else f"{float(expected_seconds):.3f}", "wanted": "" if not wanted_seconds else f"{float(wanted_seconds):.3f}",
              "wanted_until": "" if not wanted_seconds else f"{now + float(wanted_seconds):.3f}", "drop_keys": json.dumps(list(drop_keys)), "url": url or ""}
    pipe = connection.pipeline(transaction=True)
    pipe.hset(job_key(job_id), mapping=fields)
    pipe.expire(job_key(job_id), RESULT_TTL_SECONDS + int(limit))
    pipe.zadd(pending_key(), {job_id: now})
    pipe.execute()
    return job_id


def pending_count(connection, reference, scan=1000):
    """How many of the oldest `scan` pending jobs call `reference` (the node prefetch's depth rule: nothing is added while enough of its jobs wait)."""
    found = [_text(item) for item in connection.zrange(pending_key(), 0, scan - 1)]
    if not found:
        return 0
    pipe = connection.pipeline(transaction=False)
    for job_id in found:
        pipe.hget(job_key(job_id), "function")
    return sum(1 for name in pipe.execute() if name is not None and _text(name) == reference)


def status(connection, job_id, ask=True):
    """
    One job's state (`api_jobs.status`).

    Args:
        ask (bool): The call is a visitor asking after the job, which keeps a job that has a deadline wanted (`wanted_until`) and starts a supervisor for a job that waits.

    Returns:
        dict | None: `id`, `state`, `result`, `error`, `error_status` and `waiting`; `None` for an unknown or expired id.
    """
    record = _record(connection, job_id)
    if record is None:
        return None
    state = record.get("state", "failed")
    if ask and state == "queued":
        if record.get("wanted"):
            connection.hset(job_key(job_id), "wanted_until", f"{time.time() + float(record['wanted']):.3f}")
        ask_supervisor(connection, record.get("url") or None)        # a supervisor that was about to end when the job arrived leaves it waiting: whoever asks starts another
    body = {"id": job_id, "state": state, "result": None, "error": None, "error_status": None, "waiting": None}
    if state == "running" and record.get("waiting"):
        body["waiting"] = json.loads(record["waiting"])
    elif state == "succeeded":
        body["result"] = json.loads(record["result"]) if record.get("result") else None
    elif state == "failed":
        body["error"] = record.get("error") or "the job ended without finishing"
        body["error_status"] = int(record["error_status"]) if record.get("error_status") else None
    return body


def waited(connection, job_id):
    """Why a job that waits has not started (what `dispatch` said last), or `None`."""
    record = _record(connection, job_id)
    return None if record is None or record.get("state") != "queued" else (record.get("waited") or None)


def cancel(connection, job_id):
    """
    Takes job `job_id` off the pending list if the supervisor has not taken it yet (PERF.105).

    Returns:
        bool: Whether it was taken off; a job already running, finished or unknown is left alone.
    """
    record = _record(connection, job_id)
    if record is None or record.get("state") != "queued":
        return False
    if not connection.zrem(pending_key(), job_id):
        return False        # the supervisor took it between the two reads
    finish(connection, job_id, "failed", error=CANCELED_TEXT)
    return True


def request_stop(connection, job_id):
    """Asks the supervisor to end job `job_id` now (an update that waits for it, an admin's stop): a job not started is taken off the pending list, one that runs is killed at the supervisor's next turn.
    Returns whether the job was found live."""
    record = _record(connection, job_id)
    if record is None or record.get("state") not in ("queued", "running"):
        return False
    if record["state"] == "queued":
        return cancel(connection, job_id)
    connection.hset(job_key(job_id), "stop", "1")
    return True


def stop_all(connection):
    """Takes every pending job off the list and asks for every running one to stop (`request_stop`): what Stop All does to the site's jobs. Returns how many jobs it found live."""
    return sum(request_stop(connection, item["id"]) for item in live(connection))


def live(connection):
    """The jobs queued or running now: `[{"id", "state", "function", "lane", "created", "started"}]`, oldest first."""
    found = []
    for job_id in [_text(item) for item in connection.zrange(pending_key(), 0, -1)] + [_text(item) for item in connection.smembers(running_key())]:
        record = _record(connection, job_id)
        if record is None or record.get("state") not in ("queued", "running"):
            continue
        found.append({"id": job_id, "state": record["state"], "function": record.get("function", ""), "lane": record.get("lane", ""),
                      "created": float(record.get("created") or 0), "started": float(record["started"]) if record.get("started") else None})
    return sorted(found, key=lambda item: item["created"])


def ask_supervisor(connection=None, url=None):
    """Starts a supervisor unless one is alive (`supervisor.ensure_started`): it is the one that hands the job out. A test patches this to run one turn of `LightJobs` in its own process."""
    from planetgen.queue import supervisor
    supervisor.ensure_started(connection, url)


# -- the end of a job ----------------------------------------------------------------------------------------------------------

def finish(connection, job_id, state, result=None, error=None, error_status=None):
    """Writes the end of job `job_id`: `succeeded` with its `result`, or `failed` with its `error` (and the HTTP status of a refusal). Kept `RESULT_TTL_SECONDS`."""
    fields = {"state": state, "finished": f"{time.time():.3f}"}
    if result is not None:
        fields["result"] = redisqueue.JsonSerializer.dumps(result).decode("utf-8")
    if error is not None:
        fields["error"] = str(error)
    if error_status is not None:
        fields["error_status"] = str(int(error_status))
    pipe = connection.pipeline(transaction=True)
    pipe.hset(job_key(job_id), mapping=fields)
    pipe.hdel(job_key(job_id), "waiting")
    pipe.expire(job_key(job_id), RESULT_TTL_SECONDS)
    pipe.srem(running_key(), job_id)
    pipe.execute()


def _last_line(text):
    lines = [line for line in (text or "").strip().splitlines() if line.strip()]
    return lines[-1] if lines else "the job ended without finishing"


# -- the job's own process -----------------------------------------------------------------------------------------------------

_current = {"id": None, "url": None}


def current_id():
    """The id of the light job this process runs, or `?` outside one."""
    return _current["id"] or "?"


def report_waiting(waiting):
    """Shows that the job this process runs waits in line for the machine (`wrapper.report_wait`: `{"run", "position", "text"}`, or `None` to clear it); `api_jobs.status` reports it as `waiting`. Does nothing
    outside a light job."""
    if _current["id"] is None:
        return
    try:
        connection = redisqueue.connect(_current["url"])
        if waiting is None:
            connection.hdel(job_key(_current["id"]), "waiting")
        else:
            connection.hset(job_key(_current["id"]), "waiting", json.dumps(waiting))
    except Exception as exc:  # noqa: BLE001 -- showing the wait is a nicety
        log.swallowed("queue/lightjobs", "report_waiting", exc)


def run(job, link):
    """
    What a forked light job process runs (`template.run_job`): the API function the job names through `api_jobs.execute` (the allow-list is checked there), with its end written to the record before the process
    exits. `job` is `{"id", "function", "args", "redis_url", "nice"}`. A refusal the work raised on purpose is a failure carrying its HTTP status; any other exception is a failure with the last line of its text.
    Returns the exit code of the process.
    """
    from planetgen.queue import api_jobs
    _current.update(id=job["id"], url=job.get("redis_url"))
    try:
        os.nice(int(job.get("nice") or 0))
    except OSError:
        pass
    connection = redisqueue.connect(job.get("redis_url"))
    try:
        outcome = api_jobs.execute(job["function"], job.get("args") or [])
        if "refused" in outcome:
            finish(connection, job["id"], "failed", error=outcome["refused"], error_status=outcome.get("status"))
        else:
            finish(connection, job["id"], "succeeded", result=outcome.get("result"))
        return 0
    except template.SupervisorGone:
        raise
    except BaseException as exc:  # noqa: BLE001 -- a job's failure is its own: written, then the process ends failed
        text = _last_line("".join(traceback.format_exception_only(type(exc), exc)))
        log.swallowed("queue/lightjobs", f"the light job {job.get('function')}", exc)
        try:
            finish(connection, job["id"], "failed", error=text)
        except Exception as writing:  # noqa: BLE001 -- the supervisor settles a job whose end was not written
            log.swallowed("queue/lightjobs", "run", writing)
        return template.EXIT_FAILED


# -- the supervisor's duty -----------------------------------------------------------------------------------------------------

class Started:
    """A light job this duty forked and has not seen end."""

    def __init__(self, job_id, job, key, claim, now):
        self.id, self.job, self.key, self.claim, self.started = job_id, job, key, claim, now
        self.limit = None
        self.peak_mb = 0.0
        self.killed = None


class LightJobs:
    """
    The duty of the supervisor that hands out the light jobs (see the module). Each turn: what ended is settled, what ran past its limit or was asked to stop is killed, what nobody waits for is dropped,
    the pending jobs are judged by `dispatch.dispatch` against the host and every claim held (the jobs this duty runs and `also_held`, the launcher's), and the admitted ones are forked.

    Args:
        connection (redis.Redis): The server the records are on.
        forker (jobhost.TemplateSet, optional): The supervisor's templates, shared with the launcher; one of this duty's own where left out (`fork(kind, runner, job, owner=)`, `poll(owner=)`, `stop()`).
        read_host (callable, optional): `read_host() -> memorygate.Host`.
        clock (callable, optional): `clock() -> Unix seconds`.
        also_held (callable, optional): `also_held() -> [memorygate.Holder]`, the claims of the jobs that other duty holds.
        held_mb (callable, optional): `held_mb(pid) -> MB` a process and its children hold; the simulator and the tests pass their own.
        runner (str, optional): What the job process runs; `extra` (dict) is added to the job it is given: the tests' runner needs the prefix of the test's keys.
    """

    def __init__(self, connection, forker=None, read_host=memorygate.read_host, clock=time.time, also_held=None, held_mb=memorygate.held_mb, runner=RUNNER, extra=None, hold=None):
        self.connection, self.forker, self.read_host, self.clock = connection, forker, read_host, clock
        self.runner, self.extra = runner, dict(extra or {})
        self.also_held = also_held or (lambda: [])
        self.hold = hold or (lambda: None)      # the memory controller's `holding` (OPS.231): the site's own lanes still go, the upkeep a visit asked for (the node prefetch) waits
        self.held_mb = held_mb
        self.running = {}
        self.learned = {}       # kind -> (when read, the limit its recorded runs give or None): OPS.251
        self.recovered_at = None
        self.supervisor = None

    def attach(self, supervisor):
        self.supervisor = supervisor

    def _act(self, what, reason=""):
        if self.supervisor is not None:
            self.supervisor.act(what, reason)

    def _forker(self):
        if self.forker is None:
            self.forker = jobhost.TemplateSet(background=True)
        return self.forker

    # -- what ended --------------------------------------------------------------------------------------------------------

    def ended(self):
        """Settles every light job whose process ended: one that wrote its end is as it wrote it; one that did not is failed, naming how its process died (or that it ran past its limit). Returns the ids."""
        if self.forker is None:
            return []
        done = []
        for key, exit_ in self.forker.poll(owner=OWNER):
            found = next((job_id for job_id, started in self.running.items() if started.job.id == exit_.id and started.key == key), None)
            if found is None:
                continue
            started = self.running.pop(found)
            started.job.close()
            self._learn(started)
            record = _record(self.connection, found)
            if record is not None and record.get("state") == "running":
                if started.killed:
                    text = started.killed
                elif exit_.signal is not None:
                    text = f"the job process ended on {signal.Signals(exit_.signal).name}"
                elif exit_.code:
                    text = f"the job process exited with code {exit_.code} and said nothing"
                else:
                    text = "the job process ended and nothing said how"
                finish(self.connection, found, "failed", error=text)
                self._act(f"settled light job {found}", text)
            else:
                self.connection.srem(running_key(), found)
            done.append(found)
        return done

    def _learn(self, started):
        """What the job's processes held at most is what the next claim of its kind starts from (`memorygate.learn`), as a worker's was. How long it took is kept for the kind's time limit (OPS.251), unless the clock ended it."""
        if not started.killed:
            from planetgen.queue import jobprofiles
            jobprofiles.record(started.claim.kind, seconds=self.clock() - started.started, peak_mb=started.peak_mb or None)
        if started.peak_mb > 0:
            try:
                memorygate.learn(self.connection, started.claim.kind, started.peak_mb)
            except Exception as exc:  # noqa: BLE001 -- a figure not kept is the table's figure next time
                log.swallowed("queue/lightjobs", "_learn", exc)

    def limits(self):
        """Kills the job that ran past its time limit or that was asked to stop (`request_stop`), and notes the most the others hold. Returns the ids killed."""
        killed = []
        now = self.clock()
        for job_id, started in list(self.running.items()):
            started.peak_mb = max(started.peak_mb, self.held_mb(started.job.pid))
            record = _record(self.connection, job_id)
            if started.killed or record is None:
                continue
            if record.get("stop") == "1":
                started.killed = "the job was stopped"
            elif started.limit and now - started.started > started.limit:
                started.killed = f"the job ran past its time limit of {int(started.limit)} s and was ended"
            else:
                continue
            try:
                started.job.kill()
            except (ProcessLookupError, OSError):
                pass
            killed.append(job_id)
            self._act(f"ended light job {job_id}", started.killed)
        return killed

    def recover(self):
        """The jobs `running` that this duty did not start (a supervisor that ended took its memory with it): one whose process is gone with no end written (or that never got a process) is failed. Returns the ids."""
        now = self.clock()
        if self.recovered_at is not None and now - self.recovered_at < RECOVER_SECONDS:
            return []
        self.recovered_at = now
        settled = []
        for job_id in [_text(item) for item in self.connection.smembers(running_key())]:
            if job_id in self.running:
                continue
            record = _record(self.connection, job_id)
            if record is None or record.get("state") != "running":
                self.connection.srem(running_key(), job_id)
                continue
            pid = int(record["pid"]) if record.get("pid") else None
            if pid is not None and _alive(pid):
                continue
            if pid is None and now - float(record.get("started") or now) < NO_PID_SECONDS:
                continue
            text = "the supervisor that started the job ended and its process is gone"
            finish(self.connection, job_id, "failed", error=text)
            settled.append(job_id)
            self._act(f"settled light job {job_id}", text)
        return settled

    # -- what is waiting ---------------------------------------------------------------------------------------------------

    def drop(self):
        """Drops each pending job nobody is waiting for any more (its `wanted_until` has passed): it fails with `DROPPED_TEXT` and its `drop_keys` are deleted. Returns the ids."""
        now = self.clock()
        dropped = []
        for job_id in [_text(item) for item in self.connection.zrange(pending_key(), 0, -1)]:
            record = _record(self.connection, job_id)
            if record is None:
                self.connection.zrem(pending_key(), job_id)
                continue
            if not record.get("wanted_until") or float(record["wanted_until"]) > now:
                continue
            if not self.connection.zrem(pending_key(), job_id):
                continue            # canceled, or taken, between the reads
            keys = json.loads(record.get("drop_keys") or "[]")
            if keys:
                self.connection.delete(*keys)
            finish(self.connection, job_id, "failed", error=DROPPED_TEXT)
            dropped.append(job_id)
            self._act(f"dropped light job {job_id}", f"nobody asked after {record.get('function', 'it')} for {float(record['wanted']):g} s")
        return dropped

    def claim_of(self, record):
        """The claim of a pending job: the API function's figure (`claims.api_table`), raised to what the gate has learned of the kind (`memorygate.claim_mb`); the default figure for a function with none."""
        from planetgen.queue import claims
        reference = record.get("function", "")
        lane = record.get("lane") or LANE_INTERACTIVE
        try:
            mb, cpu = claims.api_table(reference)
            kind = f"api:{reference}"
        except claims.NoClaim:
            mb, cpu, kind = claims.DEFAULT_MB, 1.0, f"fn:{reference}"
        try:
            mb = memorygate.claim_mb(self.connection, kind, mb)
        except Exception as exc:  # noqa: BLE001 -- the table's figure is the claim
            log.swallowed("queue/lightjobs", "claim_of", exc)
        return memorygate.Claim(kind, float(mb), float(cpu), lane)

    def tickets(self):
        """`(tickets, records)` of the oldest pending jobs: the `scheduler.Ticket`s `dispatch` judges, and each job's record by ticket name."""
        from planetgen.queue import scheduler
        records, tickets = {}, []
        for job_id in [_text(item) for item in self.connection.zrange(pending_key(), 0, WINDOW - 1)]:
            record = _record(self.connection, job_id)
            if record is None or record.get("state") != "queued":
                self.connection.zrem(pending_key(), job_id)
                continue
            name = f"light-{job_id}"
            records[name] = (job_id, record)
            tickets.append(scheduler.Ticket(name, self.claim_of(record), float(record["created"]), 0, "light"))
        return tickets, records

    def holders(self):
        return [memorygate.Holder(f"light-{job_id}", started.claim, started.job.pid, "light", started.claim.mb, started.claim.mb) for job_id, started in self.running.items()]

    def start(self, job_id, record, claim):
        """Takes the job off the pending list, marks it running and forks it. Returns the `Started`, or `None` for a job another turn, or a cancel, took."""
        if not self.connection.zrem(pending_key(), job_id):
            return None
        now = self.clock()
        pipe = self.connection.pipeline(transaction=True)
        pipe.hset(job_key(job_id), mapping={"state": "running", "started": f"{now:.3f}"})
        pipe.hdel(job_key(job_id), "waited")
        pipe.sadd(running_key(), job_id)
        pipe.execute()
        job = {"id": job_id, "function": record["function"], "args": json.loads(record["args"]), "redis_url": record.get("url") or None, "nice": PROCESS_NICE.get(claim.lane, 0), **self.extra}
        try:
            key, forked = self._forker().fork(KIND, self.runner, job, owner=OWNER)
        except jobhost.NotReady as exc:                 # its template is still being built: handed out at a later turn
            self._back(job_id, record, str(exc))
            return None
        except Exception as exc:  # noqa: BLE001 -- this job fails, the others go on
            log.swallowed("queue/lightjobs", "start", exc)
            finish(self.connection, job_id, "failed", error=f"the job could not be started: {type(exc).__name__}: {exc}"[:500])
            return None
        self.connection.hset(job_key(job_id), "pid", forked.pid)
        started = Started(job_id, forked, key, claim, now)
        started.limit = self._limit_of(claim.kind, float(record["limit"]) if record.get("limit") else None)
        self.running[job_id] = started
        return started

    def _limit_of(self, kind, typed):
        """The job's time limit: the kind's own (`learned_limit` of what it took before), else the typed one the job was added with. What was read of a kind is kept `LEARNED_KEEP_SECONDS`; a database that does not answer leaves the typed limit."""
        found = self.learned.get(kind)
        now = self.clock()
        if found is None or now - found[0] > LEARNED_KEEP_SECONDS:
            try:
                from planetgen.queue import jobprofiles
                found = (now, learned_limit(jobprofiles.seconds_of(kind)))
            except Exception as exc:  # noqa: BLE001 -- the typed limit holds
                log.swallowed("queue/lightjobs", "_limit_of", exc)
                found = (now, None)
            self.learned[kind] = found
        return found[1] if found[1] is not None else typed

    def _back(self, job_id, record, reason):
        """A job that could not be forked yet goes back to the pending list with its first place and a reason."""
        pipe = self.connection.pipeline(transaction=True)
        pipe.hset(job_key(job_id), mapping={"state": "queued", "waited": reason[:200]})
        pipe.hdel(job_key(job_id), "started")
        pipe.srem(running_key(), job_id)
        pipe.zadd(pending_key(), {job_id: float(record["created"])})
        pipe.execute()

    def _note(self, records, decided):
        for name, reason in decided.reasons.items():
            found = records.get(name)
            if found is None:
                continue
            job_id, record = found
            if record.get("waited") != reason[:200]:
                self.connection.hset(job_key(job_id), "waited", reason[:200])

    def turn(self):
        """One turn: ends, limits, drops, tickets, dispatch, starts. Returns the ids of the jobs started."""
        self.ended()
        self.limits()
        self.recover()
        self.drop()
        tickets, records = self.tickets()
        held = self.hold()
        if held and tickets:
            waiting = [ticket for ticket in tickets if ticket.claim.lane not in memorygate.SHORT_LANES]
            self._note(records, types.SimpleNamespace(reasons={ticket.name: held for ticket in waiting}))
            tickets = [ticket for ticket in tickets if ticket.claim.lane in memorygate.SHORT_LANES]
        if not tickets:
            return []
        decided = dispatch.dispatch(tickets, self.read_host(), [*self.holders(), *self.also_held()], self.clock(), {})
        self._note(records, decided)
        started = []
        for name in decided.admitted:
            job_id, record = records[name]
            made = self.start(job_id, record, next(ticket.claim for ticket in tickets if ticket.name == name))
            if made is not None:
                started.append(job_id)
        return started

    def busy(self):
        """Whether there is anything for the supervisor to watch: a job waiting or running."""
        return bool(self.connection.zcard(pending_key())) or bool(self.running) or bool(self.connection.scard(running_key()))

    def clear(self):
        """The supervisor ends: every job's pipe closes (a light job finishes what it is doing and exits when it is done), and the templates stop."""
        for started in self.running.values():
            started.job.close()
        self.running = {}
        stop = getattr(self.forker, "stop", None)
        if stop is not None:
            stop()


def _alive(pid):
    """Whether process `pid` exists."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True
