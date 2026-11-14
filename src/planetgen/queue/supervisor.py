# planetgen/queue/supervisor.py

"""
The supervisor (OPS.207, S1 of the job engine's stage 2): one process per host that, in the stages that follow, makes every decision about the queue
(`docs/design/operations/job-engine.md`, `/mnt/project-files/reports/job-engine-design-2026-11-13.md` section 3.2). Here it only lives: it takes the host's
lock, publishes that it is alive and since when, logs its turns, and exits after an idle spell. Nothing it does changes a run.

* One per host, by a kernel file lock (`flock`) on `supervisor.lock` in the state folder (`maintenance.state_dir`). The kernel frees the lock the moment
  the process ends, however it ends, so a dead supervisor never blocks the next.
* Started on demand and detached (`ensure_started`) by any submit (`workers.ensure_workers`), any status poll (the Queue page) and the maintenance
  warm-up's start; a sign of life in Redis (`alive_key`) says it need not be, and a short claim (`starting_key`) keeps a burst of callers from starting two.
* Its loop takes everything it touches as arguments (`Supervisor`: a clock, a store, a process table, the question "is there anything to watch?"), so the
  engine simulator (`/mnt/project-files/engine/stage0/enginesim.py`, OPS.212) can drive the real loop with a fake clock and fake store.
"""

import fcntl
import json
import os
import sys
import time

from planetgen.admin import maintenance
from planetgen.queue import redisqueue, runline, runrecord
from planetgen.util import log

LOCK_NAME = "supervisor.lock"
"""str: The file in the state folder the supervisor holds a lock on while it lives."""

TURN_SECONDS = 2.0
"""float: How long the supervisor waits between its turns."""

BUSY_TURN_SECONDS = 1.0
"""float: The wait between turns while a run is live, so a dead or quiet unit is judged as often as the starter used to (OPS.208)."""

IDLE_SECONDS = 300.0
"""float: How long the supervisor goes on with nothing to watch before it exits; the next submit, poll or timer starts it again."""

ALIVE_SECONDS = 15
"""int: How long the sign of life lasts after a turn writes it: several turns, so one slow turn is not a death."""

STARTING_SECONDS = 20
"""int: How long a caller's claim to be starting the supervisor holds off the other callers."""

RECENT_KEPT = 20
"""int: How many of its latest acts the supervisor keeps for the Queue page (OPS.214)."""

RECENT_SECONDS = 7 * 86400
"""int: How long the list of its latest acts is kept after the last one."""

LOG_EVERY_TURNS = 30
"""int: The turns between the supervisor's own lines in the log (a turn that changes nothing is not worth a line each)."""


def lock_path():
    """The file the supervisor's lock is taken on."""
    return os.path.join(maintenance.state_dir(), LOCK_NAME)


def take_lock(path=None):
    """
    Takes the host's supervisor lock without waiting.

    Returns:
        file or None: The open file that holds the lock (keep it open for as long as the lock is wanted), or `None` when another process holds it.
    """
    path = path or lock_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    handle = open(path, "a+", encoding="ascii")
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        handle.close()
        return None
    return handle


def alive_key():
    """The hash the supervisor publishes its sign of life in: `pid`, `since` (Unix seconds), `turn_at`, `turns`; it expires `ALIVE_SECONDS` after a turn."""
    return runrecord._k("supervisor", "alive")


def recent_key():
    """The list of the supervisor's latest acts, newest first, each a JSON `{"at", "what", "reason"}` (OPS.214): what it did and why, in words."""
    return runrecord._k("supervisor", "recent")


def nudges_key():
    """The hash of the asks for workers (OPS.209): queue name -> JSON `{"url", "mysql_host", "count"}` of the latest ask; the supervisor takes them all at each turn."""
    return runrecord._k("supervisor", "nudges")


def starting_key():
    """The string a caller sets to say it is starting the supervisor, so the callers that follow do not."""
    return runrecord._k("supervisor", "starting")


class RedisStore:
    """What the supervisor reads and writes in Redis: its sign of life."""

    def __init__(self, connection):
        self.connection = connection

    def publish(self, pid, since, now, turns):
        """Writes the sign of life: who, since when, the time of this turn and the count of turns."""
        pipe = self.connection.pipeline()
        pipe.hset(alive_key(), mapping={"pid": pid, "since": since, "turn_at": now, "turns": turns})
        pipe.expire(alive_key(), ALIVE_SECONDS)
        pipe.execute()

    def retract(self):
        """Takes the sign of life down, as the supervisor does when it ends."""
        self.connection.delete(alive_key())

    def note(self, now, what, reason=""):
        """Adds an act to the list the Queue page shows (OPS.214): what the supervisor did and why."""
        pipe = self.connection.pipeline()
        pipe.lpush(recent_key(), json.dumps({"at": now, "what": what, "reason": reason}))
        pipe.ltrim(recent_key(), 0, RECENT_KEPT - 1)
        pipe.expire(recent_key(), RECENT_SECONDS)
        pipe.execute()

    def busy(self):
        """Whether there is anything for a supervisor to watch: a live run, a run waiting for the machine, or an ask for workers not acted on yet."""
        return bool(self.connection.hlen(nudges_key())) or bool(self.connection.hlen(runline._starts_key())) or any(runrecord.active_run(self.connection, kind) for kind in (runrecord.GENERATION, runrecord.MAINTENANCE, runrecord.CHECK)) \
            or bool(self.connection.zcard(runline._line_key()))


class SystemClock:
    """The wall clock; the simulator gives a fake with the same two calls."""

    def now(self):
        return time.time()

    def sleep(self, seconds):
        time.sleep(seconds)


class ProcessTable:
    """The host's processes as the supervisor asks about them; the simulator gives a fake."""

    def alive(self, pid):
        """Whether process `pid` exists."""
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        return True


RUN_KINDS = (runrecord.GENERATION, runrecord.MAINTENANCE, runrecord.CHECK)
"""tuple: The kinds of run that can be live; each holds its own active key (`runrecord.active_run`)."""


def running_stage(conn, run_id):
    """
    The stage of the run whose units are being worked now (its plan is stored and not every unit is done), or `None`: no stage is in flight (the run is between stages, or its next stage has not
    been cut yet). At most one stage of a run is in flight.
    """
    count = len(json.loads(runrecord._s(conn.hget(runrecord.run_key(run_id), "stages")) or "[]"))
    for index in range(count):
        plan = runrecord.read_plan(conn, run_id, index)
        if plan is None:
            return None
        if plan["units"] and conn.scard(runrecord.done_key(run_id, index)) < plan["units"]:
            return index
    return None


def live_runs(connection):
    """The ids of the runs that are live now, one at most for each kind."""
    found = []
    for kind in RUN_KINDS:
        run_id = runrecord.active_run(connection, kind)
        if run_id and runrecord._s(connection.hget(runrecord.run_key(run_id), "state")) in runrecord.LIVE_STATES:
            found.append(run_id)
    return found


class RunWatches:
    """
    OPS.208: the watch of every live run's units, which was the starter's (`unitwatch.UnitWatch`, `unitstuck.Stuck`, OPS.138 and OPS.140) and is the supervisor's now. At each `turn` it finds
    the stage in flight of each live run, judges the beat of its running units (a change of state logged once, the dead and the quiet acted on after two checks), and drops the watch of a
    stage that ended or of a run that is over. A run whose starting process is gone is watched all the same (W2 of the audit).

    Args:
        connection (redis.Redis): The server the runs are on.
        make (callable, optional): `make(connection, run_id, stage, kind)` builds the watch of one stage (`unitwatch.UnitWatch` with the run's gap history and rates where left out); the
            simulator passes its own.
    """

    def __init__(self, connection, make=None):
        self.connection = connection
        self.make = make or self._unit_watch
        self.watches = {}

    @staticmethod
    def _unit_watch(connection, run_id, stage, kind):
        from planetgen.queue import commandgaps, unitrates, unitwatch
        return unitwatch.UnitWatch(connection, run_id, stage, kind, history_of=commandgaps.history_of(connection, run_id), rate_of=unitrates.rate_of(connection, run_id))

    def _kind_of(self, run_id, stage):
        """The command of the stage, from the steps the run recorded."""
        key = json.loads(runrecord._s(self.connection.hget(runrecord.run_key(run_id), "stages")))[stage]["key"]
        return next(item["kind"] for item in runrecord.read_steps(self.connection, run_id) if item["key"] == key)

    def live_runs(self):
        """The ids of the runs that are live now, one at most for each kind."""
        return live_runs(self.connection)

    def turn(self):
        """
        Judges the units of every stage in flight, once.

        Returns:
            list: The `(run_id, stage)` of each stage watched this turn.
        """
        wanted = {}
        for run_id in self.live_runs():
            stage = running_stage(self.connection, run_id)
            if stage is not None:
                wanted[(run_id, stage)] = None
        for key in [key for key in self.watches if key not in wanted]:
            self._drop(key)
        for key in wanted:
            try:
                if key not in self.watches:
                    self.watches[key] = self.make(self.connection, key[0], key[1], self._kind_of(*key))
                self.watches[key].check()
            except Exception as exc:  # noqa: BLE001 -- one run's watch failing must not stop the watch of the others, or the supervisor
                log.swallowed("queue/supervisor", "RunWatches.turn", exc)
        return sorted(wanted)

    def _drop(self, key):
        watch = self.watches.pop(key)
        try:
            watch.clear()
        except Exception as exc:  # noqa: BLE001 -- the verdicts are for the page
            log.swallowed("queue/supervisor", "RunWatches._drop", exc)

    def clear(self):
        """Drops every watch (the supervisor ends)."""
        for key in list(self.watches):
            self._drop(key)


SUPERVISED_SECONDS = 15
"""int: How long a run counts as watched by the supervisor after a turn that saw it (`runrecord.supervised_key`): several turns, so one slow turn is not a death."""

REPAIR_GRACE_SECONDS = 30.0
"""float: How long a stage may stand finished without its end having queued the next opening before the supervisor ends it (OPS.186): the unit that ends a stage is normally at it at once, and
a conclusion can take a while, so a lease (`runrecord.ending_key`) still held is waited for however long it is."""


class Reconciler:
    """
    OPS.210: each turn, and first thing when the supervisor starts, compares what the run records say with what is so and corrects it, and says so in the supervisor's acts (OPS.214).

    * Every live run is marked as watched (`runrecord.supervised_key`), so a new run does not take over a run whose starting process is gone (`runrecord.holder_is_dead`): the run
      goes on without it, its stages advancing from the units, and ends by itself.
    * A stage that finished (its last unit counted it) and whose end never queued the next opening (that unit's process was killed, or Redis failed it), found after `REPAIR_GRACE_SECONDS`, is ended
      by `units.repair_unended_stages` (OPS.186): the conclusion not yet done is done, the opening not yet queued is queued, each once.

    Args:
        connection (redis.Redis): The server the runs are on.
        clock (callable, optional): `clock() -> Unix seconds`; the wall clock where left out.
        repair (callable, optional): `repair(connection, run_id) -> [stages]`; `units.repair_unended_stages` where left out. The simulator passes its own.
    """

    def __init__(self, connection, clock=time.time, repair=None):
        self.connection = connection
        self.clock = clock
        self.repair = repair
        self.supervisor = None

    def attach(self, supervisor):
        """The supervisor that runs this duty, so what it corrects goes into the supervisor's acts."""
        self.supervisor = supervisor

    def _act(self, what, reason):
        if self.supervisor is not None:
            self.supervisor.act(what, reason)

    def turn(self):
        """
        Marks the live runs watched and ends the stages that stand finished without their end.

        Returns:
            list: The `(run_id, stage)` of each stage it ended.
        """
        from planetgen.queue import units
        repair = self.repair or units.repair_unended_stages
        repaired = []
        for run_id in live_runs(self.connection):
            self.connection.set(runrecord.supervised_key(run_id), "1", ex=SUPERVISED_SECONDS)
            for stage in units.unended_stages(self.connection, run_id):
                finished = runrecord._s(self.connection.hget(runrecord.stage_key(run_id, stage), "ended"))
                if finished and self.clock() - float(finished) < REPAIR_GRACE_SECONDS:
                    continue
                if self.connection.exists(runrecord.ending_key(run_id, stage)):
                    continue            # a process is ending it now
                try:
                    done = repair(self.connection, run_id)
                except Exception as exc:  # noqa: BLE001 -- a stage that could not be ended now is ended at the next turn, or fails the run through the unit's own failure
                    log.swallowed("queue/supervisor", "Reconciler.turn", exc)
                    continue
                repaired.extend((run_id, index) for index in done)
                if done:
                    self._act(f"ended stage {stage + 1} of run {run_id}", "it had finished and its end had not queued the next stage")
        return repaired

    def clear(self):
        """Nothing to drop."""


class WorkerStarts:
    """
    OPS.209: the starting of workers, which every enqueuer, web thread and starter did for itself (behind a host-wide lock) and is the supervisor's alone now. Each turn it takes the asks for workers
    (`nudge`) and acts on each with `workers.start_workers`, which holds the rules: the caps, the lanes, the memory limit, shedding and lending. Two asks for one queue before a turn are one:
    the later wins.

    Args:
        connection (redis.Redis): The server the asks are on.
        start (callable, optional): `start(queue_name, connection, url, mysql_host, count)`; `workers.start_workers` where left out. The simulator passes its own.
    """

    def __init__(self, connection, start=None):
        self.connection = connection
        self.start = start

    def take(self):
        """The asks waiting, taken and cleared in one step: a list of `(queue_name, ask)`."""
        pipe = self.connection.pipeline(transaction=True)
        pipe.hgetall(nudges_key())
        pipe.delete(nudges_key())
        asks, _ = pipe.execute()
        return [(runrecord._s(name), json.loads(runrecord._s(raw))) for name, raw in asks.items()]

    def turn(self):
        """
        Acts on every ask for workers.

        Returns:
            int: How many workers it started.
        """
        from planetgen.queue import workers
        start = self.start or workers.start_workers
        started = 0
        for queue_name, ask in self.take():
            try:
                started += start(queue_name, self.connection, ask.get("url"), ask.get("mysql_host"), ask.get("count")) or 0
            except Exception as exc:  # noqa: BLE001 -- one queue's start failing leaves the others, and its next ask
                log.swallowed("queue/supervisor", "WorkerStarts.turn", exc)
        return started

    def clear(self):
        """Nothing to drop: an ask not acted on stays on Redis for the next supervisor."""


class RunStarts:
    """
    OPS.211: the starting of runs, which each process that wanted one did for itself and is the supervisor's alone now. A process that wants a run leaves its ask (`runline.start_when_free`: a
    place in the waiting line for a run that waits, a request of its own for one that does not) and follows the run the supervisor started for it. Each turn this answers the asks that do not wait
    (a run that finds the machine taken is refused, with the run that holds it named) and starts the first run of the line when the machine is free; a run that finds it taken stays in the line.

    Args:
        connection (redis.Redis): The server the asks are on.
        start (callable, optional): `start(connection, **arguments)`; `runrecord.start_run` where left out. The simulator passes its own.
    """

    def __init__(self, connection, start=None):
        self.connection = connection
        self.start = start
        self.supervisor = None

    def attach(self, supervisor):
        self.supervisor = supervisor

    def _act(self, what, reason):
        if self.supervisor is not None:
            self.supervisor.act(what, reason)

    def _begin(self, spec):
        """Starts the run an ask describes; raises `runrecord.Busy` when the machine is taken."""
        arguments = dict(spec)
        arguments["stages"] = [runrecord.Stage(**stage) for stage in arguments["stages"]]
        return (self.start or runrecord.start_run)(self.connection, **arguments)

    def turn(self):
        """
        Answers every ask that does not wait and starts the first waiting run if the machine is free.

        Returns:
            int: How many runs it started.
        """
        started = 0
        for ticket, spec in runline.direct_starts(self.connection).items():
            if not runline.take_direct(self.connection, ticket):
                continue
            try:
                run_id = self._begin(spec)
            except runrecord.Busy as busy:
                runline.decide(self.connection, ticket, busy=str(busy))
            except Exception as exc:  # noqa: BLE001 -- the ask is answered with the failure, and the others are still answered
                log.swallowed("queue/supervisor", "RunStarts.turn", exc)
                runline.decide(self.connection, ticket, failed=f"the run could not be started: {exc}")
            else:
                runline.decide(self.connection, ticket, run_id=run_id)
                self._act("started a run", f"{spec.get('title') or 'a run'} ({run_id}), asked for without waiting")
                started += 1
        waiting = runline.waiting_start(self.connection)
        if waiting is not None:        # the supervisor of the host is the only one that starts a run, so the ask is not claimed
            ticket, spec = waiting
            try:
                run_id = self._begin(spec)
            except runrecord.Busy:
                pass        # it keeps its place; the machine is taken
            except Exception as exc:  # noqa: BLE001 -- the run that asked is told, and the line goes on
                log.swallowed("queue/supervisor", "RunStarts.turn", exc)
                runline.decide(self.connection, ticket, failed=f"the run could not be started: {exc}")
                runline.leave(self.connection, ticket)
            else:
                runline.decide(self.connection, ticket, run_id=run_id)
                runline.leave(self.connection, ticket)
                self._act("started a run", f"{spec.get('title') or 'a run'} ({run_id}), first in the line and the machine free")
                started += 1
        return started

    def clear(self):
        """Nothing to drop: an ask not answered stays on Redis for the next supervisor."""


class Supervisor:
    """
    The supervisor's loop with its dependencies passed in (OPS.207): `clock` (`now()`, `sleep(seconds)`), `store` (`publish(pid, since, now, turns)`,
    `retract()`, `busy()`), `table` (`alive(pid)`), the process id it publishes, and `stop`, a callable that says the process was asked to end.
    `lines` is called with each line the supervisor logs, in place of the log, for a test.
    """

    def __init__(self, clock, store, table, pid=None, stop=lambda: False, idle_seconds=IDLE_SECONDS, turn_seconds=TURN_SECONDS, lines=None, duties=()):
        self.clock = clock
        self.store = store
        self.table = table
        self.pid = os.getpid() if pid is None else pid
        self.stop = stop
        self.idle_seconds = idle_seconds
        self.turn_seconds = turn_seconds
        self.lines = lines
        self.duties = tuple(duties)      # what it does each turn (OPS.208 to OPS.210): objects with `turn()`, `clear()` and, if they say what they did, `attach(supervisor)`
        for duty in self.duties:
            attach = getattr(duty, "attach", None)
            if attach is not None:
                attach(self)
        self.since = None
        self.turns = 0
        self.idle_since = None

    def say(self, text):
        if self.lines is not None:
            self.lines(text)
        else:
            log.debug(text, area="queue")

    def act(self, what, reason=""):
        """Says what the supervisor did and why: in the log and, where the store keeps them, in the list the Queue page shows (OPS.214). Later items call it for each thing they do."""
        self.say(f"The supervisor {self.pid}: {what}" + (f" ({reason})" if reason else "") + ".")
        note = getattr(self.store, "note", None)
        if note is not None:
            note(self.clock.now(), what, reason)

    def turn(self):
        """One turn: publishes the sign of life, then notes whether there is anything to watch. Nothing else yet (OPS.208 and the items after it add the watch).

        Returns:
            bool: Whether there is anything to watch.
        """
        now = self.clock.now()
        if self.since is None:
            self.since = now
            self.act("started", "something asked for one")
        self.turns += 1
        self.store.publish(self.pid, self.since, now, self.turns)
        busy = bool(self.store.busy())
        for duty in self.duties:
            try:
                duty.turn()
            except Exception as exc:  # noqa: BLE001 -- a turn that failed is the next turn's to do again, and one duty failing leaves the others
                log.swallowed("queue/supervisor", "Supervisor.turn", exc)
        if self.turns == 1 or self.turns % LOG_EVERY_TURNS == 0:
            self.say(f"The supervisor {self.pid}, turn {self.turns}: {'something to watch' if busy else 'nothing to watch'}.")
        return busy

    def run(self, max_turns=None):
        """
        Runs turns until it is asked to stop, has been idle for `idle_seconds`, or has run `max_turns`.

        Returns:
            str: Why it ended: `stopped`, `idle` or `turns`.
        """
        while True:
            busy = self.turn()
            now = self.clock.now()
            if busy:
                self.idle_since = None
            elif self.idle_since is None:
                self.idle_since = now
            if self.stop():
                reason = "stopped"
            elif self.idle_since is not None and now - self.idle_since >= self.idle_seconds:
                reason = "idle"
            elif max_turns is not None and self.turns >= max_turns:
                reason = "turns"
            else:
                self.clock.sleep(BUSY_TURN_SECONDS if busy else self.turn_seconds)
                continue
            for duty in self.duties:
                duty.clear()
            self.store.retract()
            self.act("ended", {"stopped": "it was asked to stop", "idle": f"nothing to watch for {self.idle_seconds:.0f} s",
                               "turns": "the turn limit"}[reason] + f", after {self.turns} turn(s)")
            return reason


def overview(connection):
    """
    What the Queue page and `GET /api/admin/work/supervisor` show of the supervisor (OPS.214): `alive`, and while it is, its `pid`, `since` and `turn_at` (Unix seconds) and `turns`; and
    `recent`, its latest acts newest first, each `at`, `what`, `reason`.
    """
    fields = {runrecord._s(k): runrecord._s(v) for k, v in connection.hgetall(alive_key()).items()}
    recent = []
    for raw in connection.lrange(recent_key(), 0, RECENT_KEPT - 1):
        try:
            recent.append(json.loads(runrecord._s(raw)))
        except ValueError:
            continue
    body = {"alive": bool(fields), "recent": recent}
    if fields:
        body.update(pid=int(fields["pid"]), since=float(fields["since"]), turn_at=float(fields["turn_at"]), turns=int(fields["turns"]))
    return body


def nudge(connection, queue_name, url=None, mysql_host=None, count=None):
    """
    Asks the supervisor for workers for the queue `queue_name` (OPS.209: `workers.ensure_workers`) and starts the supervisor if none is alive. Never raises into the caller: an ask that was not
    left is the next poll's to leave again.
    """
    try:
        connection.hset(nudges_key(), queue_name, json.dumps({"url": url, "mysql_host": mysql_host, "count": count}))
        connection.expire(nudges_key(), RECENT_SECONDS)
    except Exception as exc:  # noqa: BLE001 -- the engine's own bookkeeping never fails the submit that asked
        log.swallowed("queue/supervisor", "nudge", exc)
        return
    ensure_started(connection, url)


def is_alive(connection):
    """Whether a supervisor's sign of life is up."""
    return bool(connection.exists(alive_key()))


def supervisor_argv(url=None, python=None):
    """The command line that starts a supervisor."""
    return [python or sys.executable, "-m", "planetgen.cli.supervisor", *(["--url", url] if url else [])]


def _spawn(argv):
    """Starts the supervisor process, detached, in the repository's folder (the test suite replaces this: a test never starts a process of the engine)."""
    from planetgen.queue import workers
    redisqueue.start_detached(argv, workers.REPO_DIR)


def ensure_started(connection=None, url=None, start=None):
    """
    Starts a supervisor, detached, unless one is alive or another caller is starting one now. Cheap for the caller that finds one alive (one Redis
    read), so a submit, a status poll and a timer each call it. Never raises: a supervisor that could not be started is the next call's to start.

    Args:
        start (callable, optional): What starts the process, given the command line; default `redisqueue.start_detached` in the repository's folder.

    Returns:
        bool: Whether this call started one.
    """
    try:
        connection = connection or redisqueue.connect(url)
        if is_alive(connection) or not connection.set(starting_key(), os.getpid(), nx=True, ex=STARTING_SECONDS):
            return False
        argv = supervisor_argv(url)
        (start or _spawn)(argv)
        return True
    except Exception as exc:  # noqa: BLE001 -- the engine's own bookkeeping never fails the submit or the page that asked
        log.swallowed("queue/supervisor", "ensure_started", exc)
        return False
