"""
OPS.212: the supervisor under a fake clock with faults, in the manner of the stage 0 skeleton (`/mnt/project-files/engine/stage0/enginesim.py`: a seed draws every event, a failing seed
replays exactly, nothing sleeps). Where the skeleton drives the three pure rules for memory, this drives the supervisor's duties for runs (`RunStarts`, `Reconciler`) and its loop
(`Supervisor`) on the suite's real Redis with a clock of its own.

What a seed draws, each turn:

* asks for runs: one that waits (a place in the line) or one that does not (checks and maintenance), at a random moment;
* the work of the run that holds the machine: its unit finishes and the end of the stage either queues the next stage or its process is killed in it (`stage_unended` is left);
* the starter of a run stops (its `alive` key is gone);
* the supervisor is killed and a new one takes its place, at any turn (between any two duties);
* Redis refuses a command (`FlakyRedis`: a `ConnectionError` before the command is applied).

The invariants (`Sim.check`):

* every ask is answered once at most, and each answered ask has one run: no run was started twice for it, and no ask started a run it was not answered with;
* a run waiting in the line is started in the order it came, one at a time, and a run that does not wait is never left in the line;
* a finished stage that stood unended past the grace is ended exactly once, and a stage that is ended is never ended again;
* a run whose starter stopped is said so once;
* when the faults stop, every ask has been answered and every stage ended: nothing waits for ever.
"""

import json
import random

import redis

from planetgen.queue import planrows, runline, runrecord, supervisor, units

TURN = 2.0
STAGES = 3


class Clock:
    """The fake clock of the supervisor and its duties."""

    def __init__(self):
        self.t = 1_000_000.0

    def now(self):
        return self.t

    def time(self):
        return self.t

    def sleep(self, seconds):
        self.t += seconds


class FlakyRedis:
    """A client that refuses a command now and then (a `ConnectionError`, before the command is applied), while `chance` is above 0."""

    def __init__(self, real, rng):
        self._real, self._rng, self.chance = real, rng, 0.0
        self.refused = 0

    def __getattr__(self, name):
        attribute = getattr(self._real, name)
        if not callable(attribute) or name in ("pipeline", "register_script"):
            return attribute

        def call(*args, **kwargs):
            if self.chance and self._rng.random() < self.chance:
                self.refused += 1
                raise redis.exceptions.ConnectionError("simulated refusal")
            return attribute(*args, **kwargs)

        return call


class Sim:
    """One seed's run of the supervisor's world."""

    def __init__(self, real, seed, turns=120, faults=True):
        self.rng = random.Random(seed)
        self.real, self.turns, self.faults = real, turns, faults
        self.conn = FlakyRedis(real, self.rng)
        self.clock = Clock()
        self.lines = []
        self.begun = set()
        self.ends = []                      # (run_id, stage) in the order stages were ended
        self.asks = {}                      # ticket -> {"waits": bool, "order": n}
        self.arrivals = 0
        self.started_for = {}               # ticket -> run ids it was answered with
        self.supervisor = None
        self.killed = 0
        self.new_supervisor()

    def repair(self, connection, run_id):
        """What the reconciler calls to end the stages left unended: ends them the way the unit that completes a stage does, and counts it."""
        found = units.unended_stages(self.real, run_id)
        for stage in found:
            self.end_stage(run_id, stage)
        return found

    def end_stage(self, run_id, stage):
        """The end of a stage: the next stage's opening queued once (its `release` job done in the ledger) and the stage begun; the last stage ends the run."""
        self.ends.append((run_id, stage))
        planrows.mark_ended(self.real, run_id, stage, planrows.RELEASE)
        if stage + 1 < STAGES:
            self.begun.add((run_id, stage + 1))
            runrecord.begin_stage(self.real, run_id, stage + 1, now=self.clock.t)
        else:
            runrecord.finish_run(self.real, run_id, "done", now=self.clock.t)

    def new_supervisor(self):
        """A supervisor in place of the one that was killed (its duties are fresh: nothing but what is on Redis survives)."""
        duties = (supervisor.Reconciler(self.conn, clock=self.clock.time, repair=self.repair), supervisor.RunStarts(self.conn))
        self.supervisor = supervisor.Supervisor(self.clock, supervisor.RedisStore(self.conn), supervisor.ProcessTable(), pid=self.rng.randrange(10, 10_000),
                                                lines=self.lines.append, duties=duties)

    def ask(self):
        """A process asks for a run: in the line, or at once (the ask is left the way `runline.start_when_free` leaves it)."""
        waits = self.rng.random() < 0.7
        start = {"title": f"Run {self.arrivals}", "stages": [{"key": f"s{n}", "label": f"S{n}"} for n in range(STAGES)]}
        try:
            if waits:
                ticket = runline.join(self.real, runline.RUN, start["title"], fields={"start": json.dumps(start)}, force=True)
            else:
                ticket = f"direct{self.arrivals}"
                self.real.hset(runline._starts_key(), ticket, json.dumps(start))
        except redis.exceptions.RedisError:
            return
        self.asks[ticket] = {"waits": waits, "order": self.arrivals}
        self.arrivals += 1

    def work(self):
        """The machine's run: its unit finishes, and the process that finished it either ends the stage or is killed in it."""
        run_id = runrecord.active_run(self.real)
        if run_id is None or self.rng.random() > 0.6:
            return
        stage = int(float(runrecord._s(self.real.hget(runrecord.run_key(run_id), "current")) or 0))
        if (run_id, stage) not in self.begun:
            self.begun.add((run_id, stage))
            runrecord.begin_stage(self.real, run_id, stage, now=self.clock.t)       # the follower opens the stage
            return
        if self.real.exists(runrecord.advanced_key(run_id, stage)):
            return                          # finished already; its end is the unit's, or the reconciler's if the unit was killed
        if runrecord.unit_done(self.real, run_id, stage, 0, now=self.clock.t).advance:
            if self.faults and self.rng.random() < 0.3:
                return                      # killed in the end: the stage stands finished and unended
            self.end_stage(run_id, stage)

    def step(self):
        """One turn of the world, then one turn of the supervisor (which a fault may have replaced or refused)."""
        if self.rng.random() < 0.12:
            self.ask()
        self.work()
        run_id = runrecord.active_run(self.real)
        if self.faults and run_id and self.rng.random() < 0.05:
            self.real.delete(runrecord.alive_key(run_id))       # its starter stops
        elif run_id and self.rng.random() < 0.9:
            runrecord.touch(self.real, run_id)                  # its starter is following it
        if self.faults and self.rng.random() < 0.05:
            self.killed += 1
            self.new_supervisor()
        self.conn.chance = 0.08 if self.faults else 0.0
        try:
            self.supervisor.turn()
        except redis.exceptions.RedisError:
            pass                                                # the supervisor's loop survives a turn that failed (`Supervisor.run`)
        self.conn.chance = 0.0
        self.clock.sleep(TURN)
        self.collect()

    def collect(self):
        """Reads the answers: which run each ask was answered with."""
        for ticket in self.asks:
            found = runline.verdict(self.real, ticket)
            if found and found.get("run_id") and found["run_id"] not in self.started_for.setdefault(ticket, []):
                self.started_for[ticket].append(found["run_id"])

    def run(self):
        """The faults, and then quiet turns until every ask is answered and every run is over, or the turns run out."""
        for _ in range(self.turns):
            self.step()
        self.faults = False
        for _ in range(self.turns * 4):
            if not self.pending():
                break
            self.step()
        return self

    def pending(self):
        """Whether anything is waiting: an ask without an answer, a run that is live."""
        return any(runline.verdict(self.real, ticket) is None for ticket in self.asks) or runrecord.active_run(self.real) is not None \
            or bool(runline.of_type(self.real, runline.RUN))

    def check(self):
        """The invariants; raises `AssertionError` naming what broke."""
        started = [run_id for runs in self.started_for.values() for run_id in runs]
        assert len(started) == len(set(started)), f"a run answers two asks: {self.started_for}"
        assert all(len(runs) <= 1 for runs in self.started_for.values()), f"an ask was answered with two runs: {self.started_for}"
        assert len(self.ends) == len(set(self.ends)), f"a stage was ended twice: {sorted(self.ends)}"
        order = [self.asks[ticket]["order"] for ticket, runs in self.started_for.items() if runs and self.asks[ticket]["waits"]]
        started_order = sorted(((self.real.hget(runrecord.run_key(runs[0]), "created"), self.asks[ticket]["order"]) for ticket, runs in self.started_for.items()
                                if runs and self.asks[ticket]["waits"]), key=lambda pair: float(runrecord._s(pair[0])))
        assert [n for _created, n in started_order] == sorted(order), f"the line was not served in the order it came: {started_order}"
        for ticket, ask in self.asks.items():
            answer = runline.verdict(self.real, ticket)
            assert answer is not None, f"ask {ticket} was never answered ({self.clock.t}): the line {runline.of_type(self.real, runline.RUN)}"
            assert not ask["waits"] or ticket not in [entry["ticket"] for entry in runline.of_type(self.real, runline.RUN)], f"{ticket} is still in the line"
        for run_id in set(started):
            assert runrecord.run_state(self.real, run_id) == "done", f"run {run_id} did not end: {runrecord.run_state(self.real, run_id)}"
            events = [fields for _id, fields in runrecord.read_events(self.real, run_id) if fields.get("kind") == "run.starter_stopped"]
            assert len(events) <= 1, f"the starter's stop was said {len(events)} times for {run_id}"
