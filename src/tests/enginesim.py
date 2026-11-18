"""
The dispatch simulator (OPS.227; the stage 0 skeleton `reports/engine/stage0/enginesim.py` with the watch and the hung jobs left out, which are OPS.193 and OPS.194's).

It drives `memorygate.decide` and `scheduler.judge` through `dispatch.dispatch` with a fake clock and a fake 16 GB host: one stage of bulk units, a maintenance warm-up beside it, a visitor's jobs arriving
through it and one cut that needs the bulk lane to itself, with claims that are too small (the kernel kills the biggest process; the retry has a claim grown by `memorygate.GROW_FACTOR`) and kills
drawn from the seed. A run is a function of its seed, so a failing seed replays exactly. Nothing here touches Redis, MySQL or a process.
"""

import math
import random
from dataclasses import dataclass, field

from planetgen.queue import dispatch, memorygate, scheduler, tokens
from planetgen.queue.memorygate import Claim, Holder, Host

PLACES = 2
BASE_MB = 2500.0
TICK = 1.0
BIG_MB = 7000.0


@dataclass
class Job:
    name: str
    claim: Claim
    seconds: float
    peak_mb: float
    queued: float = 0.0
    kill_at: float = None
    state: str = "ready"
    started: float = None
    finished: float = None
    done_times: int = 0
    waits: list = field(default_factory=list)

    def held_mb(self, now):
        """Memory ramps to the peak over the first fifth of the attempt."""
        return self.peak_mb * min(1.0, (now - self.started) / max(1.0, 0.2 * self.seconds))


class DeciderCrash(Exception):
    """What a decider that fails in a turn raises (OPS.236): the supervisor's duty fails, the turn admits nothing and the next turn does it again (`supervisor.Supervisor.turn`)."""


def scenario(rng, units, cut_at=None, named=False, big=False, small=0.02):
    """One stage of bulk units, a warm-up beside it, a visitor's jobs through it and the cut that needs the lane to itself (and, with `big`, a claim of nearly half the host); `small` is the share of the units whose claim is too small."""
    jobs = []
    for i in range(units):
        mb = rng.choice((600, 900, 900, 1400))
        peak = mb * rng.uniform(0.55, 0.98)
        if rng.random() < small:
            peak = mb * rng.uniform(1.3, 1.7)                   # a claim that is too small
        job = Job(f"unit-{i}", Claim("scatter-mass", mb, 1.0, memorygate.LANE_BULK), rng.lognormvariate(4.0, 0.5), peak)
        if rng.random() < 0.01:
            job.kill_at = rng.uniform(1, job.seconds)
        jobs.append(job)
    for i in range(units // 4):
        jobs.append(Job(f"warm-{i}", Claim("warm-systems", 500, 1.0, memorygate.LANE_MAINTENANCE), rng.lognormvariate(3.0, 0.4), 400))
    for i in range(units // 10):
        job = Job(f"page-{i}", Claim("fn:page", 300, 1.0, memorygate.LANE_INTERACTIVE), rng.uniform(0.5, 4.0), 250)
        job.queued = rng.uniform(0, units * 4.0)
        jobs.append(job)
    if named:           # OPS.232: a tenth of the units hold the one writer's token, three warm-ups build one warm set and one restarted unit runs alone
        for job in jobs[: units // 10 * 10 : 10]:
            job.claim = Claim(job.claim.kind, job.claim.mb, job.claim.cpu, job.claim.lane, False, (tokens.MYSQL_WRITE,))
        for job in [job for job in jobs if job.name.startswith("warm-")][:3]:
            job.claim = Claim(job.claim.kind, job.claim.mb, job.claim.cpu, job.claim.lane, False, (tokens.warm_set("tables"),))
        alone = Job("restarted-unit", Claim("scatter-mass", 900, 1.0, memorygate.LANE_BULK, False, (tokens.ALONE,)), 30.0, 700)
        alone.queued = units * 0.5
        jobs.append(alone)
    if big:             # OPS.236: a claim no unit's neighbour can leave room for without the rule keeping it
        jobs.append(Job("big-claim", Claim("scatter-mass", BIG_MB, 2.0, memorygate.LANE_BULK), 60.0, BIG_MB * 0.8, queued=units * 0.5))
    cut = Job("cut-next-stage", Claim("open:scatter-luminosity", 3000, 1.0, memorygate.LANE_BULK, exclusive=True), 90.0, 2600)
    cut.queued = units * 1.0 if cut_at is None else cut_at
    jobs.append(cut)
    return jobs


def old_rule(waiting, host, held, now, used, window, places):
    """What main's order does with the whole backlog in one ranking, the reservation ending at the gate's bound (the control the findings were measured against)."""
    held, admitted = list(held), []
    ranked = scheduler.rank(waiting, now, used)[:window]
    pending = list(ranked)
    for ticket in ranked:
        others = [other for other in pending if other.name != ticket.name]
        if scheduler.judge(ticket, others, host, held, now, used, places=places).admit:
            admitted.append(ticket.name)
            pending = others
            held.append(Holder(ticket.name, ticket.claim, 0, "sim", 0.0, 0.0))
    return dispatch.Turn(admitted, {})


def run(seed, units=100, rule=None, window=dispatch.WINDOW, total_mb=16000.0, cores=8, horizon=100000.0, cut_at=None, named=False, big=False, crash=0.0, small=0.02, pressure=()):
    """
    One simulated run with `rule` (the dispatch rule by default). Returns `(measures, broken)`: what was measured and the invariants that failed.

    `crash` is the share of the turns whose decider fails (`DeciderCrash`, drawn from the seed apart from the scenario's draws): such a turn admits nothing, as the supervisor's does. The invariants
    `pressure` is `(from, to, MB)` windows in which something outside the engine (the database's buffer pool, a visitor's tool) holds that much of the host: the host says so, and the kernel kills the
    biggest process when the jobs and the rest come to more than the host. A job the kernel killed runs again alone on a claim grown by `memorygate.GROW_FACTOR` (OPS.236). The invariants
    checked in every turn are that the claims of the jobs running never come to more than the host less its reserve (but for one claim that is larger and runs alone), that a job's token is held by one job, that the job that runs alone runs alone,
    and that the kernel kills only when a job went past its claim; at the end that each job finished once.
    """
    rule = rule or (lambda waiting, host, held, now, used, window, places: dispatch.dispatch(waiting, host, held, now, used, window, places))
    rng = random.Random(seed)
    jobs = scenario(rng, units, cut_at, named, big, small)
    crashes = random.Random(f"crash-{seed}")
    by_name = {job.name: job for job in jobs}
    now, used, broken = 0.0, {}, []
    measures = {"oom_kills": 0, "fault_kills": 0, "turns": 0, "crashed_turns": 0}
    budget = total_mb - memorygate.reserve_mb(Host(total_mb, total_mb, cores, 0.0))
    while now < horizon and any(job.state != "done" for job in jobs):
        now += TICK
        for job in [job for job in jobs if job.state == "running"]:
            spent = now - job.started
            if job.kill_at is not None and spent >= job.kill_at:
                job.kill_at = None
                job.state, job.queued = "ready", now
                measures["fault_kills"] += 1
            elif spent >= job.seconds:
                job.state = "done"
                job.finished = now
                job.done_times += 1
                used[job.claim.kind] = used.get(job.claim.kind, 0.0) + job.seconds * job.claim.cpu
        running = [job for job in jobs if job.state == "running"]
        outside = sum(mb for start, end, mb in pressure if start <= now < end)
        while running and BASE_MB + outside + sum(job.held_mb(now) for job in running) > total_mb:
            victim = max(running, key=lambda job: job.held_mb(now))
            victim.state, victim.queued = "ready", now
            victim.claim = Claim(victim.claim.kind, max(victim.claim.mb, victim.peak_mb) * memorygate.GROW_FACTOR, victim.claim.cpu, victim.claim.lane, victim.claim.exclusive,
                                 tuple(dict.fromkeys((*victim.claim.tokens, tokens.ALONE))))
            measures["oom_kills"] += 1
            if not outside and all(job.peak_mb <= job.claim.mb for job in running):
                broken.append(f"t={now:.0f}: the kernel killed {victim.name} though every job was inside its claim")
            running = [job for job in running if job is not victim]
        running = [job for job in jobs if job.state == "running"]
        for token in {t for job in running for t in job.claim.tokens if t != tokens.ALONE}:
            if sum(token in job.claim.tokens for job in running) > 1:
                broken.append(f"t={now:.0f}: {token} was held by two jobs at once")
        if any(tokens.ALONE in job.claim.tokens for job in running) and any(job.claim.lane in (memorygate.LANE_BULK, memorygate.LANE_MAINTENANCE) and tokens.ALONE not in job.claim.tokens for job in running):
            broken.append(f"t={now:.0f}: a job ran beside the one that runs alone")
        host = Host(total_mb, total_mb - BASE_MB - outside - sum(job.held_mb(now) for job in running), cores, min(len(running) / cores, 1.0))
        held = [Holder(job.name, job.claim, 0, "sim", job.held_mb(now), job.peak_mb) for job in running]
        waiting = [scheduler.Ticket(job.name, job.claim, job.queued, 0, "sim") for job in jobs if job.state == "ready" and job.queued <= now]
        if waiting:
            measures["turns"] += 1
            try:
                if crashes.random() < crash:
                    raise DeciderCrash(f"t={now:.0f}")
                admitted = rule(waiting, host, held, now, used, window, PLACES).admitted
            except DeciderCrash:
                measures["crashed_turns"] += 1
                admitted = []
            for name in admitted:
                job = by_name[name]
                job.waits.append(now - job.queued)
                job.state, job.started = "running", now
            claimed = sum(job.claim.mb for job in jobs if job.state == "running")
            if claimed > budget + 1e-6 and sum(job.state == "running" for job in jobs) > 1:     # a claim over the whole budget runs on its own, as the gate's solo claim does
                broken.append(f"t={now:.0f}: the claims admitted come to {claimed:.0f} MB, over the budget of {budget:.0f} MB")
    for job in jobs:
        if job.done_times != 1:
            broken.append(f"{job.name} finished {job.done_times} times")
    site = sorted(wait for job in jobs if job.claim.lane == memorygate.LANE_INTERACTIVE for wait in job.waits) or [0.0]
    measures["wait_p99_site"] = site[int(0.99 * (len(site) - 1))]
    calm = sorted(wait for job in jobs if job.claim.lane == memorygate.LANE_INTERACTIVE and not any(start <= job.queued < end for start, end, _mb in pressure) for wait in job.waits) or [0.0]
    measures["wait_p99_site_calm"] = calm[int(0.99 * (len(calm) - 1))]       # a visitor's job that comes while something outside holds the memory waits for the memory, which is right
    measures["wait_cut"] = max(by_name["cut-next-stage"].waits or [math.nan])
    for name in ("big-claim", "restarted-unit", "cut-next-stage"):
        job = by_name.get(name)
        if job is not None:         # how long after everything of the bulk and maintenance lanes queued before it had finished it started: later arrivals do not hold it back for ever
            ahead = [other.finished for other in jobs if other is not job and other.claim.lane in (memorygate.LANE_BULK, memorygate.LANE_MAINTENANCE) and other.queued <= job.queued and other.finished is not None]
            measures[f"late_{name}"] = max(job.waits or [math.nan]) - (max(ahead, default=job.queued) - job.queued)
    measures["wait_max_bulk"] = max(wait for job in jobs if job.claim.lane == memorygate.LANE_BULK for wait in job.waits)
    measures["makespan"] = now
    return measures, broken
