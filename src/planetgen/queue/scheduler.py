# planetgen/queue/scheduler.py

"""
The fair-share scheduler of the work queue (OPS.130, Boss 3:44 PM EST on 2026-11-08: "we want to use a task scheduler and sharing system time so
warmup or other commands do not block the site so all tasks that can run get a chance to run").

The gate (`memorygate`) answers one question, *does this job fit the host now?* (`memorygate.decide`, a pure function). When several jobs
wait for room the gate alone gives it to whichever worker polls first, so a warm-up that queued thousands of units could take every freed
megabyte from a visitor's job, and a big claim could wait behind small ones for ever. This file answers the other question, *which waiting job goes
next?*, on the same view of the host and the claims held (`Scheduler.decide` is asked in place of `memorygate.decide`, under the gate's lock).

**A waiting job is a ticket.** A worker that asks the gate for room writes a ticket (`planetgen:gate:tickets`: the claim, when it began to wait, its
process; renewed at every poll, dropped when the job is admitted or the worker is gone). The tickets are the queue the scheduler orders.

**The order** is the ticket's priority, highest first:

    priority = weight(class) + AGING_RATE x seconds waited - FAIR_PENALTY x share of the recent system time the job's kind has used

* *Class.* Work that serves the site goes first: the interactive lane (`LANE_WEIGHT`), the visitor's node prefetch (`SITE_KINDS`), the long lane, then the
  generation run's units and the maintenance lane's warm-ups and checks, then everything else. The reserved places of the interactive and long
  lanes keep their memory and cores against the other lanes (`memorygate.kept_mb`, `memorygate.kept_cores`), so the site's next job always has room.
* *Fair share.* Within a class, a job kind that has used less of the host's time recently goes before one that has used more. The use is each
  finished job's seconds times its cores, added to its kind's total (`charge`, from `memorygate.release`) and halved every `USAGE_HALF_LIFE` seconds, so
  a warm-up's thousands of units cannot keep a check or a run's cut from a turn.
* *Aging.* The rate is the gap between the heaviest and the lightest class over `AGING_BOUND_SECONDS`, so a job that has waited that long outranks a job of
  any class that has just arrived: **no job is kept waiting by the order for more than `AGING_BOUND_SECONDS`.** (How long it waits for room is the time
  the claims in its way take to end; the gate's own bound, `memorygate.MAX_WAIT_SECONDS`, lets it in by the table cap after that.)

**A job goes when it is the highest ticket that fits.** A ticket is held back (`Scheduler.decide`) while a higher ticket that fits now is waiting
(it takes the room on its next poll), and, once a higher ticket has waited `RESERVE_AFTER_SECONDS` without room, while that ticket is still inside the
gate's bound: the room that frees up is kept for it and smaller, lower jobs stop slipping past it (a large claim is not starved by small ones).
Otherwise a lower ticket that fits goes ahead of a higher one that does not (backfill), so the host is not left idle for the sake of order.

Nothing here changes what a run generates (the work seed does not depend on the order or the number of jobs running at once).
"""

import json
import socket
import time
from typing import NamedTuple

from planetgen.queue import memorygate
from planetgen.queue.memorygate import Claim, Decision, Holder

LANE_WEIGHT = {
    memorygate.LANE_INTERACTIVE: 100.0,
    memorygate.LANE_LONG: 80.0,
    memorygate.LANE_BULK: 40.0,
    memorygate.LANE_MAINTENANCE: 30.0,
    memorygate.LANE_OTHER: 30.0,
}
"""dict: The class weight of each lane. A first setting: the order of the classes is the rule (the site first, a run's units before the upkeep), the numbers set how
far apart they are against `FAIR_PENALTY` and `AGING_RATE`."""

SITE_KINDS = {"fn:planetgen.web.lib.tileprefetch.build_nodes": 90.0}
"""dict: A kind of job on a lane of upkeep that serves a visitor all the same (the nodes of the page they opened), with the weight it goes by."""

AGING_BOUND_SECONDS = 600.0
"""float: After this long a waiting job outranks a job of any class that has just arrived (the worst case the order itself adds)."""

FAIR_PENALTY = 10.0
"""float: The most weight a kind loses for having used all of the recent system time (half of it for half)."""

AGING_RATE = (max(LANE_WEIGHT.values()) - min(LANE_WEIGHT.values()) + FAIR_PENALTY) / AGING_BOUND_SECONDS
"""float: The weight a waiting job gains each second: the gap between the heaviest and lightest class and the whole fair-share penalty, over the bound, so the
bound holds even for the lightest class that has used the most against a fresh job of the heaviest class that has used nothing."""

USAGE_HALF_LIFE = 3600.0
"""float: Seconds after which the system time a kind has used counts half."""

RESERVE_AFTER_SECONDS = 120.0
"""float: How long a ticket waits for room before the room that frees up is kept for it."""

TICKET_SECONDS = 30.0
"""float: How long a ticket lives without its worker renewing it (a poll is `memorygate.POLL_SECONDS`)."""


class Ticket(NamedTuple):
    """A job waiting for room: who (`name`, `pid`, `host`), the `claim` and when it began to wait (`since`, epoch seconds)."""

    name: str
    claim: Claim
    since: float
    pid: int
    host: str


def weight(claim):
    """The class weight of a claim: its kind's, if the kind serves the site (`SITE_KINDS`), else its lane's."""
    return SITE_KINDS.get(claim.kind, LANE_WEIGHT.get(claim.lane, min(LANE_WEIGHT.values())))


def priority(ticket, now, used):
    """The ticket's priority (see the module): the class weight, plus what waiting has earned, less the fair-share penalty of its kind. `used` is `usage`."""
    total = sum(used.values())
    share = used.get(ticket.claim.kind, 0.0) / total if total > 0.0 else 0.0
    return weight(ticket.claim) + AGING_RATE * max(0.0, now - ticket.since) - FAIR_PENALTY * share


def usage(conn, now=None):
    """The system time each kind has used lately, core-seconds, with the half-life applied (`{kind: seconds}`)."""
    now = time.time() if now is None else now
    found = {}
    for kind, raw in conn.hgetall(memorygate._key("usage")).items():
        try:
            row = json.loads(memorygate._decode(raw))
            found[memorygate._decode(kind)] = float(row["v"]) * 0.5 ** (max(0.0, now - float(row["t"])) / USAGE_HALF_LIFE)
        except (ValueError, KeyError, TypeError):
            continue
    return found


def charge(conn, kind, core_seconds, now=None):
    """Adds `core_seconds` to what `kind` has used (a finished job's seconds times its cores). Never raises: the order is a nicety, the job is done."""
    now = time.time() if now is None else now
    try:
        raw = conn.hget(memorygate._key("usage"), kind)
        old = 0.0
        if raw is not None:
            row = json.loads(memorygate._decode(raw))
            old = float(row["v"]) * 0.5 ** (max(0.0, now - float(row["t"])) / USAGE_HALF_LIFE)
        conn.hset(memorygate._key("usage"), kind, json.dumps({"v": old + max(0.0, core_seconds), "t": now}))
    except Exception:  # noqa: BLE001 -- bookkeeping
        pass


def tickets(conn, now=None):
    """The jobs waiting for room, as `Ticket`s. A ticket whose time ran out, or whose process is gone from this host, is dropped here."""
    now = time.time() if now is None else now
    found = []
    for name, raw in conn.hgetall(memorygate._key("tickets")).items():
        name = memorygate._decode(name)
        try:
            row = json.loads(memorygate._decode(raw))
            claim = Claim(row["kind"], float(row["mb"]), float(row["cpu"]), row["lane"], bool(row.get("exclusive")))
            ticket = Ticket(name, claim, float(row["since"]), int(row["pid"]), row["host"])
            expires = float(row["expires"])
        except (ValueError, KeyError, TypeError):
            conn.hdel(memorygate._key("tickets"), name)
            continue
        if expires < now or (ticket.host == _here() and not memorygate.alive(ticket.pid)):
            conn.hdel(memorygate._key("tickets"), name)
            continue
        found.append(ticket)
    return found


def _here():
    return socket.gethostname()


def rank(waiting, now, used):
    """`waiting` (tickets) highest priority first; the earlier ticket first when two are equal."""
    return sorted(waiting, key=lambda ticket: (-priority(ticket, now, used), ticket.since, ticket.name))


def judge(me, others, host, held, now, used, places=None):
    """
    The scheduler's rule, pure: may the ticket `me` be admitted now, given the other waiting tickets, the host and the claims held? The gate's own rule
    (`memorygate.decide`) must pass first; then no higher ticket may fit now, and none that has waited `RESERVE_AFTER_SECONDS` (inside the gate's bound)
    may be held up by this job's taking room.
    """
    base = memorygate.decide(me.claim, host, held, own_mb=0.0, places=places, name=me.name)
    if not base.admit:
        return base
    mine = priority(me, now, used)
    for other in rank(others, now, used):
        if other.name == me.name or priority(other, now, used) <= mine:
            continue
        waited = now - other.since
        fits = memorygate.decide(other.claim, host, held, places=places, name=other.name)
        if fits.admit:
            after = memorygate.decide(other.claim, host, held + [Holder(me.name, me.claim, me.pid, me.host, 0.0, 0.0)], places=places, name=other.name)
            if after.admit:
                continue                # room for both: the lower job does not hold the higher one up
            return Decision(False, f"{other.claim.kind} ({other.claim.lane}, ahead by priority) fits now and goes first", base.need_mb, base.free_mb)
        if RESERVE_AFTER_SECONDS <= waited < memorygate.MAX_WAIT_SECONDS:
            return Decision(False, f"room is kept for {other.claim.kind} ({other.claim.lane}), which has waited {waited:,.0f} s", base.need_mb, base.free_mb)
    return base


class Scheduler:
    """The scheduler `memorygate.acquire` is given: it keeps this worker's ticket and answers `judge` in place of the gate's plain rule."""

    def __init__(self, places=None, clock=time.time):
        self.places = places
        self.clock = clock

    def decide(self, conn, claim, name, host, held, own_mb, waited, pid):
        """Writes (or renews) the ticket of `name`, answers the scheduler's rule for it and, when the job may go, drops the ticket."""
        now = self.clock()
        waiting = tickets(conn, now)
        mine = next((ticket for ticket in waiting if ticket.name == name), None)
        since = mine.since if mine is not None else now - waited
        ticket = Ticket(name, claim, since, pid, _here())
        conn.hset(memorygate._key("tickets"), name, json.dumps({
            "kind": claim.kind, "mb": claim.mb, "cpu": claim.cpu, "lane": claim.lane, "exclusive": claim.exclusive, "since": since, "pid": pid,
            "host": ticket.host, "expires": now + TICKET_SECONDS}))
        others = [other for other in waiting if other.name != name]
        if not others:
            decision = memorygate.decide(claim, host, held, own_mb=own_mb, places=self.places, name=name)
        else:
            decision = judge(ticket, others, host, held, now, usage(conn, now), places=self.places)
            if decision.admit and own_mb > 0.0:
                decision = memorygate.decide(claim, host, held, own_mb=own_mb, places=self.places, name=name)       # the credit for what the worker holds
        if decision.admit:
            conn.hdel(memorygate._key("tickets"), name)
        return decision

    def done(self, conn, name):
        """The job is done waiting (admitted, let in by the table cap, or the worker is leaving): its ticket goes."""
        try:
            conn.hdel(memorygate._key("tickets"), name)
        except Exception:  # noqa: BLE001 -- the ticket expires by itself
            pass


def snapshot(conn, now=None):
    """What the scheduler holds, for a page or a report: the tickets in the order they would go, with their priority, and the use of each kind."""
    now = time.time() if now is None else now
    used = usage(conn, now)
    return {"waiting": [(ticket, round(priority(ticket, now, used), 2)) for ticket in rank(tickets(conn, now), now, used)], "usage": used}
