# planetgen/queue/dispatch.py

"""
The dispatch rule over the whole backlog (OPS.227, handle D2 of `reports/job-engine-items-stage-3-6-2026-11-13.md`; design section 3.5 and the stage 0 findings of section 10.5).

On main every worker asks the gate for room alone (`memorygate.acquire`, `scheduler.Scheduler`), so the scheduler only ever sees the jobs that already hold a worker. The job engine's supervisor
sees every job that is ready and makes every hand-out; this is the one function it asks, once a turn: *which of the jobs that are ready go now?* It is pure (the host, the claims held and the
clock are given), so the stage 0 simulator, the tests and the supervisor all run the same code. Nothing calls it from a worker yet: stage 4's later items (OPS.228 on) move a lane at a time.

**The order.**

1. *The site's lanes first* (`memorygate.SHORT_LANES`: the interactive and long lanes), oldest first, each against the places kept for them (`memorygate.decide`, which gives these lanes
   half the host's reserve and keeps their places clear of the other lanes' claims). They do not age against the bulk work and are not ranked with it: ranked together, a visitor's job waited
   1,823 s at the 99th percentile in the simulator, because every unit queued at the start of a stage outranks a new job after about 450 s of aging; served first, under 1 s in every seed.
2. *Every other lane* ranked by `scheduler.rank` (class weight, aging, fair share) and judged by `scheduler.judge`, so a job goes when it is the highest that fits, a lower one that fits goes
   ahead of a higher one that does not (backfill), and a solo or exclusive claim is let by.

**A reservation lasts until the job it was made for runs.** Once a job has waited `scheduler.RESERVE_AFTER_SECONDS` without room, the room that frees up is kept for it and lower jobs stop slipping
past it. A worker that gives up after `memorygate.MAX_WAIT_SECONDS` and goes in by the table cap ends the reservation on main; the supervisor has no such exit (a job that never runs is not let in), so
the reservation has no end (`scheduler.judge(..., reserve_until=math.inf)`). The cut that needs the bulk lane to itself waited as long as the whole stage (3,218 s on average) without it.

**Only the top of the order is judged** (`WINDOW`, 64): a turn over 4,321 waiting jobs took 6 ms judging the top 64 and 665 ms judging all.

**Each job's last reason is kept.** `Turn.reasons` holds the reason of every job judged and not let go; `note_reasons` writes them onto the jobs' ledger rows (`attempts.note_wait`, which writes a
changed reason at once and the same one at most every `attempts.WAIT_NOTE_SECONDS`), where the Queue page reads them.

Nothing here changes what a run generates (the work seed does not depend on the order or on how many jobs run at once).
"""

import math
from typing import NamedTuple

from planetgen.queue import attempts, memorygate, scheduler
from planetgen.queue.memorygate import Holder
from planetgen.queue.scheduler import Ticket

WINDOW = 64
"""int: How many waiting jobs of each group are judged in a turn (the rest wait their turn: judging 4,000 took 0.7 s)."""

HOST_NAME = "dispatch"
"""str: The host name a holder admitted in this turn carries."""


class Turn(NamedTuple):
    """What a turn decided: the names of the jobs to start, in the order they were admitted, and `reasons`, the reason each job judged and not admitted was given."""

    admitted: list
    reasons: dict


def site_first(waiting):
    """`waiting` (tickets) of the site's lanes, oldest first (equal ages by name), and the others: `(site, rest)`."""
    site = sorted((ticket for ticket in waiting if ticket.claim.lane in memorygate.SHORT_LANES), key=lambda ticket: (ticket.since, ticket.name))
    return site, [ticket for ticket in waiting if ticket.claim.lane not in memorygate.SHORT_LANES]


def dispatch(waiting, host, held, now, used, window=WINDOW, places=None):
    """
    One turn of the rule (see the module): which of `waiting` (`scheduler.Ticket`s of every job that is ready) start now, given the `host`, the claims `held` (`memorygate.Holder`s of the jobs
    running), the clock `now` and what each kind has used lately (`used`, `scheduler.usage`). A job admitted in this turn is a holder that has taken none of its claim yet, so the next is judged
    against it. Returns a `Turn`.
    """
    held = list(held)
    site, rest = site_first(waiting)
    admitted, reasons = [], {}
    for number, group in enumerate((site, scheduler.rank(rest, now, used))):
        top = group[:window]
        pending = list(top)
        for ticket in top:
            others = [other for other in pending if other.name != ticket.name]
            if number == 0:
                decision = memorygate.decide(ticket.claim, host, held, places=places, name=ticket.name)
            else:
                decision = scheduler.judge(ticket, others, host, held, now, used, places=places, reserve_until=math.inf)
            if decision.admit:
                admitted.append(ticket.name)
                reasons.pop(ticket.name, None)
                pending = others
                held.append(Holder(ticket.name, ticket.claim, 0, HOST_NAME, 0.0, 0.0))
            else:
                reasons[ticket.name] = decision.reason
        for ticket in group[window:]:
            reasons[ticket.name] = f"{len(group) - window} job(s) are ahead of the {window} judged in this turn; this one waits its turn"
    return Turn(admitted, reasons)


def note_reasons(conn, turn, where, now=None):
    """
    Writes the reason of each job `turn` did not admit onto its ledger row (`attempts.note_wait`). `where` maps a job's name to `(run_id, row)`; a name that is not in it (a job with no row) is
    left out. Returns how many rows were written.
    """
    written = 0
    for name, reason in turn.reasons.items():
        place = where.get(name)
        if place is not None and attempts.note_wait(conn, place[0], place[1], reason, now=now):
            written += 1
    return written
