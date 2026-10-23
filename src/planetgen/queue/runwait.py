# planetgen/queue/runwait.py

"""
A run that finds the machine taken waits in line instead of failing (PERF.153, the interim of PERF.154): Boss,
2026-10-22: "a task that arrives is never lost or failed only because another is running". Used by
`wrapper.run` for a `GENERATION` run, whoever starts it: a Generate page job's stage that finds the machine held by
an API job (a one-sector regenerate, a neighbourhood), and the API job that finds a galaxy stage live.

The line is one sorted set in Redis, `planetgen:run:waiters` (first come, first served: the score is the arrival
time). A waiter holds its place with a sign of life (`ALIVE_SECONDS`), so one that was killed drops out of the line
by itself. Only the first in line tries the machine, so a late comer never jumps the line. It gives up when
`queue.wait_minutes` have passed (0 never waits) and fails naming the run it waited for.

PERF.154's single waiting list for every start path replaces this module.
"""

import secrets
import time

from planetgen.queue import runrecord

ALIVE_SECONDS = 30
"""int: How long a waiter's place is kept without a sign of life (it gives one at every poll)."""

POLL_SECONDS = 1.0
"""float: How often a waiter looks again."""


def _waiters_key():
    return runrecord._k("run", "waiters")


def _alive_key(ticket):
    return runrecord._k("run", "waiters", ticket)


def _prune(conn):
    """Drops from the line each waiter whose sign of life has run out."""
    for member in conn.zrange(_waiters_key(), 0, -1):
        ticket = runrecord._s(member)
        if not conn.exists(_alive_key(ticket)):
            conn.zrem(_waiters_key(), member)


def line(conn):
    """The tickets of the runs waiting for the machine, first in line first."""
    _prune(conn)
    return [runrecord._s(member) for member in conn.zrange(_waiters_key(), 0, -1)]


def start_when_free(conn, take, limit_seconds, on_wait=None, kind=runrecord.GENERATION, clock=time.time,
                    sleep=time.sleep):
    """
    Takes the machine with `take()` (it returns the run id, or raises `runrecord.Busy`), waiting in line for it when
    another run holds it.

    Args:
        conn: The Redis client.
        take (callable): Starts the run; raises `runrecord.Busy` when the machine is taken.
        limit_seconds (float): The longest to wait; 0 or less does not wait (the first `Busy` is raised).
        on_wait (callable, optional): `on_wait(holder, position)` each time the wait changes (the run that holds the
            machine, and the place in line, 1 being next), and `on_wait(None, None)` once the wait is over.
        kind (str): The kind of run, to name the holder.

    Returns:
        str: What `take` returned.

    Raises:
        runrecord.Busy: The machine stayed taken for `limit_seconds`; the message names the run.
    """
    if limit_seconds <= 0 or not line(conn):
        try:
            return take()
        except runrecord.Busy:
            if limit_seconds <= 0:
                raise
    ticket = secrets.token_hex(8)
    conn.zadd(_waiters_key(), {ticket: clock()})
    conn.set(_alive_key(ticket), "1", ex=ALIVE_SECONDS)
    deadline = clock() + limit_seconds
    shown = None
    try:
        while True:
            conn.set(_alive_key(ticket), "1", ex=ALIVE_SECONDS)
            waiting = line(conn)
            position = waiting.index(ticket) + 1 if ticket in waiting else 1
            holder = runrecord.active_run(conn, kind)
            if position == 1:
                try:
                    return take()
                except runrecord.Busy as busy:
                    holder = runrecord.active_run(conn, kind) or holder
                    last = busy
            else:
                last = runrecord.Busy(f"run {holder} is running")
            if clock() >= deadline:
                raise runrecord.Busy(f"{last}; waited {limit_seconds / 60:g} minutes for it")
            if on_wait is not None and (holder, position) != shown:
                shown = (holder, position)
                on_wait(holder, position)
            sleep(POLL_SECONDS)
    finally:
        conn.zrem(_waiters_key(), ticket)
        conn.delete(_alive_key(ticket))
        if shown is not None and on_wait is not None:
            on_wait(None, None)
