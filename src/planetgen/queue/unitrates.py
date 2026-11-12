# planetgen/queue/unitrates.py

"""
How fast the units of each command work, in their own items a second (PERF.147), kept so a unit's bar has an ETA from its first second (`taskprogress.progress`, `stored_rate`).

* **In the run.** When a unit finishes it notes the items it did and the seconds it worked (not counting the time it waited on a lock or the gate, `unitbeat.UnitBeat.totals`). The run
  keeps the sums per command (`runrecord.rates_key`, expiring with the run), so the rate of a command in this run is what the units that finished so far did.
* **Across runs.** A run that ends `done` folds its rates into the server's (`persisted_key`, beside `stagetimes` and `commandgaps`): half of the stored rate stays (`KEEP_WEIGHT`), so
  one busy night does not become the rule. A run that failed or was stopped adds nothing.
* **For the bar.** `rate_of(conn, run_id)` is the `rate_of` of a `unitwatch.UnitWatch`: for a command, the rate of this run's finished units when there are any (the current run counts
  for more than the stored figures), else the server's.
"""

from planetgen.queue import runrecord
from planetgen.util import log

KEEP_WEIGHT = 0.5
"""float: The share of a stored rate that stays when a run's rate is folded in."""

MIN_SECONDS = 1.0
"""float: A unit that worked for less than this says nothing about a rate and is not noted."""


def persisted_key():
    return runrecord._k("unit-rates")


def note(conn, run_id, kind, done, seconds):
    """Notes that a unit of command `kind` did `done` items in `seconds` of work; the run keeps the sums. A unit that did nothing, or worked for less than `MIN_SECONDS`, is not noted."""
    if not (done and done > 0 and seconds and seconds >= MIN_SECONDS):
        return
    key = runrecord.rates_key(run_id)
    pipe = conn.pipeline()
    pipe.hincrbyfloat(key, f"{kind}:done", float(done))
    pipe.hincrbyfloat(key, f"{kind}:seconds", float(seconds))
    pipe.expire(key, runrecord.KEEP_SECONDS)
    pipe.execute()


def read(conn, run_id):
    """`{command: items a second}` the run's finished units did so far."""
    sums = {runrecord._s(field): float(runrecord._s(value)) for field, value in conn.hgetall(runrecord.rates_key(run_id)).items()}
    return {name[:-len(":done")]: sums[name] / sums[f"{name[:-len(':done')]}:seconds"]
            for name in sums if name.endswith(":done") and sums.get(f"{name[:-len(':done')]}:seconds", 0.0) > 0.0}


def persisted(conn):
    """`{command: items a second}` the server has learned."""
    return {runrecord._s(kind): float(runrecord._s(value)) for kind, value in conn.hgetall(persisted_key()).items()}


def fold(conn, run_id):
    """Folds the rates of a run that ended `done` into the server's (see the module docstring); returns the new `persisted`."""
    stored = persisted(conn)
    for kind, rate in read(conn, run_id).items():
        old = stored.get(kind)
        stored[kind] = rate if old is None else KEEP_WEIGHT * old + (1.0 - KEEP_WEIGHT) * rate
        conn.hset(persisted_key(), kind, round(stored[kind], 6))
    return stored


def rate_of(conn, run_id):
    """A `rate_of(kind) -> items a second or None` for a `unitwatch.UnitWatch`: this run's so far, else the server's (read again at each call, so a unit that finished a moment ago counts)."""
    def rate(kind):
        try:
            return read(conn, run_id).get(kind) or persisted(conn).get(kind)
        except Exception as exc:       # noqa: BLE001 -- a Redis that cannot answer for one poll leaves the bar with the running rate alone
            log.swallowed("queue/unitrates", "rate_of.rate", exc)
            return None
    return rate
