# planetgen/queue/commandgaps.py

"""
The longest time a command's counter has stood still between two changes (OPS.135), kept so the judge's quiet rule has a history (`beatjudge.quiet_limit_seconds`: 5 times the
longest gap seen for the command, never under 5 minutes).

* **In the run.** When a unit finishes it notes the longest gap its counter had (`unitbeat.UnitBeat.longest_gap`: the seconds between two changes of `done`, from the start of the unit to its
  end, less the time spent waiting on a lock or the gate, which is not standing still). The run keeps the largest per command (`runrecord.gaps_key`, a hash of command to seconds, expiring
  with the run), so the record of a finished run holds the longest counter gap of each command it ran.
* **Across runs.** A run that ends `done` folds its gaps into the server's (`persisted_key`, beside `stagetimes`): a gap longer than the one stored replaces it, a shorter one pulls it down by
  half (`KEEP_WEIGHT`), so one slow night does not loosen the rule for ever. A run that failed or was stopped adds nothing: its stuck units never finished, and what finished is not the
  longest to learn from.
* **For the judge.** `history_of(conn, run_id)` is the `history_of` of a `unitwatch.UnitWatch`: for a command, the larger of the server's gap and this run's so far.
"""

from planetgen.queue import beatjudge, runrecord

KEEP_WEIGHT = 0.5
"""float: The share of a stored gap that stays when a shorter one is folded in."""

_MAX = """
local old = tonumber(redis.call('HGET', KEYS[1], ARGV[1]) or 0)
if tonumber(ARGV[2]) > old then redis.call('HSET', KEYS[1], ARGV[1], ARGV[2]) end
redis.call('EXPIRE', KEYS[1], tonumber(ARGV[3]))
return 1
"""


def persisted_key():
    return runrecord._k("command-gaps")


def note(conn, run_id, kind, gap):
    """Notes that a unit of command `kind` had a counter gap of `gap` seconds; the run keeps the largest. A gap of nothing is not noted."""
    if gap is None or not gap > 0.0:
        return
    conn.register_script(_MAX)(keys=[runrecord.gaps_key(run_id)], args=[kind, round(float(gap), 3), runrecord.KEEP_SECONDS])


def read(conn, run_id):
    """`{command: seconds}` the run holds."""
    return {runrecord._s(kind): float(runrecord._s(value)) for kind, value in conn.hgetall(runrecord.gaps_key(run_id)).items()}


def persisted(conn):
    """`{command: seconds}` the server has learned."""
    return {runrecord._s(kind): float(runrecord._s(value)) for kind, value in conn.hgetall(persisted_key()).items()}


def fold(conn, run_id):
    """Folds the gaps of a run that ended `done` into the server's (see the module docstring); returns the new `persisted`."""
    stored = persisted(conn)
    for kind, gap in read(conn, run_id).items():
        old = stored.get(kind)
        stored[kind] = gap if old is None else max(gap, KEEP_WEIGHT * old + (1.0 - KEEP_WEIGHT) * gap)
        conn.hset(persisted_key(), kind, round(stored[kind], 3))
    return stored


def history_of(conn, run_id):
    """A `history_of(kind) -> beatjudge.History` for a `unitwatch.UnitWatch`: the larger of the server's gap and this run's so far (read again at each call, so a unit that finished
    a moment ago already counts)."""
    def history(kind):
        try:
            gap = max(persisted(conn).get(kind, 0.0), read(conn, run_id).get(kind, 0.0))
        except Exception:       # noqa: BLE001 -- a Redis that cannot answer for one poll leaves the judge with the floor
            gap = 0.0
        return beatjudge.History(longest_gap=gap)
    return history
