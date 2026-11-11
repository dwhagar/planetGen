# planetgen/queue/beatjudge.py

"""
The heartbeat judge (OPS.137; reports/heartbeat-research-2026-11-09.md): one pure function from a unit's beat record, the history of its command and a
clock to a state, so the follower (OPS.138) can say what a running unit is doing and a person can see it. No Redis, no clock of its own (`now` is passed,
on the Redis server clock, which is also the clock the beats are written with), so every rule is driven by a fake clock in a test, the way `runeta` is.

What is judged is a **progress** heartbeat, not a liveness one. The horse writes a beat every `BEAT_SECONDS` from a thread (OPS.136) that keeps beating while
the work loop is frozen, and the work loop sets a counter (`done` of `estimate`), a `phase` and a `waiting_on` tag. "Beat fresh, counter frozen" is then
stuck, not dead.

The states, checked in this order:

* **dead**: no beat for `DEAD_BEATS` intervals (30 s). The process is gone.
* **waiting**: `waiting_on` is set (the gate, the sky lock, MySQL). Not stuck. It has its own, longer limit (`WAITING_LIMIT_SECONDS`); past it the state stays
  waiting and `over_limit` is set, because what it waits for may be a good job (a sky cell can take half an hour): the holder is what gets judged, not the waiter.
* **quiet** (stuck): the beat is fresh and the counter has not changed (not counting time spent waiting, which ends with the wait) for longer than `max(5 min, 5 x the longest gap between counter changes seen for this
  command)`, doubled while the host is starved (memory gate pressure, swap, load far over the cores).
* **looping**: the counter moves but is not converging. `done` is over 3 times `estimate` (severity `flag`) or over 10 times (`suspect`); or `done` has gone
  backward inside one phase; or the same phase has begun again more than `RESTART_LIMIT` times (`suspect`).
* **slow**: the counter moves, at less than a third of the expected rate. A flag with a projected finish.
* **working**: anything else.

A unit with no counter (no estimate, and `done` never changed) cannot be judged for progress: it is working until it is dead, with `counter` false in the
verdict. A unit that does not report is blind, not stuck.

**Only dead and quiet are `actionable`.** Boss's rule (research report, "Restart or fail"): a job whose counter is moving is never killed unless he says so, at
any rate and any size over its estimate, so looping and slow are flags for a person and never an action. The judge says what a unit is; that two checks in a
row must agree before anything is done is `agreed`.

    verdict = judge(beat, history, now)                  # Verdict(state="quiet", ...)
    if beatjudge.agreed(previous, verdict): ...          # two checks in a row say the same thing about a unit that may be acted on

The thresholds come from recorded gaps, with floors so a quick start is never flagged (the same shape as `runrecord.is_stale` and `runeta`'s stall rule).
"""

import math
from dataclasses import dataclass
from typing import Optional

BEAT_SECONDS = 10.0
"""float: How often a horse writes its beat (OPS.136)."""

DEAD_BEATS = 3
"""int: Beats missed in a row that make a unit dead."""

QUIET_FLOOR_SECONDS = 300.0
"""float: A counter unchanged for less than this is never quiet, whatever the history says."""

QUIET_GAP_FACTOR = 5.0
"""float: ... and quiet starts at this many times the longest gap between counter changes seen for the command."""

STARVED_FACTOR = 2.0
"""float: The quiet limit is multiplied by this while the host is starved: a slow host looks like a frozen unit."""

WAITING_LIMIT_SECONDS = 1800.0
"""float: A wait on the gate, the sky lock or MySQL longer than this is flagged (`over_limit`), never judged quiet."""

LOOP_FLAG_FACTOR = 3.0
LOOP_SUSPECT_FACTOR = 10.0
"""float: `done` over this many times `estimate` is looping, as a flag and as a suspect."""

RESTART_LIMIT = 3
"""int: The same phase beginning again more often than this is looping (suspect)."""

SLOW_FRACTION = 1.0 / 3.0
"""float: A rate under this share of the expected one is slow."""

SLOW_MIN_SECONDS = 60.0
"""float: A unit younger than this is not judged slow (a start-up says little about a rate)."""

WORKING, SLOW, WAITING, QUIET, LOOPING, DEAD = "working", "slow", "waiting", "quiet", "looping", "dead"
STATES = (WORKING, SLOW, WAITING, QUIET, LOOPING, DEAD)
"""tuple: Every state, in the order the admin page lists them."""

ACTIONABLE = frozenset({QUIET, DEAD})
"""frozenset: The states a restart or a failure may follow from (never a moving counter)."""

FLAG, SUSPECT = "flag", "suspect"


@dataclass(frozen=True)
class Beat:
    """
    What one running unit's record holds (OPS.136). Times are seconds on the Redis server clock.

    Attributes:
        beat (float): When the horse's thread last wrote the record.
        started (float): When the unit began.
        done (float): The work finished, in the unit's own items.
        estimate (float): The work expected, 0 or `None` for a unit that has no counter.
        moved (float, optional): When `done` last changed; `None` while it never has.
        phase (str): What it is doing ("drawing tier 3").
        waiting_on (str): `gate`, `sky lock`, `mysql` or empty.
        restarts (int): How many times this phase has begun again.
        waited (float): Seconds the unit has spent waiting (gate, sky lock, MySQL) and stopped since `moved`; a wait that has ended is not time the counter stood still.
        resident_mb (float, optional): The memory it holds, carried for the page; not judged.
    """
    beat: float
    started: float
    done: float = 0.0
    estimate: float = 0.0
    moved: Optional[float] = None
    phase: str = ""
    waiting_on: str = ""
    restarts: int = 0
    waited: float = 0.0
    resident_mb: Optional[float] = None


@dataclass(frozen=True)
class History:
    """
    What earlier work of the same command says (OPS.135).

    Attributes:
        longest_gap (float): The longest seconds between two counter changes seen for this command, 0 where none is recorded.
        expected_seconds (float, optional): How long a unit like this takes, from the earlier runs; with the unit's `estimate` it gives the expected rate.
    """
    longest_gap: float = 0.0
    expected_seconds: Optional[float] = None


@dataclass(frozen=True)
class Verdict:
    """
    The judge's answer.

    Attributes:
        state (str): One of `STATES`.
        severity (str): `FLAG` or `SUSPECT` where the state is a flag (looping and slow); empty otherwise.
        reason (str): One plain sentence with the numbers, for the log line and the page.
        counter (bool): Whether the unit has a counter to judge.
        rate (float, optional): Items finished a second since the unit began, where there is a counter that moved.
        projected_finish_seconds (float, optional): Seconds left at that rate, for slow and looping.
        still_seconds (float): How long the counter has not changed (0 where it never could be judged).
        over_limit (bool): A wait past `WAITING_LIMIT_SECONDS`.
    """
    state: str
    severity: str = ""
    reason: str = ""
    counter: bool = True
    rate: Optional[float] = None
    projected_finish_seconds: Optional[float] = None
    still_seconds: float = 0.0
    over_limit: bool = False

    @property
    def actionable(self):
        """Whether a restart or a failure may follow from this state: dead and quiet only, never a unit whose counter is moving."""
        return self.state in ACTIONABLE


def _finite(value):
    return value is not None and isinstance(value, (int, float)) and math.isfinite(value)


def dead_after_seconds():
    """Seconds without a beat after which a unit is dead."""
    return DEAD_BEATS * BEAT_SECONDS


def quiet_limit_seconds(longest_gap, starved=False):
    """
    How long a counter may stand still before the unit is quiet: `max(5 min, 5 x longest_gap)`, doubled while the host is starved.

    Args:
        longest_gap (float): The longest recorded gap between counter changes for this command, seconds (0 or `None` where none is recorded).
        starved (bool): Whether the host is starved.
    """
    gap = longest_gap if _finite(longest_gap) and longest_gap > 0 else 0.0
    limit = max(QUIET_FLOOR_SECONDS, QUIET_GAP_FACTOR * gap)
    return limit * STARVED_FACTOR if starved else limit


def _has_counter(beat):
    return (_finite(beat.estimate) and beat.estimate > 0) or beat.moved is not None


def _rate(beat, now):
    """Items a second since the unit began, or `None` where it has done nothing or has not run for any time."""
    ran = now - beat.started
    if not _finite(beat.done) or beat.done <= 0 or not _finite(ran) or ran <= 0:
        return None
    return beat.done / ran


def _projected(beat, rate):
    if rate is None or rate <= 0 or not (_finite(beat.estimate) and beat.estimate > beat.done):
        return None
    return (beat.estimate - beat.done) / rate


def _expected_rate(beat, history):
    if not (_finite(history.expected_seconds) and history.expected_seconds > 0 and _finite(beat.estimate) and beat.estimate > 0):
        return None
    return beat.estimate / history.expected_seconds


def _looping(beat, before, rate):
    """The looping verdict, or `None`: over the estimate, going backward in one phase, or the same phase restarting."""
    if _finite(beat.estimate) and beat.estimate > 0 and _finite(beat.done):
        ratio = beat.done / beat.estimate
        if ratio > LOOP_SUSPECT_FACTOR:
            return Verdict(LOOPING, SUSPECT, f"{beat.done:g} done of {beat.estimate:g} expected, {ratio:.1f} times the estimate and still counting", rate=rate)
        if ratio > LOOP_FLAG_FACTOR:
            return Verdict(LOOPING, FLAG, f"{beat.done:g} done of {beat.estimate:g} expected, {ratio:.1f} times the estimate", rate=rate)
    if before is not None and before.phase == beat.phase and _finite(before.done) and _finite(beat.done) and beat.done < before.done:
        return Verdict(LOOPING, SUSPECT, f"the count went back from {before.done:g} to {beat.done:g} inside {beat.phase or 'one phase'}", rate=rate)
    if beat.restarts > RESTART_LIMIT:
        return Verdict(LOOPING, SUSPECT, f"{beat.phase or 'the same phase'} has begun again {beat.restarts} times", rate=rate)
    return None


def judge(beat, history, now, before=None, starved=False):
    """
    What a unit is doing, from its beat, its command's history and the clock.

    Args:
        beat (Beat): The unit's latest record.
        history (History): The command's recorded gaps and expected time (a `History()` where nothing is recorded).
        now (float): The time, on the clock the beats are written with (the Redis server's).
        before (Beat, optional): An earlier beat of the same unit; with it, `done` going backward inside one phase is looping.
        starved (bool): Whether the host is starved (the memory gate under pressure, swap, load far over the cores): the quiet limit doubles.

    Returns:
        Verdict: The state, why, and whether it may be acted on.
    """
    silent = now - beat.beat
    if not _finite(silent) or silent > dead_after_seconds():
        return Verdict(DEAD, reason=f"no beat for {silent:.0f} s (limit {dead_after_seconds():.0f} s)" if _finite(silent) else "no usable beat",
                       counter=_has_counter(beat))
    has_counter = _has_counter(beat)
    last_change = beat.moved if beat.moved is not None else beat.started
    spent_waiting = beat.waited if _finite(beat.waited) and beat.waited > 0 else 0.0
    still = max(now - last_change - spent_waiting, 0.0) if _finite(last_change) else 0.0
    rate = _rate(beat, now) if has_counter else None
    if beat.waiting_on:
        waited = max(now - last_change, 0.0)
        over = waited > WAITING_LIMIT_SECONDS
        return Verdict(WAITING, FLAG if over else "", f"waiting on {beat.waiting_on} for {waited:.0f} s" + (f" (over {WAITING_LIMIT_SECONDS:.0f} s)" if over else ""),
                       counter=has_counter, rate=rate, still_seconds=waited, over_limit=over)
    if not has_counter:
        return Verdict(WORKING, reason="no counter reported, so only a missing beat can be judged", counter=False)
    limit = quiet_limit_seconds(history.longest_gap, starved)
    if still > limit:
        return Verdict(QUIET, reason=f"beat fresh but the count has stood at {beat.done:g} for {still:.0f} s (limit {limit:.0f} s)", rate=rate, still_seconds=still)
    looping = _looping(beat, before, rate)
    if looping is not None:
        return Verdict(looping.state, looping.severity, looping.reason, rate=rate, still_seconds=still,
                       projected_finish_seconds=_projected(beat, rate) if looping.severity == FLAG else None)
    expected = _expected_rate(beat, history)
    if expected is not None and rate is not None and now - beat.started >= SLOW_MIN_SECONDS and rate < SLOW_FRACTION * expected:
        projected = _projected(beat, rate)
        tail = f", about {projected:.0f} s left at that rate" if projected is not None else ""
        return Verdict(SLOW, FLAG, f"{rate:.3g} a second against about {expected:.3g} expected{tail}", rate=rate, projected_finish_seconds=projected,
                       still_seconds=still)
    return Verdict(WORKING, rate=rate, still_seconds=still)


def agreed(previous, current):
    """
    Whether two checks in a row say the same thing about a unit that may be acted on (the follower acts only then, never on one late beat).

    Args:
        previous (Verdict): The check before.
        current (Verdict): The check now.
    """
    return previous is not None and current is not None and current.actionable and previous.state == current.state
