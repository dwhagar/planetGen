# planetgen/queue/pressure.py

"""
The one memory controller of the supervisor (OPS.231, D6 of the job engine): it stops admitting on lasting memory pressure, and then asks the youngest bulk job to stop at its next piece.

The memory gate (`memorygate.decide`) admits a job when the claims it holds leave room. That is a count of what the supervisor has promised, and the kernel can still run short: a job grows past its
claim, something else on the host takes memory. This is the signal that says so, and what the supervisor does about it before the kernel kills a job (the kill stays the last line; a killed job is
retried alone with a grown claim, as it was).

**The signal.** Where the kernel reports pressure (`/proc/pressure/memory`, PSI) the controller reads it: *some*, the share of the last ten seconds in which a task waited for memory, and *full*, the
share in which every task did. Where it does not (a container, macOS) it reads the host's available memory against what the gate keeps back (`memorygate.reserve_mb`). Either way
the reading is a level: 0 calm, 1 pressed (`HOLD_SOME`, or available memory under `PRESSURE_SHARE` of the reserve), 2 starved (`YIELD_SOME` or `YIELD_FULL`, or under `STARVED_SHARE`).

**What it does.**

* *Hold.* Pressed for `HOLD_AFTER_SECONDS` in a row: nothing of the bulk or maintenance lanes is admitted (`holding` gives the reason; the launcher and the light jobs' duty keep the site's lanes going and say
  why the rest waits). It lets go when the host has been calm for `RELEASE_AFTER_SECONDS`.
* *Yield.* Starved for `YIELD_AFTER_SECONDS` in a row, held all that while: it asks the youngest bulk job (the one that has done the least) to stop at its next piece (`Launcher.yield_youngest`, the pipe's
  `stop` request). A command that resumes (`commands.CommandSpec.resumes`) keeps the pieces it finished, and is asked before one that does not; the row goes back to *ready* and the stop does not count against its retry (`ledger.STOPPED`). One job at a time, `YIELD_COOLDOWN_SECONDS` apart, so the first has given its
  memory back before the next is judged.

Every change of state and every ask is said (`Supervisor.act`, kept in the list the Queue page shows). The edges are design stand-ins, not measurements: H4 of the job engine has them learned from recorded
runs. Nothing here changes what a run generates (the work seed does not depend on how many jobs run at once). The control is off with the gate (`memorygate.DISABLE_ENV`).
"""

import time
from typing import NamedTuple

from planetgen.queue import memorygate
from planetgen.util import log
from planetgen.util.format import format_duration_seconds, format_megabytes

PSI_FILE = "/proc/pressure/memory"
"""str: The kernel's memory pressure file (PSI)."""

HOLD_SOME = 10.0
"""float: `some` (percent of the last ten seconds a task waited for memory) at or over which the host is pressed. Stand-in: a design edge, not a measurement."""

YIELD_SOME = 25.0
"""float: `some` at or over which the host is starved. Stand-in."""

YIELD_FULL = 5.0
"""float: `full` (percent in which every task waited) at or over which the host is starved. Stand-in."""

PRESSURE_SHARE = 0.75
"""float: Without PSI, the host is pressed when its available memory is under this share of `memorygate.reserve_mb`. Stand-in."""

STARVED_SHARE = 0.5
"""float: Without PSI, the host is starved when its available memory is under this share of the reserve. Stand-in."""

HOLD_AFTER_SECONDS = 10.0
"""float: How long the host is pressed in a row before the controller stops admitting. Stand-in."""

YIELD_AFTER_SECONDS = 30.0
"""float: How long the host is starved in a row, held all that while, before the youngest bulk job is asked to stop. Stand-in."""

YIELD_COOLDOWN_SECONDS = 60.0
"""float: The least time between two asks: long enough for the job asked to stop to give its memory back. Stand-in."""

RELEASE_AFTER_SECONDS = 30.0
"""float: How long the host is calm in a row before the controller admits again. Stand-in."""

CALM, PRESSED, STARVED = 0, 1, 2
LEVELS = {CALM: "calm", PRESSED: "pressed", STARVED: "starved"}


class Reading(NamedTuple):
    """One look at the host: where it came from (`psi` or `memory`), the two PSI figures (`None` without PSI), the available memory against the reserve (MB) and the `level`."""

    source: str
    some: float
    full: float
    available_mb: float
    reserve_mb: float
    level: int

    def words(self):
        if self.source == "psi":
            return f"memory pressure: tasks waited for memory {self.some:g}% of the last ten seconds, all of them {self.full:g}%"
        return f"{format_megabytes(self.available_mb)} of memory available, the gate keeps {format_megabytes(self.reserve_mb)} back"


def parse_psi(text):
    """`(some, full)` avg10 of the text of the kernel's pressure file, or `None` when it does not read so."""
    found = {}
    for line in text.splitlines():
        parts = line.split()
        if not parts or parts[0] not in ("some", "full"):
            continue
        fields = dict(item.split("=", 1) for item in parts[1:] if "=" in item)
        try:
            found[parts[0]] = float(fields["avg10"])
        except (KeyError, ValueError):
            return None
    return (found["some"], found.get("full", 0.0)) if "some" in found else None


def read(read_host=memorygate.read_host, path=PSI_FILE, edges=None):
    """The `Reading` now: PSI where the kernel has it, the host's available memory otherwise. A test passes its own `read_host` and `path`. `edges` are the `pressurelearn.Edges` it judges by (OPS.249); `None` is the typed ones above."""
    hold_some, yield_some, yield_full, pressure_share, starved_share = (HOLD_SOME, YIELD_SOME, YIELD_FULL, PRESSURE_SHARE, STARVED_SHARE) if edges is None else edges[:5]
    host = read_host()
    reserve = memorygate.reserve_mb(host)
    try:
        with open(path, encoding="utf-8") as handle:
            figures = parse_psi(handle.read())
    except OSError:
        figures = None
    if figures is not None:
        some, full = figures
        level = STARVED if some >= yield_some or full >= yield_full else PRESSED if some >= hold_some else CALM
        return Reading("psi", some, full, host.available_mb, reserve, level)
    level = STARVED if host.available_mb < starved_share * reserve else PRESSED if host.available_mb < pressure_share * reserve else CALM
    return Reading("memory", 0.0, 0.0, host.available_mb, reserve, level)


class Controller:
    """
    The supervisor's duty that watches memory (see the module). `launcher` is a callable that returns the `Launcher` whose youngest bulk job it asks to stop (or `None` while there is none: the launcher is made on the first turn that needs it); `look` returns a `Reading` (`read`, or a test's own); `clock` gives the time.
    `holding()` is what the launcher and the light jobs' duty ask before they admit: the reason the host is held, or `None`.
    """

    def __init__(self, launcher=None, look=read, clock=time.time, why="the host is short of memory", recorder=None):
        self.launcher, self.look, self.clock = launcher or (lambda: None), look, clock
        self.recorder = recorder        # OPS.249: `pressurelearn.Recorder`: it keeps the minutes and gives the learned edges to the default look
        self.edges = None               # the edges the last look was judged by
        self.why = why    # the reason the hold gives: the database's signal (`dbsignal`, OPS.248) is this same machine with its own look
        self.supervisor = None
        self.held_since = None      # when the hold began
        self.pressed_since = None   # when the level last rose to pressed or more (None while calm)
        self.starved_since = None
        self.calm_since = None
        self.asked_at = None
        self.last = None            # the last Reading

    def attach(self, supervisor):
        self.supervisor = supervisor

    def _act(self, what, reason=""):
        if self.supervisor is not None:
            self.supervisor.act(what, reason)

    def _look(self):
        """The reading now: by the learned edges (`pressurelearn`) when the default look is used and a recorder gives them; a look of the caller's own is used as it is. A change of the edges is said."""
        if self.recorder is None or self.look is not read:
            return self.look()
        edges = self.recorder.edges()
        if self.edges is not None and edges.sensitivity != self.edges.sensitivity:
            self._act("moved the memory edges", edges.words())
        self.edges = edges
        return read(edges=edges)

    def holding(self):
        """Why nothing of the bulk or maintenance lanes is admitted, or `None`."""
        if self.held_since is None:
            return None
        return f"admission is held: {self.why} ({self.last.words() if self.last else 'pressure'}); the site's jobs still go"

    def turn(self):
        """One look: moves the state by what the host reads, and asks a job to stop when the rule says so. Returns the job asked (its id) or `None`."""
        if not memorygate.enabled():
            self.held_since = None
            return None
        now = self.clock()
        try:
            reading = self._look()
        except Exception as exc:  # noqa: BLE001 -- a reading that cannot be had leaves the state as it is
            log.swallowed("queue/pressure", "turn", exc)
            return None
        self.last = reading
        if self.recorder is not None and self.look is read:
            self.recorder.note(reading, self.held_since is not None)
        if reading.level >= PRESSED:
            self.calm_since = None
            self.pressed_since = self.pressed_since if self.pressed_since is not None else now
        else:
            self.pressed_since = None
            self.calm_since = self.calm_since if self.calm_since is not None else now
        self.starved_since = (self.starved_since if self.starved_since is not None else now) if reading.level >= STARVED else None
        if self.held_since is None:
            if self.pressed_since is not None and now - self.pressed_since >= HOLD_AFTER_SECONDS:
                self.held_since = now
                self._act("stopped admitting bulk and upkeep jobs", f"{reading.words()} for {format_duration_seconds(now - self.pressed_since)}")
            return None
        if self.calm_since is not None and now - self.calm_since >= RELEASE_AFTER_SECONDS:
            self._act("admitting again", f"the host has been calm for {format_duration_seconds(now - self.calm_since)}; held for {format_duration_seconds(now - self.held_since)}")
            self.held_since = self.asked_at = None
            return None
        if self.starved_since is None or now - self.starved_since < YIELD_AFTER_SECONDS or self.launcher() is None:
            return None
        if self.asked_at is not None and now - self.asked_at < YIELD_COOLDOWN_SECONDS:
            return None
        asked = self.launcher().yield_youngest()
        if asked is not None:
            self.asked_at = now
            if self.recorder is not None and self.look is read:
                self.recorder.note(reading, True, asked=True)
            job_id, command = asked
            self._act(f"asked job {job_id} ({command}) to stop at its next piece", f"{reading.words()} for {format_duration_seconds(now - self.starved_since)}; it keeps its finished pieces and starts again when there is room")
        return None if asked is None else asked[0]

    def clear(self):
        """The supervisor ends: nothing to let go of."""
        self.held_since = None
