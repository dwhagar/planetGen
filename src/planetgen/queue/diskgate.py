# planetgen/queue/diskgate.py

"""
The disk controller of the supervisor (OPS.288, 8.1): a full disk holds the run with a plain reason instead of failing it with "its process was gone".

Three places fill with a long run and each stops it in its own way when it is full: the database's data directory (MySQL fails every write), the folder of the Redis server (its append-only file cannot
grow, and a Redis that cannot write refuses writes), and the engine's state folder (logs, unit files, the packed sector lists). The first of them was measured before a run (`generation.stats.database_disk`
refuses a run that cannot fit); none was watched while the run went on, so a run that filled the disk ended as a process that vanished.

**The signal.** Once a few seconds the controller measures the free space of each of the three places (the same drive counted once) against a **reserve**: the least that is left for the engine to write its own
hold with, `max(RESERVE_MIN_BYTES, RESERVE_SHARE of the drive)`, never more than `RESERVE_MAX_SHARE` of a small drive. Free space under twice the reserve is *pressed*, under the reserve *starved*.

**What it does** is the memory controller's own hold and release (`pressure.Controller`, reused: the same timing, the same acts in the Queue page's list, and the same `holding()` that the launcher and the light
jobs' duty ask before they admit, so nothing of the bulk or upkeep lanes starts while a drive is short). On top of that, while the disk is starved, the live generation run is **held** (`held.drain_run` with
`held_for="disk"`): its units stop at their next piece, keep what they committed, and the run says `DISK_REASON`: "the disk is full: free space, then Retry held or finish-run". Nothing carries it on by
itself: freeing space is a person's business, and the hold is only the admission's to let go of. The record never says the process was gone for a full disk (`gone_text`).
"""

import os
import shutil
import time
from typing import NamedTuple

from planetgen.queue import pressure
from planetgen.util import log
from planetgen.util.format import format_bytes

RESERVE_MIN_BYTES = 512 * 1024 * 1024
"""int: The least a drive keeps free for the engine's own writes (the hold, the run record, a log line). Stand-in: a design edge, not a measurement."""

RESERVE_SHARE = 0.005
"""float: The share of a drive the reserve is, when that is more than `RESERVE_MIN_BYTES`. Stand-in."""

RESERVE_MAX_SHARE = 0.1
"""float: The most of a small drive the reserve may be. Stand-in."""

PRESSED_FACTOR = 2.0
"""float: Free space under this many reserves holds admission (nothing new starts); under one reserve the run is held. Stand-in."""

LOOK_SECONDS = 5.0
"""float: How long a measure of the places is kept. Stand-in."""

WHY = "the disk is full"
DISK_REASON = "the disk is full: free space, then Retry held or finish-run"
"""str: What a run held for a full disk says (OPS.288)."""

HELD_FOR_DISK = "disk"


class Reading(NamedTuple):
    """One look at the worst place: its name and path, its free and total bytes, the reserve it was judged against and the `level` (`pressure.CALM`, `PRESSED` or `STARVED`)."""

    where: str
    path: str
    free_bytes: int
    total_bytes: int
    reserve_bytes: int
    level: int

    def words(self):
        if not self.path:
            return "no drive could be measured"
        return (f"{format_bytes(self.free_bytes)} free of {format_bytes(self.total_bytes)} on the drive of {self.where} ({self.path}), "
                f"the engine keeps {format_bytes(self.reserve_bytes)} back")


def reserve_bytes(total_bytes):
    """The reserve of a drive of `total_bytes`: `max(RESERVE_MIN_BYTES, RESERVE_SHARE of it)`, at most `RESERVE_MAX_SHARE` of a small drive."""
    return int(min(max(RESERVE_MIN_BYTES, RESERVE_SHARE * total_bytes), RESERVE_MAX_SHARE * total_bytes))


def level_of(free_bytes, total_bytes):
    """`pressure.STARVED` under the reserve, `PRESSED` under `PRESSED_FACTOR` reserves, else `CALM`."""
    reserve = reserve_bytes(total_bytes)
    return pressure.STARVED if free_bytes < reserve else pressure.PRESSED if free_bytes < PRESSED_FACTOR * reserve else pressure.CALM


def places(control=None, connection=None):
    """
    `[(name, path)]` of the places a run writes to that this machine can measure: the database's data directory (asked of the server, and only when the server runs here), the folder of the Redis server
    (when it runs here) and the state folder. A place that cannot be found is left out; none of them is a reason to fail.
    """
    from planetgen.queue.work import _is_local_host
    from planetgen.util import settings
    found = []
    try:
        if control is not None:
            from planetgen.db import store
            if _is_local_host(store.control_mysql_config().host):
                found.append(("the database", control.execute("SELECT @@datadir AS d").fetchone()["d"]))
    except Exception as exc:  # noqa: BLE001 -- a server that does not say is not measured
        log.swallowed("queue/diskgate", "places", exc)
    try:
        if connection is not None:
            host = connection.connection_pool.connection_kwargs.get("host", "")
            if not host or _is_local_host(host):
                found.append(("the run record's Redis server", connection.config_get("dir")["dir"]))
    except Exception as exc:  # noqa: BLE001 -- a server that does not say (or a refused CONFIG) is not measured
        log.swallowed("queue/diskgate", "places", exc)
    found.append(("the state folder", settings.state_dir()))
    return found


def measure(found, usage=None):
    """The `Reading` of the worst of `found` (`[(name, path)]`): the least free space against its reserve. Two places on one drive count once; a path that cannot be read is left out."""
    worst, seen = None, set()
    for name, path in found:
        try:
            probe = path
            while probe and not os.path.exists(probe):
                parent = os.path.dirname(probe)
                if parent == probe:
                    break
                probe = parent
            drive = os.stat(probe).st_dev
            if drive in seen:
                continue
            seen.add(drive)
            figures = (usage or shutil.disk_usage)(probe)
        except OSError as exc:
            log.swallowed("queue/diskgate", "measure", exc)
            continue
        reading = Reading(name, probe, figures.free, figures.total, reserve_bytes(figures.total), level_of(figures.free, figures.total))
        if worst is None or (reading.level, -reading.free_bytes / max(1, reading.reserve_bytes)) > (worst.level, -worst.free_bytes / max(1, worst.reserve_bytes)):
            worst = reading
    return worst or Reading("", "", 0, 0, 0, pressure.CALM)


class Look:
    """
    The controller's `look`: measures the places each time it is called, and finds them again (`places`) while none was found. `usage` and `find` are what a test replaces.
    """

    def __init__(self, control=None, connection=None, usage=None, find=None):
        self.control, self.connection, self.usage = control, connection, usage
        self.find = find or self._find
        self.found = None

    def _find(self):
        """The places, asked of the database through a connection of its own when `control` makes them (and closed again)."""
        if not callable(self.control):
            return places(self.control, self.connection)
        control = self.control()
        try:
            return places(control, self.connection)
        finally:
            control.close()

    def __call__(self):
        if not self.found:
            self.found = self.find()
        return measure(self.found, self.usage)


class Controller(pressure.Controller):
    """
    The disk's controller (see the module): the memory controller's hold and release with the disk for its `look`, and the live generation run held while a drive is starved. `connection` is the run records' Redis.
    `drain(connection, run_id, reason)` holds one run (`held.drain_run` with `held_for="disk"`); a test passes its own.
    """

    def __init__(self, connection, look, clock=time.time, drain=None):
        super().__init__(launcher=None, look=look, clock=clock, why=WHY)
        self.connection, self.drain = connection, drain
        self.held_runs = []

    def turn(self):
        asked = super().turn()
        if self.held_since is not None and self.last is not None and self.last.level >= pressure.STARVED:
            self._hold_runs()
        return asked

    def _hold_runs(self):
        """Holds the live generation run (`DISK_REASON`); each run is held once, and said."""
        from planetgen.queue import held, runrecord
        try:
            run_id = runrecord.active_run(self.connection, runrecord.GENERATION)
            if run_id is None or run_id in self.held_runs:
                return
            drain = self.drain or (lambda conn, run, reason: held.drain_run(conn, run, reason, held_for=HELD_FOR_DISK))
            if drain(self.connection, run_id, DISK_REASON):
                self.held_runs.append(run_id)
                self._act(f"held run {run_id}", f"{DISK_REASON} ({self.last.words()})")
        except Exception as exc:  # noqa: BLE001 -- a run that cannot be held now is held at the next turn
            log.swallowed("queue/diskgate", "_hold_runs", exc)

    def clear(self):
        super().clear()
        self.held_runs = []


def controller(connection, control=None):
    """The disk's controller for a supervisor: its look measures the database's drive (through `control`, a connection or a callable that gives one), the Redis server's and the state folder's."""
    return Controller(connection, Look(control=control, connection=connection))


def full_text():
    """What a run that ended with its process gone says instead when a drive is starved (`DISK_REASON`), or `None` when none is: the record must not blame a process for a full disk."""
    try:
        from planetgen.db import store
        control = store.get_control_connection(store.control_mysql_config())
        try:
            reading = measure(places(control))
        finally:
            control.close()
        return f"{DISK_REASON} ({reading.words()})" if reading.level >= pressure.STARVED else None
    except Exception as exc:  # noqa: BLE001 -- the usual text stands where the disk cannot be read
        log.swallowed("queue/diskgate", "full_text", exc)
        return None
