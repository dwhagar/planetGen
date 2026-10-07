# stellarObjects/generationStats.py

"""
How fast this server generates and how much space it takes (PERF.10),
and the size and time estimate every bulk generation shows before it
writes anything (PERF.3).

Speed is kept per stellar density, because a dense sector costs much
more than a sparse one. Densities (`relative_density`, 1.0 = the Sun's
neighborhood) fall into log-scale buckets, `BUCKETS_PER_DECADE` to a
decade from `MIN_DENSITY` up, with no top edge: a denser sector than any
seen before simply opens a new bucket. Every sector a run fills adds its
time to its bucket as a decaying average (`DECAY`, the weight of each
new sector), so the numbers follow the server as it gets faster or
busier, for as long as the control database exists. A reset galaxy keeps
them: they describe the machine, not the galaxy.

Times are the wall time a worker spent on one task, so with several
workers a run takes about the sum of its tasks' times divided by the
worker count; the contention of running side by side is already in the
measured times.

Size is measured from the galaxy database itself (its tables' data and
index bytes over its star systems), after every run, per galaxy
database.

Everything here fails open: without the control database (not created
yet, or no grant on it) the estimate uses `DEFAULT_SECONDS_PER_SYSTEM`
and `DEFAULT_BYTES_PER_SYSTEM`, and nothing is recorded.
"""

import math
import os
import time
from dataclasses import dataclass, field

from planetgen.util import log

MIN_DENSITY = 0.01
"""float: The bottom edge of bucket 0; anything sparser counts as 0."""

BUCKETS_PER_DECADE = 2
"""int: Buckets per factor of ten in density (half-decades)."""

DECAY = 0.05
"""float: The weight of each new sector in its bucket's decaying average
(about the last 20 sectors count for most of it)."""

SIZE_MARGIN = 0.10
"""float: Added to every size estimate (Boss: "+10%")."""

MAX_DISK_SHARE = 0.25
"""float: A run may not take more than this share of the disk's total."""

MIN_FREE_BYTES = 5 * 10 ** 9
"""int: A run may not leave less than this free (5 GB)."""

DEFAULT_SECONDS_PER_SYSTEM = 0.2
"""float: Worker time per star system before anything was measured.
Measured 2026-10-01 on a 4-core test server with MariaDB on the same
machine: 0.04 s with one worker, 0.17 s each with three (they wait on
the database); rounded up."""

DEFAULT_BYTES_PER_SYSTEM = 70 * 1024
"""int: Database bytes per star system before the galaxy has any.
Measured 2026-10-01: 10,545 systems in 692 MB (66 KB each, over half
of it moons)."""

STATS_ENV_VAR = "PLANETGEN_GENERATION_STATS"
"""str: `0` keeps `generate.py` from reading or recording any stats (the
test suite sets it, so test runs never touch a real server's numbers)."""

FLUSH_SECONDS = 30
"""int: How often a run writes its new averages back while it runs."""

KINDS = ("sector", "scatter")
"""tuple: What a bucket times: a sector's fill, or (PERF.9) a layer of
the bright-star scatter, by the layer's mean density."""


MAX_BUCKET = BUCKETS_PER_DECADE * 12
"""int: The top bucket (density 10^10 and up, infinite included): a
sector that dense is full long before its count, so they all time alike."""


def bucket_index(density):
    """The bucket `density` falls in (0 for anything below
    `MIN_DENSITY`, or not a number; `MAX_BUCKET` at most)."""
    if density is None or not density > MIN_DENSITY:
        return 0
    if not math.isfinite(density):
        return MAX_BUCKET
    index = math.floor(BUCKETS_PER_DECADE * (math.log10(density) - math.log10(MIN_DENSITY)) + 1e-9)
    return int(min(index, MAX_BUCKET))


def bucket_bounds(index):
    """`(low, high)` density edges of bucket `index`."""
    return (MIN_DENSITY * 10 ** (index / BUCKETS_PER_DECADE),
            MIN_DENSITY * 10 ** ((index + 1) / BUCKETS_PER_DECADE))


@dataclass
class Bucket:
    """One density bucket's decaying averages.

    Attributes:
        kind (str): `"sector"` or `"scatter"`.
        index (int): `bucket_index`.
        samples (int): Tasks counted, ever.
        seconds_per_task (float): Wall time of one task (a sector).
        seconds_per_system (float): That time over the systems it made
            (a sector with none counts as one).
        systems_per_task (float): Systems a task made.
        stars_per_system (float): Stars per system.
        max_density (float): The densest task seen in this bucket.
    """

    kind: str
    index: int
    samples: int = 0
    seconds_per_task: float = 0.0
    seconds_per_system: float = 0.0
    systems_per_task: float = 0.0
    stars_per_system: float = 0.0
    max_density: float = 0.0

    def add(self, density, seconds, systems, stars):
        """Adds one finished task to the averages."""
        per_system = seconds / max(systems, 1)
        star_ratio = stars / systems if systems else None
        if self.samples == 0:
            self.seconds_per_task = seconds
            self.seconds_per_system = per_system
            self.systems_per_task = systems
            self.stars_per_system = star_ratio or 1.0
        else:
            self.seconds_per_task += DECAY * (seconds - self.seconds_per_task)
            self.seconds_per_system += DECAY * (per_system - self.seconds_per_system)
            self.systems_per_task += DECAY * (systems - self.systems_per_task)
            if star_ratio is not None:
                self.stars_per_system += DECAY * (star_ratio - self.stars_per_system)
        self.samples += 1
        self.max_density = max(self.max_density, density or 0.0)

    def as_dict(self):
        low, high = bucket_bounds(self.index)
        return {
            "kind": self.kind, "bucket": self.index, "density_low": low, "density_high": high,
            "samples": self.samples, "seconds_per_task": self.seconds_per_task,
            "seconds_per_system": self.seconds_per_system, "systems_per_task": self.systems_per_task,
            "stars_per_system": self.stars_per_system, "max_density": self.max_density,
        }


_COLUMNS = ("samples", "seconds_per_task", "seconds_per_system", "systems_per_task", "stars_per_system",
            "max_density")


class GenerationStats:
    """
    The stored buckets and sizes, loaded once and written back by
    `flush`. Without a control database it holds nothing and writes
    nothing.

    Args:
        control_config (MySQLConfig, optional): The control database.
    """

    def __init__(self, control_config=None):
        self.control_config = control_config
        self.buckets = {}
        self.sizes = {}
        self._dirty = set()
        self._dirty_sizes = set()
        self._last_flush = time.monotonic()
        self.available = False
        if control_config is not None:
            self._load()

    def _connect(self):
        from stellarObjects import _db

        return _db.get_control_connection(self.control_config)

    def _load(self):
        try:
            conn = self._connect()
        except Exception as exc:  # noqa: BLE001 -- no control database: defaults only
            log.debug(f"Generation stats: no control database ({exc}); using defaults.")
            return
        try:
            self.read(conn)
            self.available = True
        except Exception as exc:  # noqa: BLE001 -- control schema older than v6
            log.debug(f"Generation stats: can't read them ({exc}); run update.sh. Using defaults.")
        finally:
            conn.close()

    def read(self, conn):
        """Loads every stored bucket and size from an open control
        database connection (raises when the tables aren't there)."""
        for row in conn.execute(f"SELECT kind, bucket, {', '.join(_COLUMNS)} FROM generation_stats").fetchall():
            bucket = Bucket(row["kind"], int(row["bucket"]),
                            **{name: (int(row[name]) if name == "samples" else float(row[name] or 0.0))
                               for name in _COLUMNS})
            self.buckets[(bucket.kind, bucket.index)] = bucket
        for row in conn.execute("SELECT database_name, bytes_per_system, systems, total_bytes"
                                " FROM generation_size").fetchall():
            self.sizes[row["database_name"]] = {
                "bytes_per_system": float(row["bytes_per_system"]), "systems": int(row["systems"]),
                "total_bytes": int(row["total_bytes"]),
            }
        return self

    # -- recording --------------------------------------------------------

    def record(self, kind, density, seconds, systems=0, stars=0):
        """Adds one finished task (a sector filled, a scatter layer) to
        its density's bucket, writing back every `FLUSH_SECONDS`."""
        if seconds is None or not math.isfinite(seconds) or seconds < 0:
            return
        index = bucket_index(density)
        bucket = self.buckets.setdefault((kind, index), Bucket(kind, index))
        bucket.add(density, seconds, systems, stars)
        self._dirty.add((kind, index))
        if time.monotonic() - self._last_flush >= FLUSH_SECONDS:
            self.flush()

    def measure_size(self, galaxy_conn, database):
        """Measures `database`'s bytes per star system from its own
        tables (`information_schema`), when it holds any systems, and
        keeps it for the next estimate."""
        try:
            systems = galaxy_conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"]
            total = galaxy_conn.execute(
                "SELECT COALESCE(SUM(data_length + index_length), 0) AS b FROM information_schema.tables"
                " WHERE table_schema = DATABASE()"
            ).fetchone()["b"]
        except Exception as exc:  # noqa: BLE001 -- a size we can't read just isn't updated
            log.debug(f"Generation stats: can't measure the database's size ({exc}).")
            return None
        if not systems:
            return None
        size = {"bytes_per_system": float(total) / systems, "systems": int(systems), "total_bytes": int(total)}
        self.sizes[database] = size
        self._dirty_sizes.add(database)
        return size

    def flush(self):
        """Writes every changed bucket and size back (a no-op without the
        control database)."""
        self._last_flush = time.monotonic()
        if not self.available or not (self._dirty or self._dirty_sizes):
            self._dirty.clear()
            self._dirty_sizes.clear()
            return
        try:
            conn = self._connect()
        except Exception as exc:  # noqa: BLE001
            log.debug(f"Generation stats: not saved ({exc}).")
            return
        try:
            with conn:
                for key in sorted(self._dirty):
                    b = self.buckets[key]
                    low, high = bucket_bounds(b.index)
                    conn.execute(
                        f"INSERT INTO generation_stats (kind, bucket, density_low, density_high, {', '.join(_COLUMNS)},"
                        " updated_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, NOW(6)) ON DUPLICATE KEY UPDATE "
                        + ", ".join(f"{name} = VALUES({name})" for name in _COLUMNS) + ", updated_at = NOW(6)",
                        (b.kind, b.index, low, high) + tuple(getattr(b, name) for name in _COLUMNS),
                    )
                for database in sorted(self._dirty_sizes):
                    size = self.sizes[database]
                    conn.execute(
                        "INSERT INTO generation_size (database_name, bytes_per_system, systems, total_bytes,"
                        " measured_at) VALUES (?, ?, ?, ?, NOW(6)) ON DUPLICATE KEY UPDATE"
                        " bytes_per_system = VALUES(bytes_per_system), systems = VALUES(systems),"
                        " total_bytes = VALUES(total_bytes), measured_at = NOW(6)",
                        (database, size["bytes_per_system"], size["systems"], size["total_bytes"]),
                    )
            self._dirty.clear()
            self._dirty_sizes.clear()
        except Exception as exc:  # noqa: BLE001 -- stats never stop a run
            log.debug(f"Generation stats: not saved ({exc}).")
        finally:
            conn.close()

    # -- reading ----------------------------------------------------------

    def seconds_per_system(self, kind, density):
        """The measured worker seconds per star system at `density`: its
        own bucket's, else the nearest measured bucket's, else
        `DEFAULT_SECONDS_PER_SYSTEM`."""
        index = bucket_index(density)
        measured = [b for (k, _i), b in self.buckets.items() if k == kind and b.samples]
        if not measured:
            return DEFAULT_SECONDS_PER_SYSTEM
        nearest = min(measured, key=lambda b: (abs(b.index - index), -b.samples))
        return nearest.seconds_per_system

    def bytes_per_system(self, database):
        size = self.sizes.get(database)
        return size["bytes_per_system"] if size else DEFAULT_BYTES_PER_SYSTEM

    def rows(self, kind=None):
        """Every bucket as a dict, densest last."""
        return [b.as_dict() for (k, i), b in sorted(self.buckets.items()) if kind is None or k == kind]


@dataclass
class DiskSpace:
    """The disk the galaxy database lives on."""

    path: str
    total_bytes: int
    free_bytes: int


def database_disk(galaxy_conn, mysql_host):
    """
    The disk holding the MySQL server's data directory, when that server
    runs on this machine and the directory can be looked at; `None`
    otherwise (a database on another machine can't be measured from
    here, so nothing is refused for space then).
    """
    from stellarObjects.workQueue import _is_local_host

    if not _is_local_host(mysql_host):
        return None
    try:
        path = galaxy_conn.execute("SELECT @@datadir AS d").fetchone()["d"]
        stats = os.statvfs(path)
    except Exception as exc:  # noqa: BLE001 -- e.g. no os.statvfs on Windows
        log.debug(f"Generation estimate: can't read the database's disk ({exc}).")
        try:
            import shutil

            usage = shutil.disk_usage(path)
            return DiskSpace(path, usage.total, usage.free)
        except Exception:  # noqa: BLE001
            return None
    return DiskSpace(path, stats.f_blocks * stats.f_frsize, stats.f_bavail * stats.f_frsize)


@dataclass
class Estimate:
    """
    What a bulk run is expected to make and take.

    Attributes:
        sectors (int): Sectors to fill.
        systems (float): Expected star systems.
        stars (float): Expected stars.
        bytes (int): Expected database growth, `SIZE_MARGIN` included.
        seconds (float): Expected wall time with `workers`.
        workers (int): Worker processes.
        disk (DiskSpace or None): The database's disk, when measurable.
        measured (bool): Whether the speed came from this server's own
            stored runs (otherwise defaults).
    """

    sectors: int = 0
    systems: float = 0.0
    stars: float = 0.0
    bytes: int = 0
    seconds: float = 0.0
    workers: int = 1
    disk: DiskSpace = None
    measured: bool = False
    refusal: str = field(default=None)

    def as_dict(self):
        return {
            "sectors": self.sectors, "systems": round(self.systems), "stars": round(self.stars),
            "bytes": self.bytes, "seconds": round(self.seconds, 1), "workers": self.workers,
            "measured": self.measured, "refused": self.refusal is not None, "refusal": self.refusal,
            "disk": None if self.disk is None else {
                "path": self.disk.path, "total_bytes": self.disk.total_bytes, "free_bytes": self.disk.free_bytes,
            },
            "summary": self.summary(),
        }

    def summary(self):
        """One line for a person."""
        text = (f"About {format_bytes(self.bytes)} and {format_duration(self.seconds)} for {self.sectors:,} "
                f"sector{'s' if self.sectors != 1 else ''} (~{round(self.systems):,} star systems, "
                f"~{round(self.stars):,} stars, {self.workers} worker{'s' if self.workers != 1 else ''}"
                f"{'' if self.measured else ', speed not yet measured on this server'}).")
        if self.disk is not None:
            text += f" {format_bytes(self.disk.free_bytes)} free of {format_bytes(self.disk.total_bytes)}."
        return text


MAX_SYSTEMS_PER_SECTOR = 1e9
"""float: The most systems one sector counts for in an estimate (an
absurd or infinite `--density` fills its cube and stops well short)."""


def estimate(sectors, stats, database, workers=1, kind="sector"):
    """
    The size and time `sectors` will take, and whether it's refused.

    Args:
        sectors (iterable): `(density, expected_systems)` per sector to
            fill.
        stats (GenerationStats): The stored speeds and sizes.
        database (str): The galaxy database's name (its size per system).
        workers (int): Worker processes the run will use.
        kind (str): Which speeds to use.

    Returns:
        Estimate: `disk` and `refusal` still unset (`check_disk`).
    """
    result = Estimate(workers=max(1, int(workers or 1)))
    worker_seconds = 0.0
    stars_ratio = _stars_per_system(stats, kind)
    for density, expected in sectors:
        expected = min(max(float(expected or 0.0), 0.0), MAX_SYSTEMS_PER_SECTOR)
        result.sectors += 1
        result.systems += expected
        worker_seconds += stats.seconds_per_system(kind, density) * max(expected, 1.0)
    result.stars = result.systems * stars_ratio
    result.bytes = int(math.ceil(round(result.systems * stats.bytes_per_system(database) * (1 + SIZE_MARGIN), 6)))
    result.seconds = worker_seconds / min(result.workers, max(result.sectors, 1))
    result.measured = any(b.samples for (k, _i), b in stats.buckets.items() if k == kind)
    return result


def _stars_per_system(stats, kind):
    measured = [b for (k, _i), b in stats.buckets.items() if k == kind and b.samples and b.stars_per_system]
    if not measured:
        return 1.3
    return sum(b.stars_per_system * b.samples for b in measured) / sum(b.samples for b in measured)


def check_disk(result, disk):
    """
    Sets `result.disk` and, when the run would take more than
    `MAX_DISK_SHARE` of the disk or leave less than `MIN_FREE_BYTES`
    free, `result.refusal` (why, and how much it needs).

    Returns:
        str or None: The refusal.
    """
    result.disk = disk
    result.refusal = None
    if disk is None:
        return None
    if result.bytes > disk.total_bytes * MAX_DISK_SHARE:
        result.refusal = (
            f"Refused: this would take about {format_bytes(result.bytes)}, more than a quarter of the "
            f"database disk ({format_bytes(disk.total_bytes)} at {disk.path}). Generate fewer sectors "
            f"(at most about {format_bytes(int(disk.total_bytes * MAX_DISK_SHARE))} worth)."
        )
    elif disk.free_bytes - result.bytes < MIN_FREE_BYTES:
        result.refusal = (
            f"Refused: this would take about {format_bytes(result.bytes)} and leave "
            f"{format_bytes(max(disk.free_bytes - result.bytes, 0))} free on the database disk ({disk.path}); "
            f"at least {format_bytes(MIN_FREE_BYTES)} must stay free. It needs "
            f"{format_bytes(result.bytes + MIN_FREE_BYTES - disk.free_bytes)} more free space, or fewer sectors."
        )
    return result.refusal


def format_bytes(count):
    """`"1.4 GB"`, `"820 MB"`, `"12 KB"` (powers of 1000)."""
    count = float(count or 0)
    if not math.isfinite(count):
        return "more than any disk holds"
    for unit, scale in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("KB", 1e3)):
        if count >= scale:
            value = count / scale
            return f"{value:.1f} {unit}" if value < 10 else f"{value:,.0f} {unit}"
    return f"{int(count)} bytes"


def format_duration(seconds):
    """`"3 days 4 h"`, `"2 h 05 m"`, `"4 m 10 s"`, `"12 s"`."""
    if seconds is not None and not math.isfinite(seconds):
        return "longer than anyone will wait"
    seconds = max(int(round(seconds or 0)), 0)
    days, rest = divmod(seconds, 86400)
    hours, rest = divmod(rest, 3600)
    minutes, secs = divmod(rest, 60)
    if days:
        return f"{days} day{'s' if days != 1 else ''} {hours} h"
    if hours:
        return f"{hours} h {minutes:02d} m"
    if minutes:
        return f"{minutes} m {secs:02d} s"
    return f"{secs} s"
