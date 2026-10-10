# planetgen/generation/stats.py

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

Size is measured from the galaxy database itself (the data and index
bytes of the tables a fill writes, over its star systems), after every
run, per galaxy database. The galaxy-wide tables a plan or the
bright-star scatter writes (`GALAXY_WIDE_TABLES`) are left out: they
don't grow with the systems a run adds, and in a freshly planned galaxy
they would make the first few systems look enormous (PERF.26).

Everything here fails open: without the control database (not created
yet, or no grant on it) the estimate uses `DEFAULT_SECONDS_PER_SYSTEM`
and `DEFAULT_BYTES_PER_SYSTEM`, and nothing is recorded.
"""

import math
import os
import shutil
import time
from dataclasses import dataclass, field

from planetgen.galaxy.version_key import version_key as current_version_key
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

DEFAULT_BYTES_PER_SYSTEM = 48 * 1024
"""int: Database bytes per star system before the galaxy has any.
Measured 2026-10-07 (PERF.26) on MariaDB 10.11 as the growth of every
table across a fill: 52 KB each in a sparse ring (550 systems), 48 KB
in a mid-disk ring (2,288), 37 KB in a dense inner ring (3,120), 44 KB
in the test suite's small fill. Moons are about 40%, then planets,
`nearest_systems` and rogue planets. The old 70 KB (2026-10-01) counted
the whole database, galaxy-wide tables and empty ones included."""

GALAXY_WIDE_TABLES = frozenset({
    "bright_stars", "galaxy_column", "galaxy_layer", "galaxy_shape", "generation_runs", "id_blocks",
    "schema_migrations", "sector_stats",
})
"""frozenset: Tables a plan or the bright-star scatter fills for the
whole galaxy at once, left out of the size per system."""

STATS_ENV_VAR = "PLANETGEN_GENERATION_STATS"
"""str: `0` keeps `planetgen` from reading or recording any stats (the
test suite sets it, so test runs never touch a real server's numbers)."""

FLUSH_SECONDS = 30
"""int: How often a run writes its new averages back while it runs."""

KINDS = ("sector", "scatter", "phenomena")
"""tuple: What a bucket times: a sector's fill, (PERF.9) a layer of the
bright-star scatter by the layer's mean density, or (PERF.32) a layer of
the phenomenon scatter. A benchmark's kinds start with `BENCH_PREFIX`."""

BENCH_PREFIX = "bench:"
"""str: Starts the kind of every row a benchmark (PERF.31) writes, so none
feeds a live estimate; `estimate` and `seconds_per_system` only read the
kind they are asked for."""


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
        kind (str): One of `KINDS`.
        index (int): `bucket_index`.
        workers (int): The worker count of the runs it averages (rates
            at four workers are not rates at one).
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
    workers: int = 1
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
            "kind": self.kind, "workers": self.workers, "bucket": self.index,
            "density_low": low, "density_high": high,
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
        """dict: `(kind, workers, bucket index)` -> `Bucket`."""
        self.sizes = {}
        self._dirty = set()
        self._dirty_sizes = set()
        self._purged = False
        self._last_flush = time.monotonic()
        self.available = False
        if control_config is not None:
            self._load()

    def _connect(self):
        from planetgen.db import store

        return store.get_control_connection(self.control_config)

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
        """Loads the stored buckets of this version (`version_key`, PERF.32;
        rows of any other are left alone here and deleted by the first
        `flush`) and every size from an open control database connection
        (raises when the tables aren't there)."""
        for row in conn.execute(
                f"SELECT kind, workers, bucket, {', '.join(_COLUMNS)} FROM generation_stats WHERE version_key = ?",
                (current_version_key(),)).fetchall():
            bucket = Bucket(row["kind"], int(row["bucket"]), workers=int(row["workers"]),
                            **{name: (int(row[name]) if name == "samples" else float(row[name] or 0.0))
                               for name in _COLUMNS})
            self.buckets[(bucket.kind, bucket.workers, bucket.index)] = bucket
        for row in conn.execute("SELECT database_name, bytes_per_system, systems, total_bytes"
                                " FROM generation_size").fetchall():
            self.sizes[row["database_name"]] = {
                "bytes_per_system": float(row["bytes_per_system"]), "systems": int(row["systems"]),
                "total_bytes": int(row["total_bytes"]),
            }
        return self

    # -- recording --------------------------------------------------------

    def record(self, kind, density, seconds, systems=0, stars=0, workers=1):
        """Adds one finished task (a sector filled, a scatter layer) to
        its density's bucket for the run's `workers`, writing back every
        `FLUSH_SECONDS`."""
        if seconds is None or not math.isfinite(seconds) or seconds < 0:
            return
        workers = max(1, int(workers or 1))
        index = bucket_index(density)
        bucket = self.buckets.setdefault((kind, workers, index), Bucket(kind, index, workers=workers))
        bucket.add(density, seconds, systems, stars)
        self._dirty.add((kind, workers, index))
        if time.monotonic() - self._last_flush >= FLUSH_SECONDS:
            self.flush()

    def measure_size(self, galaxy_conn, database):
        """Measures `database`'s bytes per star system from its own
        tables (`information_schema`, `GALAXY_WIDE_TABLES` left out),
        when it holds any systems, and keeps it for the next estimate."""
        try:
            # MySQL 8 caches these sizes for a day by default, so a run's
            # growth wouldn't show; MariaDB (no such variable) reads live.
            galaxy_conn.execute("SET SESSION information_schema_stats_expiry = 0")
        except Exception:  # noqa: BLE001
            pass
        excluded = sorted(GALAXY_WIDE_TABLES)
        try:
            systems = galaxy_conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"]
            total = galaxy_conn.execute(
                "SELECT COALESCE(SUM(data_length + index_length), 0) AS b FROM information_schema.tables"
                f" WHERE table_schema = DATABASE() AND table_name NOT IN ({', '.join(['?'] * len(excluded))})",
                excluded,
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
                self._delete_other_versions(conn)
                for key in sorted(self._dirty):
                    b = self.buckets[key]
                    low, high = bucket_bounds(b.index)
                    conn.execute(
                        "REPLACE INTO generation_stats (kind, workers, bucket, version_key, density_low, density_high,"
                        f" {', '.join(_COLUMNS)}, updated_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, NOW(6))",
                        (b.kind, b.workers, b.index, current_version_key(), low, high)
                        + tuple(getattr(b, name) for name in _COLUMNS),
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

    def _delete_other_versions(self, conn):
        """PERF.32: a rate measured by another release, Python, OS or
        architecture says nothing about this one, so the first write of a
        new version deletes them (once per run), logging how many."""
        if self._purged:
            return
        removed = conn.execute("DELETE FROM generation_stats WHERE version_key <> ?",
                               (current_version_key(),)).rowcount
        self._purged = True
        if removed:
            log.debug(f"Generation stats: deleted {removed} rate(s) measured by another version.")

    def reset(self, conn):
        """Deletes every stored rate (the Admin Reset stats button) through
        the open control connection `conn`, and forgets the loaded ones.

        Returns:
            int: Rows deleted.
        """
        self.buckets.clear()
        self._dirty.clear()
        with conn:
            return conn.execute("DELETE FROM generation_stats").rowcount

    # -- reading ----------------------------------------------------------

    def seconds_per_system(self, kind, density, workers=1):
        """The measured worker seconds per star system at `density` for a
        run of `workers`: the nearest measured density bucket among the
        rows for that worker count; with none for it, the rows of the
        neighbouring worker counts blended by distance (the nearest one
        when only one side has any); else `DEFAULT_SECONDS_PER_SYSTEM`."""
        workers = max(1, int(workers or 1))
        counts = sorted({w for (k, w, _i), b in self.buckets.items() if k == kind and b.samples})
        if not counts:
            return DEFAULT_SECONDS_PER_SYSTEM
        if workers in counts:
            return self._nearest_bucket(kind, density, workers)
        below = [w for w in counts if w < workers]
        above = [w for w in counts if w > workers]
        if below and above:
            low, high = below[-1], above[0]
            fraction = (workers - low) / (high - low)
            return ((1 - fraction) * self._nearest_bucket(kind, density, low)
                    + fraction * self._nearest_bucket(kind, density, high))
        return self._nearest_bucket(kind, density, (below or above)[-1 if below else 0])

    def _mean_for_workers(self, kind, workers, attribute):
        """`attribute` of the `kind` buckets at `workers`, averaged by their samples; with no rows for that
        worker count the neighbouring counts' averages blended by distance (the nearest when one side has
        none); `None` with nothing recorded."""
        workers = max(1, int(workers or 1))
        counts = sorted({w for (k, w, _i), b in self.buckets.items() if k == kind and b.samples})
        if not counts:
            return None

        def mean(count):
            rows = [b for (k, w, _i), b in self.buckets.items() if k == kind and w == count and b.samples]
            return sum(getattr(b, attribute) * b.samples for b in rows) / sum(b.samples for b in rows)

        if workers in counts:
            return mean(workers)
        below = [w for w in counts if w < workers]
        above = [w for w in counts if w > workers]
        if below and above:
            fraction = (workers - below[-1]) / (above[0] - below[-1])
            return (1 - fraction) * mean(below[-1]) + fraction * mean(above[0])
        return mean((below or above)[-1 if below else 0])

    def pool_rate(self, kind, workers=1):
        """
        PERF.33: the units a second the whole pool of `workers` got through
        for `kind` (a unit is what `record` counted as `systems`: a system,
        or a star of a scatter layer), averaged over densities; `None` with
        nothing recorded.
        """
        seconds = self._mean_for_workers(kind, workers, "seconds_per_system")
        return max(1, int(workers or 1)) / seconds if seconds else None

    def seconds_per_task(self, kind, workers=1):
        """PERF.33: the wall seconds one task of `kind` takes a worker when `workers` run side by side,
        averaged over densities; `None` with nothing recorded."""
        return self._mean_for_workers(kind, workers, "seconds_per_task") or None

    def _nearest_bucket(self, kind, density, workers):
        index = bucket_index(density)
        measured = [b for (k, w, _i), b in self.buckets.items() if k == kind and w == workers and b.samples]
        nearest = min(measured, key=lambda b: (abs(b.index - index), -b.samples))
        return nearest.seconds_per_system

    def bytes_per_system(self, database):
        size = self.sizes.get(database)
        return size["bytes_per_system"] if size else DEFAULT_BYTES_PER_SYSTEM

    def rows(self, kind=None):
        """Every bucket as a dict, by kind, worker count, density."""
        return [b.as_dict() for (k, _w, _i), b in sorted(self.buckets.items()) if kind is None or k == kind]


@dataclass
class DiskSpace:
    """The drive the galaxy database's data directory lives on."""

    path: str
    total_bytes: int
    free_bytes: int
    mount: str = None
    """str or None: The mount point holding `path`."""
    source: str = "data directory"
    """str: How it was measured: the server's `data directory` seen from here, or `server report`
    (MariaDB's `information_schema.DISKS`, for a database on another machine)."""

    def where(self):
        """`"/mnt/data (data directory /mnt/data/mysql/)"`: the drive and the path measured."""
        if self.mount and self.mount != self.path:
            return f"{self.mount} (data directory {self.path})"
        return self.path


@dataclass
class Unmeasured:
    """The drive could not be measured, and why (shown instead of any number)."""

    reason: str


def _mount_point(path):
    """The mount point holding `path`, symlinks and bind mounts resolved; `None` if unknown."""
    try:
        real = os.path.realpath(path)
        while not os.path.ismount(real):
            parent = os.path.dirname(real)
            if parent == real:
                break
            real = parent
        return real
    except OSError:
        return None


def _disk_from_server(conn, datadir):
    """MariaDB lists its mounts in `information_schema.DISKS` (the `disks` plugin); the one holding
    `datadir` is the longest `Path` that prefixes it. `None` when the server has no such table."""
    try:
        rows = conn.execute("SELECT Path AS p, Total AS t, Available AS a FROM information_schema.DISKS").fetchall()
    except Exception:  # noqa: BLE001 -- MySQL, or the plugin is off
        return None
    wanted = datadir.rstrip("/")
    best = None
    for row in rows:
        path = (row["p"] or "").rstrip("/")
        if (wanted == path or wanted.startswith(path + "/") or path == "") and (best is None or len(path) > len(best[0])):
            best = (path, row)
    if best is None:
        return None
    row = best[1]
    return DiskSpace(datadir, int(row["t"]), int(row["a"]), mount=row["p"] or "/", source="server report")


def database_disk(galaxy_conn, mysql_host, database=None):
    """
    The drive holding the MySQL server's data directory, asked of the
    server itself (`SELECT @@datadir`), never the boot drive: the
    directory is resolved to its mount (symlinks, bind mounts) and measured there. A server on this machine is
    measured at that path. A server on another machine is measured only if
    that path is really its data directory seen from here (the database's
    own folder is in it) or, failing that, from the server's own report
    (MariaDB's `information_schema.DISKS`).

    Args:
        galaxy_conn (Connection): Any connection to the server.
        mysql_host (str): Its host, as configured.
        database (str, optional): The schema, to check a path on another
            machine is the same data directory.

    Returns:
        DiskSpace or Unmeasured: What to show; an `Unmeasured` says why
            (and nothing is refused for space then).
    """
    from planetgen.queue.work import _is_local_host

    try:
        datadir = galaxy_conn.execute("SELECT @@datadir AS d").fetchone()["d"]
    except Exception as exc:  # noqa: BLE001
        log.debug(f"Database disk: the server did not report its data directory ({exc}).")
        return Unmeasured("the server did not report its data directory")
    local = _is_local_host(mysql_host)
    visible = local or (bool(database) and os.path.isdir(os.path.join(datadir, database)))
    if visible:
        try:
            usage = shutil.disk_usage(datadir)
            disk = DiskSpace(datadir, usage.total, usage.free, mount=_mount_point(datadir))
            log.debug(f"Database disk: {disk.where()}, {usage.free} of {usage.total} bytes free.")
            return disk
        except OSError as exc:
            log.debug(f"Database disk: can't read {datadir} ({exc}).")
    reported = _disk_from_server(galaxy_conn, datadir)
    if reported is not None:
        log.debug(f"Database disk: {reported.where()} from the server's own report.")
        return reported
    where = "its data directory isn't visible from this machine" if not local else f"can't read {datadir}"
    return Unmeasured(f"the database server's data directory ({datadir}) is not measurable here: {where}")


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
    disk_note: str = None

    def as_dict(self):
        return {
            "sectors": self.sectors, "systems": round(self.systems), "stars": round(self.stars),
            "bytes": self.bytes, "seconds": round(self.seconds, 1), "workers": self.workers,
            "measured": self.measured, "refused": self.refusal is not None, "refusal": self.refusal,
            "disk": None if self.disk is None else {
                "path": self.disk.path, "mount": self.disk.mount, "where": self.disk.where(),
                "total_bytes": self.disk.total_bytes, "free_bytes": self.disk.free_bytes,
            },
            "disk_note": self.disk_note,
            "summary": self.summary(),
        }

    def summary(self):
        """One line for a person."""
        text = (f"About {format_bytes(self.bytes)} and {format_duration(self.seconds)} for {self.sectors:,} "
                f"sector{'s' if self.sectors != 1 else ''} (~{round(self.systems):,} star systems, "
                f"~{round(self.stars):,} stars, {self.workers} worker{'s' if self.workers != 1 else ''}"
                f"{'' if self.measured else ', speed not yet measured on this server'}).")
        if self.disk is not None:
            text += (f" {format_bytes(self.disk.free_bytes)} free of {format_bytes(self.disk.total_bytes)} "
                     f"on {self.disk.where()}.")
        elif self.disk_note:
            text += f" Free space not measured: {self.disk_note}; nothing is refused for space."
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
        worker_seconds += stats.seconds_per_system(kind, density, result.workers) * max(expected, 1.0)
    result.stars = result.systems * stars_ratio
    result.bytes = int(math.ceil(round(result.systems * stats.bytes_per_system(database) * (1 + SIZE_MARGIN), 6)))
    result.seconds = worker_seconds / min(result.workers, max(result.sectors, 1))
    result.measured = any(b.samples for (k, _w, _i), b in stats.buckets.items() if k == kind)
    return result


def _stars_per_system(stats, kind):
    measured = [b for (k, _w, _i), b in stats.buckets.items() if k == kind and b.samples and b.stars_per_system]
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
    result.refusal = None
    result.disk_note = None
    if isinstance(disk, Unmeasured):
        result.disk, result.disk_note = None, disk.reason
        return None
    result.disk = disk
    if disk is None:
        return None
    if result.bytes > disk.total_bytes * MAX_DISK_SHARE:
        result.refusal = (
            f"Refused: this would take about {format_bytes(result.bytes)}, more than a quarter of the "
            f"database disk ({format_bytes(disk.total_bytes)} at {disk.where()}). Generate fewer sectors "
            f"(at most about {format_bytes(int(disk.total_bytes * MAX_DISK_SHARE))} worth)."
        )
    elif disk.free_bytes - result.bytes < MIN_FREE_BYTES:
        result.refusal = (
            f"Refused: this would take about {format_bytes(result.bytes)} and leave "
            f"{format_bytes(max(disk.free_bytes - result.bytes, 0))} free on the database disk ({disk.where()}); "
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
