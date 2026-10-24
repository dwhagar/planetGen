# planetgen/util/memwatch.py

"""
Peak memory of a job step (PERF.183)
====================================

Watches the processes a job step runs beside, every `INTERVAL_SECONDS`, and keeps the highest figure of each group:

* `mysql`: every `mysqld` on the host;
* `workers`: every planetgen worker process (`planetgen.cli.worker`);
* `step`: the step's own process tree (its process group, which the job runner makes for it);
* `host`: the whole host's memory in use.

The figures are proportional set sizes (`Pss` from `/proc/<pid>/smaps_rollup`) where the system has them, so a work
horse forked from its worker is not counted for the worker's pages as well, else resident set sizes. The job runner
prints `PeakTracker.line` after each step in the stage log, and `planetgen.cli.smoke` compares the figures with the
budget in `planetgen.tuning`.

Standard library only, so the job runner can import it before anything else.
"""

import os
import subprocess
import threading
import time

INTERVAL_SECONDS = 2.0
"""float: How often the processes are looked at."""

WORKER_COMMAND = "planetgen.cli.worker"
"""str: What every planetgen worker process has in its command line (the same as `planetgen.queue.workers`)."""

GROUPS = ("mysql", "workers", "step", "host")
"""tuple: The groups a peak is kept for."""


def _megabytes_of(pid):
    """The memory of process `pid` in megabytes (Pss, else Rss), or `None` when it cannot be read."""
    try:
        with open(f"/proc/{pid}/smaps_rollup", "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("Pss:"):
                    return int(line.split()[1]) / 1024.0
    except (OSError, ValueError, IndexError):
        pass
    return None


def _process_list():
    """`[(pid, ppid, pgid, rss_kb, args)]` of every process, or `[]` where `ps` is missing or fails."""
    try:
        listing = subprocess.run(["ps", "-axww", "-o", "pid=,ppid=,pgid=,rss=,args="], capture_output=True, text=True,
                                 timeout=10, check=False).stdout
    except (OSError, subprocess.SubprocessError):
        return []
    rows = []
    for line in listing.splitlines():
        parts = line.split(None, 4)
        if len(parts) == 5 and all(part.isdigit() for part in parts[:4]):
            rows.append((int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), parts[4]))
    return rows


def host_memory_mb():
    """`(used, total)` of the host's memory in megabytes (`MemTotal` less `MemAvailable`), or `(None, None)`."""
    try:
        found = {}
        with open("/proc/meminfo", "r", encoding="utf-8") as handle:
            for line in handle:
                name, _sep, rest = line.partition(":")
                if name in ("MemTotal", "MemAvailable"):
                    found[name] = int(rest.split()[0]) / 1024.0
        return found["MemTotal"] - found["MemAvailable"], found["MemTotal"]
    except (OSError, ValueError, KeyError, IndexError):
        return None, None


def sample(step_group=None):
    """
    One look: `{"mysql": MB, "workers": MB, "workers_n": count, "step": MB, "host": MB}`. `step_group` is the process
    group of the step (`None` for none). A figure that cannot be read is 0.
    """
    figures = {"mysql": 0.0, "workers": 0.0, "workers_n": 0, "step": 0.0, "host": 0.0}
    for pid, _ppid, pgid, rss_kb, args in _process_list():
        command = args.split(None, 1)[0] if args else ""
        in_mysql = os.path.basename(command) == "mysqld"
        in_workers = WORKER_COMMAND in args
        in_step = step_group is not None and pgid == step_group
        if not (in_mysql or in_workers or in_step):
            continue
        megabytes = _megabytes_of(pid)
        if megabytes is None:
            megabytes = rss_kb / 1024.0
        if in_mysql:
            figures["mysql"] += megabytes
        if in_workers:
            figures["workers"] += megabytes
            figures["workers_n"] += 1
        if in_step:
            figures["step"] += megabytes
    used, _total = host_memory_mb()
    figures["host"] = used or 0.0
    return figures


class PeakTracker:
    """
    Keeps the highest figure of each group while a step runs. Use it as a context manager, with the step's process
    group once it is known (`set_group`); `peaks` is `{group: MB, "workers_n": most workers seen, "seconds": the
    time it watched, "host_total": the host's memory}`.
    """

    def __init__(self, interval=None, clock=None, sampler=None):
        self.interval = INTERVAL_SECONDS if interval is None else interval
        self._clock = clock or time.monotonic
        self._sampler = sampler or sample
        self._group = None
        self._stop = threading.Event()
        self._thread = None
        self._started = None
        self.peaks = {"mysql": 0.0, "workers": 0.0, "workers_n": 0, "step": 0.0, "host": 0.0, "seconds": 0.0,
                      "host_total": None}

    def set_group(self, group):
        """Tells it which process group is the step's."""
        self._group = group

    def _take(self):
        figures = self._sampler(self._group)
        for name in ("mysql", "workers", "workers_n", "step", "host"):
            self.peaks[name] = max(self.peaks[name], figures.get(name, 0))

    def _run(self):
        while not self._stop.is_set():
            try:
                self._take()
            except Exception:  # noqa: BLE001 -- watching must never stop the step
                pass
            self._stop.wait(self.interval)

    def __enter__(self):
        self._started = self._clock()
        self.peaks["host_total"] = host_memory_mb()[1]
        self._thread = threading.Thread(target=self._run, name="planetgen-memwatch", daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *_exc):
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=self.interval + 15)
        try:
            self._take()
        except Exception:  # noqa: BLE001
            pass
        self.peaks["seconds"] = self._clock() - self._started
        return False

    def line(self):
        """The line the stage log gets after the step: what peaked, in whole megabytes."""
        return format_line(self.peaks)


def format_line(peaks):
    """`"Memory peak: MySQL 884 MB, 6 workers 2,310 MB, the step itself 256 MB, the whole host 3,100 of 15,990 MB."`"""
    total = peaks.get("host_total")
    host = f"{peaks['host']:,.0f} of {total:,.0f} MB" if total else f"{peaks['host']:,.0f} MB"
    workers = peaks["workers_n"]
    return (f"Memory peak: MySQL {peaks['mysql']:,.0f} MB, {workers} {'worker' if workers == 1 else 'workers'} "
            f"{peaks['workers']:,.0f} MB, the step itself {peaks['step']:,.0f} MB, the whole host {host}.")
