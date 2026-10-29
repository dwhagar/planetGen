# planetgen/queue/work.py

"""
The machine's worker count and a task's seed.

The units of work themselves are RQ jobs on Redis (`planetgen.queue.units`, `planetgen.queue.wrapper`) and a run is
recorded there (`planetgen.queue.runrecord`); this file keeps `worker_count` (at most about 80% of the machine's
cores) and `task_seed` (a task's `random` seed from the run's seed and its own key, so what it generates doesn't
depend on how many workers there are or which finished first). The MySQL job tree and the `work_lease` it shared
with the queue are gone (PERF.98, PERF.100).
"""

import hashlib
import ipaddress
import math
import os
import socket

CPU_SHARE = 0.8
"""float: The share of the machine's cores the pool may use."""

WORKERS_ENV_VAR = "PLANETGEN_WORKERS"
"""str: Environment variable overriding the worker count (`--workers`
on the command line wins over it); 0 or "auto" means `CPU_SHARE` of the
cores."""

WORKER_NICENESS = 10
"""int: How much `os.nice` lowers a worker's priority on Linux and
macOS."""


def cpu_count():
    """The cores this process may run on (its affinity mask where the OS
    has one), at least 1."""
    if hasattr(os, "process_cpu_count"):
        count = os.process_cpu_count()
    elif hasattr(os, "sched_getaffinity"):
        count = len(os.sched_getaffinity(0))
    else:
        count = os.cpu_count()
    return max(1, count or 1)


def _is_local_host(host):
    if not host:
        return False
    if host.lower() in ("localhost", socket.gethostname().lower()):
        return True
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return False


def worker_count(requested=None, mysql_host=None):
    """
    How many worker processes a run uses.

    Args:
        requested (int | str | None): `--workers`, else the
            `PLANETGEN_WORKERS` environment variable; a positive number
            is used as given, and 0, "auto" or nothing means `CPU_SHARE`
            of `cpu_count()`.
        mysql_host (str, optional): The database host. When it's this
            machine, one fewer worker, so MySQL keeps a core of the 80%.

    Returns:
        int: At least 1.
    """
    if requested is None:
        requested = os.environ.get(WORKERS_ENV_VAR)
    if isinstance(requested, str):
        requested = requested.strip().lower()
        requested = 0 if requested in ("", "auto") else int(requested)
    if requested:
        return max(1, int(requested))
    count = max(1, math.floor(CPU_SHARE * cpu_count()))
    if count > 1 and _is_local_host(mysql_host):
        count -= 1
    return count


def lower_priority(niceness=WORKER_NICENESS):
    """Lowers this process's scheduling priority (by `WORKER_NICENESS`, or by `niceness`) so everything else on
    the machine comes first. Never raises."""
    try:
        os.nice(niceness)
    except Exception:  # noqa: BLE001 -- a normal-priority worker is still a worker
        pass


def task_seed(run_seed, key):
    """The `random` seed for one task: the run's seed and the task's key,
    hashed, so it doesn't depend on worker count or finish order."""
    digest = hashlib.sha256(f"{run_seed}:{key}".encode("utf-8")).digest()
    return int.from_bytes(digest[:16], "big")
