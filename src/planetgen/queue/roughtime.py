# planetgen/queue/roughtime.py

"""
A rough time for work nothing has timed yet (PERF.184). A server that has never run a stage has no recorded time, no
rate and no unit that has finished, and the estimate used to say "estimating" until the first unit came back. This is
the guess `runeta` falls back to then: what the host can run at once (its cores and memory), the units the stage
queued, and a default cost per unit. It is labelled as rough, and it is replaced as soon as there is something better:
the run's own rate, then the time and rate recorded for the stage (`stagetimes`).

The two constants are stand-ins (tier: educated guess from needed output, tweak level L3): a figure that is
sensible for the stages that exist, not a measurement. Pure functions of their arguments; `host_workers` reads the host.
"""

import math

UNIT_SECONDS = 20.0
"""float: Seconds one worker is assumed to spend on one unit of a stage nothing has timed."""

UNKNOWN_STAGE_WORKER_SECONDS = 600.0
"""float: Worker seconds assumed for a stage whose units are not known yet (not begun), divided by the workers."""

ROUGH_LOW_FACTOR = 0.4
ROUGH_HIGH_FACTOR = 3.0
"""float: The range around a rough figure, in times the figure: a guess with no measurement behind it is wide."""

OVERRUN_KEEP = 0.1
"""float: A stage that outlasts its rough guess keeps this share of the guess as its time left."""

CORE_SHARE = 0.8
"""float: The share of the cores the default worker count uses (`generation.run_common` sizes a run the same way)."""


def host_workers(cores=None, memory_bytes=None):
    """
    How many workers a run on this host would use: `CORE_SHARE` of the cores, no more than the memory holds at
    `queue.workers.WORKER_MEMORY_MB` each within `HOST_MEMORY_SHARE`, at least 1.

    Args:
        cores (int, optional): The host's cores; default the cores this process may run on.
        memory_bytes (int, optional): The host's memory; default what the OS says (`None`: no memory limit).
    """
    from planetgen.queue import work, workers
    cores = work.cpu_count() if cores is None else max(1, int(cores))
    if memory_bytes is None:
        memory_bytes = workers._physical_memory_bytes()
    count = max(1, int(cores * CORE_SHARE))
    if memory_bytes:
        count = min(count, max(1, int(memory_bytes * workers.HOST_MEMORY_SHARE / (workers.WORKER_MEMORY_MB * 2 ** 20))))
    return count


def stage_seconds(units=None, workers=None):
    """
    The rough seconds a stage takes: `units` units of `UNIT_SECONDS` each over the workers that can work at once
    (`min(workers, units)`), or `UNKNOWN_STAGE_WORKER_SECONDS` over the workers when `units` is not known yet.

    Args:
        units (int, optional): The stage's units, `None` before it is expanded.
        workers (int, optional): The workers the stage runs on; default `host_workers()`.
    """
    workers = max(1, int(workers or host_workers()))
    if not units:
        return UNKNOWN_STAGE_WORKER_SECONDS / workers
    return math.ceil(units / min(workers, units)) * UNIT_SECONDS


def rough_range(seconds):
    """`(low, high)` around a rough figure (`ROUGH_LOW_FACTOR`, `ROUGH_HIGH_FACTOR`), `None` for no figure."""
    if seconds is None:
        return None
    return seconds * ROUGH_LOW_FACTOR, seconds * ROUGH_HIGH_FACTOR
