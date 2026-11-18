# planetgen/queue/workers.py

"""
The workers of the shared queues (PERF.95, folding in PERF.40, PERF.78 and ADM.15): who starts them, how
many there are, and how a running job is stopped. Section 5.4 of docs/design/operations/redis-and-the-work-queue.md.

* Lanes (PERF.151): `redisqueue.BULK_QUEUE` (a run's units) and the maintenance lane. The site's short and long jobs are light jobs the supervisor forks (`queue/lightjobs.py`, OPS.230).
* OPS.235: RQ is gone. Nothing queues a job on Redis for a worker to take; the supervisor's launcher forks each unit, and a page job is a ledger run (`queue/pagerun.py`). What is left here is what
  sizes them: the memory one worker of a stage takes (`stage_worker_mb`) and the most a stage may run at once (`stage_worker_cap`).
"""

import os

WORKER_MEMORY_MB = 300
"""int: What one worker takes (its libraries, sampling tables and the work horse it forks), for `roughtime`'s estimate of the workers a host holds."""

HOST_MEMORY_SHARE = 0.6
"""float: The share of the host's memory the workers may take at most; the rest is for MySQL, Redis and the web server."""

BULK_WORKER_MEMORY_MB = 2100
"""int: What one worker takes at most in a stage that draws stars from the sampling tables but is not a sector stage (PERF.226): the
scatters. Raised from 1,700 to 2,100 MB (PERF.252): the bright stars of the full-galaxy run of 2026-11-08 held 1,879 MB, and 10 percent more is the margin. Such a worker holds the tables of every star population (about 1.4 GB with the datasets they read), and a worker and its
horse measured 1,456 to 1,564 MB together (`units.peak_resident_mb`, shared pages counted once) before the tables moved out of the
workers of the other stages. The kernel killed the run of 2026-10-31 (36 workers, 16,285 MB on a 15,992 MB host, MySQL first at
3:22 AM EDT) with the scatters sized by a guess of 300 MB."""

LEAN_WORKER_MEMORY_MB = 700
"""int: What one worker takes at most in a stage that draws no star (PERF.226): a worker starts without the sampling tables
(`cli/worker.warm`) and starts over to shed them before such a unit. Measured on a populated database, the worker and its horse
together: 601 MB for a link, 588 for a settle, 552 for a habitability re-score and 540 for a map-luminosity fill (1,456 to 1,564 MB
while every worker held the tables), 700 MB with the margin."""

SCATTER_WORKER_MEMORY_MB = BULK_WORKER_MEMORY_MB
"""int: What one worker takes at most while it scatters stars (PERF.220, PERF.226): a table worker's (`BULK_WORKER_MEMORY_MB`)."""

SKY_SHARE_WORKER_MEMORY_MB = 1650
"""int: What one unit of the bright-star backfill (`scatter-backfill`: a share of the sky tiers, or a run of sectors) takes at most (OPS.270). Measured 2026-11-18 on a share forked from a warm template whose sampling tables were read
back from the band-table cache (`band_table_cache`): the template holds 175 MB and the share's own private memory climbs to 1,502 MB in the first 150 s and stays there (6.5 minutes watched), because the
stellar-evolution caches of the share (`main_sequence_bright_from_gy`, `_evolved_bright_windows`, `timeline`, `branch`, `_reach`; each 130,000 to 520,000 entries) fill to their size. The sampling tables are not the cost (391 tables
are 7 MB). 1,502 MB with 10 percent more for the margin is 1,650; the figure of the other scatters, 2,100 MB, was measured when a worker built its tables cold."""

SECTOR_WORKER_MEMORY_MB = 3050
"""int: What one worker takes at most while it generates sectors (PERF.226). Measured on a populated database (a galaxy of one
shell, 4 workers), the worker's and its horse's memory together: 2,083 MB for the starting sector and 2,370 MB for a neighbourhood
sector, 10 percent more with the margin. The horse alone showed 1,487 MB (resident, which counts the worker's tables but not the
worker beside it), and a count by that figure let 10 workers in where 5 fit a 16 GB host."""

STAGE_MEMORY_SHARE = 0.5
"""float: The share of the host's memory a stage's workers may take (less than `HOST_MEMORY_SHARE`): MySQL's buffer pool, Redis,
the web server and the kernel keep the rest, so a stage that grows past its figure is slowed, not killed. 50 percent of a 16 GB
host in 700 MB workers is 11, in 2,100 MB ones 3, in 3,050 MB ones 2."""

STAGE_WORKER_MEMORY_MB = {
    "scatter-mass": SCATTER_WORKER_MEMORY_MB, "scatter-luminosity": SCATTER_WORKER_MEMORY_MB,
    "scatter-backfill": SKY_SHARE_WORKER_MEMORY_MB, "scatter-band": SCATTER_WORKER_MEMORY_MB,
    "generate-start-sector": SECTOR_WORKER_MEMORY_MB, "generate-neighborhood": SECTOR_WORKER_MEMORY_MB,
    "generate-sectors": SECTOR_WORKER_MEMORY_MB, "generate-unplaced-sector": SECTOR_WORKER_MEMORY_MB,
    "scatter-phenomena": LEAN_WORKER_MEMORY_MB, "fill-map-luminosity": LEAN_WORKER_MEMORY_MB,
    "rescore-habitability": LEAN_WORKER_MEMORY_MB, "refill-effect-ledger": LEAN_WORKER_MEMORY_MB, "draw-nebula-field": LEAN_WORKER_MEMORY_MB, "empty-nebula-field": LEAN_WORKER_MEMORY_MB, "clear-bright-stars": LEAN_WORKER_MEMORY_MB, "store-scatter-limits": LEAN_WORKER_MEMORY_MB,
    "empty-phenomenon-scatter": LEAN_WORKER_MEMORY_MB, "drop-phenomenon-indexes": LEAN_WORKER_MEMORY_MB,
    "link-neighbors": LEAN_WORKER_MEMORY_MB, "settle-sector-paths": LEAN_WORKER_MEMORY_MB,
    "warm-galaxy-nodes": LEAN_WORKER_MEMORY_MB, "warm-galaxy-deep-nodes": LEAN_WORKER_MEMORY_MB,
    "warm-sectors": LEAN_WORKER_MEMORY_MB, "warm-systems": LEAN_WORKER_MEMORY_MB, "warm-phenomena": LEAN_WORKER_MEMORY_MB,
    "refresh-block-lights": LEAN_WORKER_MEMORY_MB,      # `one_unit` (no unit-length ceiling): the roll-up is the one unit (OPS.154 ended the uncharted index units of OPS.153)
    "warm-cluster-field": BULK_WORKER_MEMORY_MB,        # not measured yet (PERF.225 cuts it into parts): the larger figure
}
"""dict: The command kinds of the stages with units, and what one worker takes at most in MB (`stage_worker_cap`). Every other
command is one unit by nature (`CommandSpec.one_unit`) or inline, so one worker and nothing to limit. `one_unit` only exempts a command
from the unit-length ceiling: a one-unit command that is cut into several units has a figure too (the block lights: the run of
2026-11-08 held 15 workers at 15,128 MB of 15,992 MB with no figure and no cap)."""

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
"""str: The planetGen checkout, where workers start."""


def _physical_memory_bytes():
    """The host's memory in bytes, or `None` where the OS does not say."""
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None


def stage_worker_mb(kind, conn=None, run_id=None):
    """What one worker of a stage of command `kind` takes at most, MB (PERF.226): its figure in `STAGE_WORKER_MEMORY_MB`, or `None`
    for a command that is one unit or inline (one worker, nothing to limit). With `conn` (OPS.157) the figure is the claim the gate makes for a unit of the kind
    (`profileclaim.claim`: from the job profile when this machine has history, held for the run `run_id`), which is the typed figure with no history."""
    typed = STAGE_WORKER_MEMORY_MB.get(kind)
    if conn is None or typed is None:
        return typed
    from planetgen.queue import profileclaim
    return profileclaim.claim(conn, f"unit:{kind}", typed, token=run_id) or typed


def stage_worker_cap(kind, conn=None, run_id=None):
    """
    The most workers a stage of command `kind` may use at once (PERF.220, PERF.226): `STAGE_MEMORY_SHARE` of the host's memory in
    workers of what that stage's worker takes (`stage_worker_mb`), at least 1. `None` for a command that is one unit or inline or
    where the OS does not say how much memory there is. What a stage generates does not depend on how many workers share it
    (`work.task_seed`), so the cap only changes the time. With `conn` the workers are sized by the claim from the job profile (OPS.157, see `stage_worker_mb`).
    """
    per_worker_mb = stage_worker_mb(kind, conn, run_id)
    total = _physical_memory_bytes()
    if per_worker_mb is None or total is None:
        return None
    return max(1, int(total * STAGE_MEMORY_SHARE / (per_worker_mb * 2 ** 20)))
