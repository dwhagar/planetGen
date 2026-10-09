# planetgen/cli/worker.py

"""
One RQ worker for planetGen's queues (PERF.24), at the lowered priority
generation workers always run at. `planetgen.queue.redisqueue.
ensure_workers` starts these in burst mode, so each exits once its
queues are empty.

Usage:
    python -m planetgen.cli.worker [--burst] [--name NAME] [--url REDIS_URL] [QUEUE ...]
"""

import argparse
import os
import re
import sys

from planetgen import _version
from planetgen.queue import redisqueue, work

WARM_MODULES = (
    "planetgen.generation.run_plan",
    "planetgen.generation.run_galaxy",
    "planetgen.generation.run_common",
    "planetgen.generation.bright_stars",
    "planetgen.queue.api_jobs",
    "planetgen.web.job_runner",
)
"""tuple: The modules the queue's tasks run in. Imported before the worker
forks, so a work horse starts with them loaded instead of importing them
again for every job (PERF.42)."""


def warm(queues=None):
    """
    PERF.42: loads in the worker, before RQ forks a horse per job, what every
    job would otherwise rebuild in its own process: the generation modules
    (about 2.4 s each job) and, for the generation queue, the bright-star
    sampling tables (`star_population._bright_table`, one per luminosity
    floor and population) at every floor a scatter or backfill uses. Nothing
    here opens a database or Redis connection, so the horses share it safely.

    Args:
        queues (iterable, optional): The queue names served; the tables are
            built when the generation queue is among them (or none are
            given).
    """
    import importlib
    for name in WARM_MODULES:
        importlib.import_module(name)
    if queues and redisqueue.GENERATION_QUEUE not in queues:
        return
    from planetgen import tuning
    from planetgen.generation import star_population
    floors = {tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL,
              *(floor for _out_to, floor in tuning.BRIGHT_STAR_BACKFILL_TIERS)}
    for floor in sorted(floors):
        for population in star_population.POPULATIONS:
            star_population._bright_table(float(floor), population)


def code_version_on_disk():
    """The release in `_version.py` as it is on disk now (the loaded module
    keeps the one this process started with)."""
    path = os.path.splitext(_version.__file__)[0] + ".py"
    try:
        with open(path, "r", encoding="utf-8") as handle:
            match = re.search(r'^__version__\s*=\s*"([^"]+)"', handle.read(), re.M)
    except OSError:
        return None
    return match.group(1) if match else None


def stale(loaded=None, on_disk=None):
    """Whether the code on disk is another release than this worker loaded."""
    on_disk = on_disk or code_version_on_disk()
    return on_disk is not None and on_disk != (loaded or _version.__version__)


def worker_class():
    """
    RQ's worker class with PERF.42's guard: a worker that pre-imports its
    code would keep serving the old release after an update (a horse used to
    import the new files on every job), so before it waits for the next job
    it starts itself over when the release on disk has changed.
    """
    base = redisqueue.worker_class()
    if not hasattr(os, "fork"):
        return base  # a spawned horse is a fresh interpreter anyway

    class WarmWorker(base):
        def dequeue_job_and_maintain_ttl(self, *args, **kwargs):
            if stale():
                self.register_death()
                os.execv(sys.executable, [sys.executable, "-m", "planetgen.cli.worker", *sys.argv[1:]])
            return super().dequeue_job_and_maintain_ttl(*args, **kwargs)

    return WarmWorker


def main(argv=None):
    parser = argparse.ArgumentParser(prog="planetgen.cli.worker", description=__doc__.strip().splitlines()[0])
    parser.add_argument("queues", nargs="*", default=list(redisqueue.QUEUES),
                        help="the queues to serve, first one first (default: all of planetGen's)")
    parser.add_argument("--burst", action="store_true", help="exit once the queues are empty")
    parser.add_argument("--name", help="the worker's name (default: RQ picks one)")
    parser.add_argument("--url", help="the Redis server (default: PLANETGEN_REDIS_URL, else config.json's redis.url)")
    args = parser.parse_args(argv)
    work.lower_priority()
    warm(args.queues)
    connection = redisqueue.connect(args.url)
    queues = [redisqueue.queue(name, connection) for name in args.queues]
    worker = worker_class()(queues, connection=connection, name=args.name)
    worker.work(burst=args.burst, with_scheduler=False)


if __name__ == "__main__":
    main()
