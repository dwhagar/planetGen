#!/usr/bin/env python3
# planetgen.cli.warm_map

"""
Builds the site's caches ahead of the first visitor (MAP.134, PERF.104): the Galaxy Map's opening
view, every list, every sector, system and phenomenon page, and the search, so the first visit
after a release (or after a cache was cleared) isn't the one who waits for it. update.sh runs it in
the background after every update; run it by hand after clearing a cache. Run it as the web server's
user so the files it writes are readable by Apache.

    python3 -m planetgen.cli.warm_map [--only FAMILY ...] [--cap N] [--no-rq] [--timeout SECONDS]

Each family (`planetgen.web.warmup.FAMILIES`) is a job of its own on the Redis queue, so each shows
on the admin runs page (`/admin/queue`) with its state and its result line; the script waits for them
and prints one line per family, with how many pages it built, how many were already cached and how
many a cap left out (`page_cache.warm_max_objects`, or `--cap`). `--no-rq` runs the families here,
one after the other. It exits non-zero if any family failed or could not start, saying which; a
Redis that does not answer, or an update still running when the jobs would be queued, is an error
printed, never a silent skip. It waits for a running update to finish before it queues the jobs.

The answers go to Redis and to the disk cache the web server reads, so they are the ones a visitor
gets; each is kept only while the data it read is unchanged (docs/design/cache-warm-up.md). When
every family succeeds it also ends the maintenance the sector-grid switch started
(`maintenance.SECTOR_SWITCH`, UX.112): the site is closed only while the maps are being rebuilt.
"""

import argparse
import sys
import time

from planetgen.admin import maintenance, updateguard
from planetgen.queue import api_jobs
from planetgen.web import warmup
from planetgen.web.helpers import db_name

POLL_SECONDS = 1.0
"""float: How often the jobs are asked how they are doing."""

DEFAULT_TIMEOUT_SECONDS = 1800
"""int: How long the script waits for the jobs (update.sh also stops it after this)."""

UPDATE_WAIT_SECONDS = 1200
"""int: How long it waits for a running update to finish before queueing the jobs."""


def parse_args(argv=None):
    parser = argparse.ArgumentParser(prog="python3 -m planetgen.cli.warm_map", description=__doc__.split("\n\n")[0])
    parser.add_argument("--only", action="append", choices=list(warmup.FAMILIES), metavar="FAMILY",
                        help=f"warm only this family (repeatable); one of: {', '.join(warmup.FAMILIES)}")
    parser.add_argument("--cap", type=int, default=None,
                        help="most sector, system and phenomenon pages to warm (default page_cache.warm_max_objects)")
    parser.add_argument("--no-rq", action="store_true", help="run the families here instead of as queued jobs")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_SECONDS,
                        help="seconds to wait for the queued jobs (default %(default)s)")
    return parser.parse_args(argv)


def _wait_for_update():
    """Waits for a running update to finish: the queue takes no job during one (OPS.53)."""
    if not updateguard.is_updating():
        return
    print("The update is still running: waiting for it to finish before queueing the warm-up jobs.", flush=True)
    deadline = time.monotonic() + UPDATE_WAIT_SECONDS
    while updateguard.is_updating():
        if time.monotonic() > deadline:
            raise warmup.WarmupFailed(f"the update was still running after {UPDATE_WAIT_SECONDS} s")
        time.sleep(2.0)


def _run_inline(db, families, cap):
    results, failures = {}, {}
    for family in families:
        try:
            results[family] = warmup.run_family_job(db, family, cap)
        except Exception as exc:  # noqa: BLE001 -- reported below, the next family still runs
            failures[family] = str(exc)
        else:
            print(results[family]["line"], flush=True)
    return results, failures


def _run_queued(db, families, cap, timeout):
    _wait_for_update()
    try:
        jobs = {family: api_jobs.submit(warmup.run_family_job, db, family, cap) for family in families}
    except api_jobs.NoQueue as exc:
        raise warmup.WarmupFailed(f"the warm-up jobs could not be queued: {exc} (--no-rq runs them here)") from exc
    print(f"Queued {len(jobs)} warm-up jobs: {', '.join(jobs)}. Their progress is on the admin runs page.", flush=True)
    results, failures = {}, {}
    deadline = time.monotonic() + timeout
    waiting = dict(jobs)
    while waiting:
        for family, job_id in list(waiting.items()):
            job = api_jobs.status(job_id)
            if job is None:
                failures[family] = "the job disappeared from the queue"
            elif job["state"] == "succeeded":
                results[family] = job["result"]
                print(job["result"]["line"], flush=True)
            elif job["state"] == "failed":
                failures[family] = job.get("error") or "the job failed"
            else:
                continue
            del waiting[family]
        if waiting and time.monotonic() > deadline:
            for family in waiting:
                failures[family] = f"still running after {timeout} s"
            break
        if waiting:
            time.sleep(POLL_SECONDS)
    return results, failures


def main(argv=None):
    args = parse_args(argv)
    families = args.only or list(warmup.FAMILIES)
    app = warmup.build_app()
    with app.app_context():
        db = db_name()
    # Registered for the update guard (OPS.53), which lists it and stops it, the caches being safe to build again.
    with updateguard.running("warm-map", "Site cache warm-up"):
        try:
            if args.no_rq:
                results, failures = _run_inline(db, families, args.cap)
            else:
                results, failures = _run_queued(db, families, args.cap, args.timeout)
        except warmup.WarmupFailed as exc:
            print(f"error: {exc}", file=sys.stderr)
            sys.exit(1)
    for family, reason in failures.items():
        print(f"error: {family} was not warmed: {reason}", file=sys.stderr)
    if failures:
        print(f"Cache warm-up failed for {len(failures)} of {len(families)} families; the first visit to what they "
              "cover builds it itself.", file=sys.stderr)
        sys.exit(1)
    print(f"Cache warm-up done: {len(results)} families.")
    if set(families) != set(warmup.FAMILIES):
        return    # a part of the warm-up does not say the maps are rebuilt
    # The Galaxy Map is back on the new sector grid (MAP.206): the site reopens (UX.112). An admin's own
    # maintenance, and one with another cause, stays.
    try:
        if maintenance.turn_off_for(maintenance.SECTOR_SWITCH, by="map rebuild"):
            print("Maintenance mode is off: the Galaxy Map is rebuilt.")
    except maintenance.MaintenanceWriteError as exc:
        # The site stays closed, so say so where it is read, and fail: the admin page shows it too.
        print(f"error: the caches are built but the site is still in maintenance mode. {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
