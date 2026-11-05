#!/usr/bin/env python3
# planetgen.cli.tiles

"""
Builds the whole Galaxy Map's node cache from scratch (PERF.233): every node of every level of the whole galaxy, into a new
cache beside the live one, swapped in when it is complete. Meant for the maintenance timer, nightly after the orbital update
(`examples/maintenance/planetgen-tiles@.service`); run it as the web server's user, so Apache can read what it writes.

    python3 -m planetgen.cli.tiles [--workers N] [--force] [--no-swap]

What it builds and how is in `planetgen/web/tilebuild.py`: the nodes of `nodes.all_node_keys` (a block of each level that holds
a generated sector or a placed object, and each generated sector), a level at a time, `tile_cache.build_workers` processes at
once, with a progress line per level. The live cache keeps serving until the swap, and a run that is killed leaves it as it
was; run it again and it picks up the nodes already built. It refuses to start unless twice the projected size is free on the
cache's disk (`--force` builds anyway). `planetgen warm` rebuilds the other caches (lists, objects, search) after it.

Exit codes: 0 done, 1 failed, 75 an update is running (nothing was built; the timer's next firing tries again), 130 stopped.
"""

import argparse
import sys

from planetgen._version import VersionAction, version_banner
from planetgen.admin import updateguard
from planetgen.db.store import add_mysql_connection_args, mysql_config_from_args
from planetgen.util import log


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Build every node of the whole Galaxy Map into a new cache and swap it in.")
    add_mysql_connection_args(parser)
    parser.add_argument("--workers", type=int, default=0, metavar="N",
                        help="Processes to build in (default: tile_cache.build_workers, 0 there is one per core).")
    parser.add_argument("--force", action="store_true",
                        help="Build although less than twice the projected size is free on the cache's disk.")
    parser.add_argument("--no-swap", action="store_true",
                        help="Leave the finished cache in the staging folder instead of swapping it in.")
    parser.add_argument("--version", action=VersionAction, banner=version_banner("planetgen.cli.tiles"))
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.apply_debug_arguments(args, parser, errors_to_stderr=True)
    if args.workers < 0:
        parser.error("--workers must be 0 or more.")
    config = mysql_config_from_args(args)

    try:
        updateguard.refuse_while_updating("The whole-map build")
    except updateguard.UpdateRunning as exc:
        log.error(f"The whole-map build did not run: {exc}")
        return updateguard.EX_TEMPFAIL

    from planetgen.web import tilebuild
    try:
        with updateguard.running("tiles", "Whole-map build"):
            tilebuild.build_all(config.database, workers=args.workers or None, force=args.force, swap=not args.no_swap)
    except tilebuild.BuildRefused as exc:
        log.error(f"The whole-map build stopped: {exc}. The live cache was not touched.")
        return 1
    except KeyboardInterrupt:
        log.error("Stopped: the live cache was not touched; run it again to pick up the nodes already built.")
        return 130
    except Exception as exc:  # noqa: BLE001 -- any failure is a failed run; the live cache stays
        log.error(f"The whole-map build failed: {exc}. The live cache was not touched.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
