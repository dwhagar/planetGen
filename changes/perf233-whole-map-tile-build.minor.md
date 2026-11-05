**Headline:** One command, `python3 -m planetgen.cli.tiles`, builds the whole Galaxy Map's tile cache from scratch, and a nightly timer runs it after the orbit update.

### Changed

- **The whole-map build (PERF.233).** `planetgen.cli.tiles` builds every node of every level of the whole galaxy (a block of each level that holds a generated sector or a placed object, and each generated sector) into a new cache beside the live one and swaps it in when it is complete, so visitors never meet an empty cache and a run that is killed leaves the old cache serving. It builds a level at a time, `tile_cache.build_workers` processes at once (0 is one per core), prints a line per level, ignores the cache budget while it runs, checks the data once more before the swap and refuses to start unless twice the projected size is free on the cache's disk.
- **`cache.max_mb` counts only the nodes visits asked for.** The nodes built ahead (the build's and the warm-up's) are no longer counted against the budget or pruned, so a full cache no longer makes every on-demand node get deleted at the next prune.
- **A nightly timer.** `examples/maintenance/planetgen-tiles@.service` and `.timer` run the build at 05:00 after the orbit update, then `planetgen warm` rebuilds the lists, objects and search; `install-maintenance-timer.sh` installs and enables them (`--skip-tiles-timer` leaves them out). Documented in INSTALL.md, `docs/reference/cli.md` and `docs/design/operations/cache-warm-up.md`.

generation-output: unchanged

TODO items resolved: none (PERF.233 is done with this; the TODO thread retires it).
