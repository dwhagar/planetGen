### Changed
- The bright-star backfill after a `galaxy` run (GEN.98) goes out from the run's edge: the backfill distance past the farthest generated sector in every direction, not a radius around the starting sector. `--backfill-from` is now `edge` (the default) or `none`; the `requested` and `all` choices and the Generate page's "backfill from every generated sector" box are gone.
