**Headline:** A script now measures what a full Galaxy Map tile build costs on a given galaxy: nodes per level, bytes and seconds a node, disk and hours.

### Changed

- **`scripts/measure_tile_cache.py` (PERF.233).** Counts the populated nodes of each level of the block ladder from the stored block rows and the sectors, builds a random sample of each level through `tilecache.fetch_nodes` into a private scratch folder (never the site's cache), and prints per level the count, the sample's mean and 95th percentile bytes and seconds a node and the projected disk and serial time, then the rate of one process against `--workers` processes in chunks of `warmup.NODE_CHUNK`, the hours a full build of every node takes and the site's `cache.max_mb`. It only reads the database. The measured figures for a test galaxy and what they say about `cache.max_mb` are in the project report `perf233-tile-cache-cost-2026-11-04.md`.

generation-output: unchanged
