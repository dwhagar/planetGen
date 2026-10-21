**Headline:** The cache warm-up now runs in five named phases and says after each one what it built, what was already stored and what the cache budget left out.

### Changed
- **`planetgen warm` runs in phases (PERF.141).** In order: `opening` (the Galaxy Map's opening views and the first page of each list), `built-levels` (the tile octree's first levels), `deep-tiles` (the zoom run in to the centre), `objects` (the sector, system and phenomenon pages) and `search`. `warm-galaxy-tiles` is split in two for it, and the new `warm-galaxy-deep-tiles` command builds the zoom run. After each phase the run prints one line: `Phase 2 of 5, built-levels: 3 built, 125 already stored, 0 left out in 2.0 s.`, or `unchanged, nothing built` when the data under it has not changed since the last run, and `stopped at the cache budget` when `cache.max_mb` cut the tile phases short. A phase name can be given instead of families (`planetgen warm deep-tiles`).

TODO items resolved / Added: PERF.141 (slice 1: the phases and their reports; the automatic start and the admin page follow) / Added: none
