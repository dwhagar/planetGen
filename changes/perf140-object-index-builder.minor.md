**Headline:** The object index fills itself: the warm-up builds it, keeps it current and counts it against the cache budget.

### Changed
- **Object index builder (PERF.140).** `planetgen warm` has a new first family, `object-index` (the command `warm-object-index`), in the `opening` phase. It writes one `object_index` row for every system, star, planet, moon, belt, comet, facility and phenomenon, and a row changes only when what it holds changed. The first run builds each kind; later runs write only the objects edited since (a planet or moon edit reaches its row through its system), and remove the rows of objects that are gone. After `planetgen orbits` or a new generator epoch the index is built whole, since those move objects between sectors without touching the rows.
- **The index is kept first in the cache budget.** What the index takes on disk comes off `cache.max_mb` before the tile cache gets any of it.
- **`planetgen check-db object-index`.** Compares the index with the tables and says how many rows are missing, extra or different; it warns when nothing has built the index yet. It is not in the default list of checks.
- **The admin queue page** shows each kind's rows, when it last ran and the size of the index under the warm-up phases (`GET /api/admin/work/warm` gains `index`).

TODO items resolved / Added: PERF.140 / Added: none
