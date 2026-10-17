# Selecting a sphere of the galaxy with a database query

Boss (2026-10-15, 4:36 PM EDT): "Often times we need to find all objects within a certain region. The most natural way to do that is with a sphere. I was wondering if we could design a way so that we can query the database as to what objects are within spherical bounds using a 3D matrix or some other model of a sphere to build a DB query that would naturally only return wanted objects rather than searching through nearby sectors." And (4:38 PM EDT): "The key is what I had in mind for the unique ID so every object has a unique key? If that's what you mean by cell key?"

Informs: GEN.223 (every object stored as a point), GEN.230 and GEN.233 (the 96-bit layout and the 262 kpc reach), NAV.43 to NAV.45 (what is within N pc), NAV.12 (nearest star), the scatter tables, and the Galaxy Map's "what is in this region" reads.
Status: research, 2026-10-15; nothing here is built. The measurements ran on real MySQL 8.0.46 with synthetic data shaped like this galaxy. Boss's rule of 5:27 PM EDT (early alpha, the schema is not a constraint) is what lets this note propose columns and indexes freely.

Evidence tags: **[M]** measured here, on the setup in section 3; **[C]** computed here from the code or from measured numbers; **[E]** extrapolated, with the rule used; **[S]** read in the repository. Scripts, raw results and logs: `/mnt/project-files/research/spherical-region-queries/`.

## Summary

**Short answers to Boss's two questions.**

1. **Yes, a sphere can be asked of the database directly, and the "3D matrix" is a nested grid of cubes.** Give every object one integer, the **cell key**: the position cut into 0.25 pc cubes, with the three cube numbers woven into one number bit by bit (a Morton, or Z-order, key). Cubes that sit inside one bigger, aligned cube have keys that are one unbroken run of numbers. To ask for a sphere the program lays the sphere over that nested grid (big cubes in the middle, small ones along the surface; about 1.5 ms in numpy), merges the cubes into at most 256 runs of keys, and the index reads exactly those runs. Everything outside the sphere's cubes is never touched. Then one exact distance test removes the cube corners [M].
2. **The cell key is not the unique ID.** The ID (80 bits today, 96 bits under GEN.230 and GEN.233) answers "which object is this, forever": birth sector plus serial, and it must never change. The cell key answers "where is it right now": it is a number derived from the position, many objects share one, and it changes when an object moves into another cube. If the ID were the key, an object that drifted out of its birth cube would be filed in the wrong place, or its ID would have to change and break every reference to it. What the two share is the idea: the ID is a birth-place key for finding an object, the cell key is a now-place key for finding what is near. In the database the pair `(cell_key, id)` is unique, which is all the index needs. The sector address (shell, band, wedge) stays too, as the generator's lookup; section 7 explains why it does not make a good sphere key.

**Recommendation** (section 6):

- One `BIGINT UNSIGNED` column **`cell_key`** on every table that holds a position: a Morton key of 21 bits per axis at 0.25 pc, covering plus or minus 262,144 pc (exactly the 262 kpc reach of GEN.233), computed by MySQL itself as a `STORED` generated column from the position columns of GEN.223, so it can never disagree with them [M: equal to the Python key on all 1,000,000 test rows].
- One **covering index** `(cell_key, x, y, z)` on those tables. It answers the whole query without opening a table row [M: 2 to 4 times faster than a plain key index at 50 pc].
- Queries carry **at most 256 key runs**, joined as a `VALUES` table with a forced join order (`STRAIGHT_JOIN` and `FORCE INDEX`), then the exact distance test in the same statement. That shape is 1.5 to 4 times faster than a long `OR` list [M], and without the forced order MySQL once chose a full scan that took 4.2 seconds instead of 70 ms [M].
- Measured at 50 pc, about 56,000 systems at the Sun's density: **42 ms at 10 million rows**. The box over indexed x, y, z took 705 ms, the production path (box over sector centres, then join) took 426 ms, and walking the sectors took 201 ms. At 10 pc: 3 ms against 139, 40 and 5. At 1 pc: 1.8 ms against 18, 16 and 0.8 [M].
- **A key does not make a dense region cheap.** Cost follows the number of systems inside the sphere, about 2 million candidate rows per second [M]. A 50 pc sphere at the core's 15 systems per cubic parsec holds 7.9 million systems [C], so dense areas need a count that stops early, a radius ladder for "nearest first", and a cap (section 6.5).
- The key matters when a **large dense region** is generated, and does little for the sparse halo and intergalactic scatter, where a plain box on indexed x, y, z is already as fast [M: section 4.5]. It costs almost nothing to keep: a move of 0.01 mpc, the update system's threshold, changes the cell in 1.2 in 10,000 updates [M], and the nightly update only touches objects that are due.

Eight working titles for TODO items are in section 9.

## 1. What the code does today

[S] `db/near.py` (`objects_within`, NAV.43) asks `sectors` for every sector whose centre lies in the box `centre ± (r + 1.5 edges)` (index `idx_sectors_center`), joins `star_systems` on `sector_id`, turns each system's sector-local position into the galaxy frame in Python (`local_to_galaxy_pc`) and measures the exact distance there. Phenomena and standalone facilities use their own `center_x_pc, center_y_pc, center_z_pc` columns with a box (`idx_*_center` composite indexes, e.g. `idx_black_holes_center`). Planets, moons, belts and comets are placed at their system's position and found through `star_system_id`.

Two properties of that design set the limits:

- **A composite index on (x, y, z) is a slab scan.** MySQL can use only the first column of the index as a range, so the box on `sectors` reads every sector row with the right x, however far away in y and z. In a generated region of radius R the slab holds about the fraction `r / R` of all rows. Measured on the 10 million row Sun dataset (R = 274 pc), a 50 pc box reads 2.1 million rows to return 56,000 (38 times too many), and a 5 pc box reads 215,000 to return 64 [M]. The slab grows with the generated area, so the cost of today's path grows as the galaxy fills.
- **An object with no sector cannot be found.** The path starts from `sectors`. GEN.223 stores every object as a point of its own, including objects with no sector; a sphere query that starts from the sector table cannot return them.

## 2. What was asked of the candidates

Typical radii: 1 to 5 pc (the nearest star, NAV.12; a neighbour list), 10 to 25 pc (what is around this system), up to 50 pc (NAV.43 to NAV.45's limit). Typical places: a dense core cell, a spiral arm, the thin disk at the Sun, a sparse halo cell. The query must be correct (return exactly the objects whose position is within r), must not depend on sectors existing, and must stay cheap to maintain when the nightly update moves objects.

The designs, with the lettering of the brief:

| | Design | One line |
|---|---|---|
| a | Box over indexed x, y, z | A range on each axis, then the exact test. Variants: one composite index (x, y, z), three single indexes. |
| b | Cell key ranges | Morton or Hilbert key of the cube holding the object; the sphere becomes key runs. Variants: run count, `OR` list or `VALUES` join, plain or covering index, key in a side table. |
| c | Uniform grid | One integer per cube of fixed size (2, 4, 8, 16 pc), one z-run per (x, y) column the sphere crosses. |
| d | MySQL spatial index | An R-tree on a `POINT(x, y)` (2D only, SRID 0), z checked afterwards. |
| e | Walking the sectors | Today: box over sector centres plus join (`sector_current`). Also: an `IN` list of every sector address the sphere touches (`sector_in`), and merged runs of the sector address as a key (`sector_ranges`). |

## 3. Method

**Data** [M]. Five synthetic datasets drawn from the repository's own density model (`galaxy/density.py`, default Milky Way shape: thin disk, thick disk, boxy bulge, two arms, the halo floor) at the model's absolute density (0.0985 systems per cubic parsec at relative density 1, `LOCAL_STELLAR_DENSITY_LY3`). Positions are real galaxy-frame points in parsecs; sector addresses are the real `sector_address_at` of the spherical grid.

| Name | Rows | What | Region | Systems per pc³ |
|---|---:|---|---|---:|
| S1M | 1,000,000 | a generated ball around the Sun (8,200 pc out) | radius 126 pc | 0.12 |
| S10M | 10,000,000 | the same, 10 times larger | radius 274 pc | 0.12 |
| C1M | 1,000,000 | a generated ball around the galactic centre | radius 25 pc | 15.2 |
| C10M | 10,000,000 | the same, 10 times larger | radius 54 pc | 15.1 |
| G1M | 1,000,000 | a thin sample of the whole model galaxy (3.3e-6 of its roughly 3e11 systems) | 30 kpc across | 4e-7 near the Sun |

A generated galaxy grows as blobs around where Boss generates, so the S and C sets are the realistic shape; G1M is the shape of the halo and intergalactic scatter. A ball of 10 million systems at the Sun's density is a small generated area (a 274 pc radius is under 2% of the disk's 15,000 pc radius); the core set has the density of the densest cells.

**Table** [M]. Each row has `id` (assigned in sector order, as the fill order does), `x, y, z` (`DOUBLE`, parsecs), a 96-byte pad (a real `star_systems` row is far wider; section 10), the keys of every design, and the indexes of every design at once, so each query picks its index with `FORCE INDEX` against identical data. 238 bytes per row; a plain key index is 24 bytes per row and a covering index (key, x, y, z) is 50 [M, 10 million rows: 238 MB and 503 MB].

**Machine and server** [M]. 4 vCPU, 15 GB, Ubuntu 24.04, MySQL 8.0.46 (the Ubuntu package, the same MySQL 8.0 as `environment/setup-planetgen-env.sh`, not MariaDB), InnoDB buffer pool 7 GB, no binary log, performance schema off, client on a Unix socket on the same machine, one connection. Warm cache. The two 10-million-row tables are about 7 GB each with all test indexes, so their runs are at the edge of the buffer pool and include some page reads.

**Queries** [M]. Centres are drawn from the dataset's own objects ("what is near this system"): 24 centres at 1M, 16 at S10M, 6 at C10M. Six radii (1, 2, 5, 10, 25, 50 pc), the tables below show five. Each design/radius cell runs every centre once after a warm-up, stopping after 8 seconds (4 seconds for C10M) with at least 3 to 4 queries. The reported time is the median of **planning in Python plus the statement** (`SELECT COUNT(*)` with the exact distance test in SQL), so a design that needs a long plan pays for it. **Every one of the thousands of counts was compared with a brute-force count in numpy and they all agree** (column `wrong` in the logs is 0 throughout).

**Key code** [M]. `keys.py`: Morton and Hilbert (Skilling's algorithm, vectorised; consecutive Hilbert keys are face neighbours, tested on a 16 by 16 by 16 lattice) at 0.25 pc, 18 bits per axis in the benchmark (a 65 kpc box; the recommendation uses 21 bits, section 6.1, which changes the constant and not the shape of the cost). The cover: start at the coarsest level whose cubes are at most 4 radii wide, then refine level by level (each cube is outside, wholly inside or crossing the sphere), keep emitting the whole inside cubes, stop when the merged run count would exceed the budget, emit the crossing cubes at the last level as runs, merge adjacent runs. The result always contains the whole sphere (tested against random points).

## 4. Results

Times are milliseconds, median, planning plus statement. The italic row is how many systems the sphere holds, so the work done is visible.

### 4.1 Sun, 10 million rows (S10M)

| Design | 1 pc | 5 pc | 10 pc | 25 pc | 50 pc |
|---|---:|---:|---:|---:|---:|
| *systems inside the sphere (median)* | *1* | *64* | *468* | *7,216* | *56,666* |
| (a) box on composite index (x, y, z) | 18 | 76 | 139 | 355 | 705 |
| (a) box on three separate indexes | 267 | 493 | 767 | 1,541 | 2,487 |
| (d) spatial index on (x, y), exact check on z | 11 | 52 | 143 | 564 | 1,724 |
| (e) today: box over sector centres, then join | 16 | 24 | 40 | 120 | 426 |
| (e) walk the cells: `IN` list of every sector the sphere touches | 0.8 | 2.1 | 5.2 | 35 | 201 |
| (c) uniform 4 pc grid, one z-run per column | 0.7 | 2.2 | 6.6 | 36 | 376 |
| (b) Morton, 256 ranges, `OR` list, plain key index | 1.6 | 2.9 | 5.0 | 35 | 261 |
| (b) Morton, 256 ranges, `OR` list, covering index | 1.7 | 2.7 | 3.8 | 24 | 164 |
| (b) Hilbert, 64 ranges, `OR` list, covering index | 3.3 | 4.1 | 4.7 | 15 | 84 |
| (b) Hilbert, 64 ranges, `VALUES` join, covering index | 3.4 | 4.0 | 4.5 | 11 | 53 |
| (b) Hilbert, 64 ranges, `VALUES` join, plain key index | 3.7 | 5.7 | 10 | 41 | 224 |
| **(b) Morton, 256 ranges, `VALUES` join, covering index (recommended)** | 1.8 | 2.7 | 3.1 | 8.9 | 42 |

### 4.2 Sun, 1 million rows (S1M)

| Design | 1 pc | 5 pc | 10 pc | 25 pc | 50 pc |
|---|---:|---:|---:|---:|---:|
| *systems inside the sphere (median)* | *1* | *62* | *502* | *7,538* | *50,676* |
| (a) box on composite index (x, y, z) | 4.4 | 16 | 32 | 79 | 184 |
| (a) box on three separate indexes | 47 | 53 | 95 | 249 | 309 |
| (d) spatial index on (x, y), exact check on z | 2.3 | 10 | 33 | 146 | 460 |
| (e) today: box over sector centres, then join | 4.0 | 6.8 | 12 | 46 | 161 |
| (e) walk the cells: `IN` list of every sector the sphere touches | 0.7 | 1.6 | 4.0 | 26 | 140 |
| (c) uniform 4 pc grid, one z-run per column | 0.6 | 1.3 | 3.7 | 28 | 279 |
| (b) Morton, 256 ranges, `OR` list, plain key index | 1.4 | 2.7 | 4.2 | 31 | 182 |
| (b) Morton, 256 ranges, `OR` list, covering index | 1.5 | 2.6 | 3.7 | 23 | 136 |
| (b) Hilbert, 64 ranges, `OR` list, covering index | 3.3 | 4.3 | 4.8 | 15 | 70 |
| (b) Hilbert, 64 ranges, `VALUES` join, covering index | 4.2 | 4.8 | 6.2 | 14 | 48 |
| (b) Hilbert, 64 ranges, `VALUES` join, plain key index | 4.4 | 5.6 | 8.0 | 34 | 113 |
| **(b) Morton, 256 ranges, `VALUES` join, covering index (recommended)** | 1.4 | 2.5 | 3.1 | 9.0 | 36 |

### 4.3 Core, 1 million rows (C1M)

| Design | 1 pc | 5 pc | 10 pc | 25 pc | 50 pc |
|---|---:|---:|---:|---:|---:|
| *systems inside the sphere (median)* | *70* | *7,620* | *51,392* | *409,826* | *1,000,000* |
| (a) box on composite index (x, y, z) | 19 | 85 | 198 | 548 | 695 |
| (a) box on three separate indexes | 80 | 241 | 304 | 569 | 714 |
| (d) spatial index on (x, y), exact check on z | 11 | 139 | 464 | 1,697 | 2,377 |
| (e) today: box over sector centres, then join | 34 | 109 | 261 | 827 | 962 |
| (e) walk the cells: `IN` list of every sector the sphere touches | 14 | 55 | 142 | 611 | 1,112 |
| (c) uniform 4 pc grid, one z-run per column | 5.0 | 35 | 127 | 1,050 | 4,842 |
| (b) Morton, 256 ranges, `OR` list, plain key index | 1.9 | 25 | 150 | 1,287 | 2,769 |
| (b) Morton, 256 ranges, `OR` list, covering index | 1.7 | 19 | 100 | 999 | 1,926 |
| (b) Hilbert, 64 ranges, `OR` list, covering index | 3.3 | 13 | 64 | 665 | 704 |
| (b) Hilbert, 64 ranges, `VALUES` join, covering index | 4.9 | 15 | 58 | 301 | 458 |
| (b) Hilbert, 64 ranges, `VALUES` join, plain key index | 9.9 | 54 | 144 | 761 | 1,154 |
| **(b) Morton, 256 ranges, `VALUES` join, covering index (recommended)** | 1.8 | 8.1 | 36 | 252 | 457 |

### 4.4 Core, 10 million rows (C10M)

| Design | 1 pc | 5 pc | 10 pc | 25 pc | 50 pc |
|---|---:|---:|---:|---:|---:|
| *systems inside the sphere (median)* | *65* | *7,770* | *53,184* | *567,209* | *3,140,468* |
| (a) box on composite index (x, y, z) | 78 | 379 | 780 | 1,658 | 4,400 |
| (e) walk the cells: `IN` list of every sector the sphere touches | 54 | 163 | 437 | 2,845 | 11,530 |
| (c) uniform 4 pc grid, one z-run per column | 8.0 | 81 | 324 | 3,266 | 19,418 |
| (b) Morton, 256 ranges, `OR` list, plain key index | 4.6 | 49 | 296 | 3,458 | 17,696 |
| (b) Hilbert, 64 ranges, `OR` list, covering index | 6.0 | 23 | 103 | 952 | 3,666 |
| (b) Hilbert, 64 ranges, `VALUES` join, covering index | 3.7 | 12 | 53 | 465 | 2,364 |
| **(b) Morton, 256 ranges, `VALUES` join, covering index (recommended)** | 1.9 | 9.3 | 43 | 380 | 1,945 |

C10M was trimmed to the competitive designs and 6 centres: the three-index box, the spatial index and today's path took tens of seconds per query at 25 and 50 pc on the 1 million row core set, and the run only needed to show the trend.

### 4.5 Thin sample of the whole galaxy, 1 million rows (G1M)

| Design | 1 pc | 5 pc | 10 pc | 25 pc | 50 pc |
|---|---:|---:|---:|---:|---:|
| *systems inside the sphere (median)* | *1* | *1* | *1* | *1* | *3* |
| (a) box on composite index (x, y, z) | 0.5 | 0.6 | 0.9 | 2.3 | 4.1 |
| (a) box on three separate indexes | 14 | 14 | 14 | 13 | 14 |
| (d) spatial index on (x, y), exact check on z | 0.4 | 0.2 | 0.2 | 0.4 | 0.7 |
| (e) today: box over sector centres, then join | 0.5 | 0.8 | 1.1 | 2.0 | 3.6 |
| (e) walk the cells: `IN` list of every sector the sphere touches | 0.3 | 0.9 | 1.7 | 11 | 64 |
| (c) uniform 4 pc grid, one z-run per column | 0.4 | 0.4 | 0.5 | 1.5 | 3.9 |
| (b) Morton, 256 ranges, `OR` list, plain key index | 1.5 | 2.4 | 2.3 | 3.2 | 3.3 |
| (b) Morton, 256 ranges, `OR` list, covering index | 1.6 | 2.3 | 2.2 | 2.9 | 3.0 |
| (b) Hilbert, 64 ranges, `OR` list, covering index | 3.1 | 3.7 | 3.8 | 4.2 | 4.3 |
| (b) Hilbert, 64 ranges, `VALUES` join, covering index | 3.3 | 3.8 | 3.8 | 4.3 | 4.2 |
| (b) Hilbert, 64 ranges, `VALUES` join, plain key index | 3.4 | 3.9 | 3.8 | 4.2 | 4.4 |
| **(b) Morton, 256 ranges, `VALUES` join, covering index (recommended)** | 1.5 | 2.5 | 2.5 | 3.3 | 3.3 |

Almost every sphere here holds one system (the centre itself), so this measures overhead only.

### 4.6 What the numbers say

- **Sparse and thin: nothing needs a key.** On G1M the composite-index box takes 0.5 to 4 ms because its slab is tiny (178 to 9,133 rows); the key designs cost 1.5 to 4 ms, almost all of it planning. The halo and intergalactic scatter (a few million rows over 524 kpc, [S] `intergalactic-and-halo-contents.md`) are in this regime; one query shape for every table is a convenience there, not a necessity.
- **Dense and large: the key wins by 1.5 to 45 times over the box (least where the sphere is nearly the whole table) and 2 to 13 times over today's path.** At 10 pc on S10M: 3.1 ms against 139 (box), 40 (today), 5.2 (cell walk). At 50 pc: 42 against 705, 426 and 201. The advantage over the box grows with the size of the generated area, because the box reads a slab and the key reads the sphere (section 5).
- **Where results get huge, nothing is cheap.** At the core, C1M's 50 pc sphere is the whole table (1,000,000 rows) and every design needs 0.5 to 5 s; at C10M it is 3.1 million rows and the best design takes 1.9 s. The recommended shape reads about 2 million candidate rows per second [M]; that, and not the index, sets the cost.
- **The tiny-radius cases favour a fine key.** At 1 pc the covering Morton query takes 1.4 to 1.9 ms on every dataset, the box 4 to 78 ms. The nearest-star question (NAV.12) is a few milliseconds whether the system is in a sparse or a dense place.
- **The spatial index (d) is the worst performer except on the thin sample.** It filters x and y only, so a sphere becomes a column through the whole disk: 4 to 156 times more rows than the sphere holds [M, section 4.7]. On a 2D projection with z handled separately the loss grows with how thick the generated slab is.
- **Hilbert versus Morton.** Hilbert needs about 40% fewer runs for the same cover (44 against 55 at 64-run budgets, r = 50). With the `OR` list and 64 runs each, that makes Hilbert 0 to 50% faster at 50 pc (70 against 74, 84 against 124, 704 against 770 and 3,666 against 4,890 ms at S1M, S10M, C1M, C10M). With the `VALUES` join, where the run count stops mattering, and 256 runs each, the two tie and Morton is slightly ahead (Morton 36, 42, 460 and 1,900 ms; Hilbert 41, 46, 469 and 2,172 ms at the same four, 50 pc). Morton is also 2 to 3 times cheaper to plan (1.4 ms against 3.3 in numpy), 3.5 times cheaper to compute for a column (3.6 s against 12.5 s for 10 million rows), and MySQL can compute it by itself (section 6.2). Morton it is. The Hilbert order stays for the generator's fill order (GEN.101), which is a different job.

### 4.7 False hits

The share of rows the index runs hand to the exact test, as candidates per true hit [M, first five centres]:

| Design | 1 pc | 5 pc | 10 pc | 25 pc | 50 pc |
|---|---:|---:|---:|---:|---:|
| Box (rows inside the cube, not what the index reads) | 1.2 | 1.8 | 1.8 | 1.9 | 1.9 |
| Spatial index on (x, y) | 156 | 60 | 31 | 14 | 7.1 |
| Sector address runs | 117 | 8.2 | 3.3 | 1.7 | 1.3 |
| Grid 4 pc | 12.6 | 3.5 | 2.0 | 1.4 | 1.2 |
| Grid 16 pc | 607 | 17.7 | 8.4 | 2.9 | 1.8 |
| Morton, 64 runs | 1.2 | 2.0 | 2.0 | 2.0 | 2.8 |
| Morton or Hilbert, 256 runs | 1.0 | 1.4 | 1.5 | 1.4 | 1.4 |
| Hilbert, 1,024 runs | 1.0 | 1.2 | 1.2 | 1.2 | 1.2 |

(S10M. Box reads, however, the whole x slab: 42,865 rows at 1 pc, 429,555 at 10 pc, 2,123,811 at 50 pc, against 1, 468 and 56,421 systems inside [M].) The cover at 256 runs is 1.4 to 1.7 times the sphere's volume; the finest cube it uses is 0.25 pc at r = 1, 1 pc at r = 5, 2 pc at 10, 4 pc at 25 and 8 pc at 50, a quarter to a sixth of the radius [M].

## 5. Why the designs behave this way

- **Box on x, y, z (a).** The optimizer can range only on the first index column, so the cost is the x slab: `N × (2r / L)`, L the width of the data in x. Three single indexes (`index_merge`) do worse: it intersects three long primary-key lists (267 ms at 1 pc on S10M for an answer of 1). Spatial index (d): same, in two dimensions.
- **Cell key ranges (b).** The cost is `planning + runs × one index probe + candidate rows ÷ about 2 million per second`, independent of how much of the table is elsewhere [M, section 4]. For the ten times larger table the key queries cost 1.0 to 1.2 times as much; the box costs 4 times as much [M].
- **Why `VALUES` beats `OR`.** With 100 to 300 `OR` terms MySQL's range optimizer plans each range (two index dives per range) before reading a row; at C1M, 10 pc, 256 runs the `OR` form takes 100 ms and the join 36 ms for the same rows [M]. Joined as a derived `VALUES` table, the runs become probes of the nested loop with no planning step. **The join order must be forced:** without `STRAIGHT_JOIN` the 10 million row query with `ORDER BY ... LIMIT 200` took 4.2 s instead of 70 ms, because the optimizer started from a full index scan [M].
- **Why covering matters.** With a plain key index each candidate costs a lookup into the clustered table, a random page. Same ranges, same join, plain key index: 224 ms at S10M, 50 pc; covering: 53 [M]. The data here is in the buffer pool and in near-spatial id order; with colder pages and ids that are not in spatial order the plain index loses more [E, not measured].
- **Grid columns (c).** One run per (x, y) column the sphere crosses: 30 to 540 runs at 4 pc, so the `OR`-list overhead dominates at 50 pc (4.8 s at C1M) and a coarse grid (16 pc) needs 600 candidates for 1 true hit at 1 pc. The multi-level cover reaches 1.4 times the sphere's volume with 256 runs; a single-level grid at the same run count reaches 2 to 4.
- **Sector address as the key (e).** A sphere crosses about 660 (shell, band) arcs at 50 pc (one run each), and the cells are 4 pc, so the cover is coarse near the surface: 117 candidates per hit at 1 pc, 3.3 at 10 pc [M]. It also needs the sector enumeration in Python (46 ms for 11,000 cells at 50 pc [M]).

## 6. The recommended design

### 6.1 The key

- **Cell**: a cube of 0.25 pc. **Key**: the three cube numbers (`floor((p + 262,144) / 0.25)` for p in x, y, z; each 21 bits) woven together, x in the highest bit of each triple: `cell_key = spread(cx) << 2 | spread(cy) << 1 | spread(cz)`. 63 bits, unsigned. The box is plus or minus 262,144 pc, which is the reach of GEN.233 (65,535 shells of 4 pc) with room to spare, so every cell the grid names has a key. When the plan stores a different reach the cell size follows: `2 × reach / 2^21`.
- **Why 0.25 pc**: the densest cell of the core holds 15 systems per pc³, so a 0.25 pc cube holds 0.23 (a cube holds about one system at the core, never many); a finer key costs run count at small radii and gains nothing, a coarser one makes the nearest-star query read its neighbourhood (a 1 pc cube holds 15 at the core).
- **`spread`** is five mask-and-shift steps on a 21-bit integer: shift 32 mask `0x1F00000000FFFF`, shift 16 `0x1F0000FF0000FF`, shift 8 `0x100F00F00F00F00F`, shift 4 `0x10C30C30C30C30C3`, shift 2 `0x1249249249249249`.

### 6.2 The database column

```sql
-- on every table with a position (GEN.223's columns, parsecs). x, y, z are DOUBLE.
cell_key BIGINT UNSIGNED
  GENERATED ALWAYS AS (
    IF(ABS(x) < 262144 AND ABS(y) < 262144 AND ABS(z) < 262144,
       (SPREAD(CAST(FLOOR((x + 262144) / 0.25) AS UNSIGNED)) << 2)
     | (SPREAD(CAST(FLOOR((y + 262144) / 0.25) AS UNSIGNED)) << 1)
     |  SPREAD(CAST(FLOOR((z + 262144) / 0.25) AS UNSIGNED)),
       NULL)) STORED,
INDEX idx_cell (cell_key, x, y, z)
```

`SPREAD(v)` is written out as the five steps above (the macro is spelled in `scripts/extra_tests.py`). **Tested** [M]: with a 65 kpc box (offset 32,768) the generated column equals the Python key on 1,000,000 of 1,000,000 rows; inserting 1,000,000 rows with the generated key and the covering index took 15 s; updating the position of 200,000 rows recomputed the key in 2.7 s; 24 of those 200,000 changed cell after a nudge of 1e-5 pc (1.2e-4, as expected: three axes times 1e-5 over 0.25). **Not tested here**: the `IF` guard for coordinates outside the box (strict mode raises on a negative `CAST ... AS UNSIGNED`; Core must test the guard) and the 21-bit constants (same code, other numbers; Core adds the equality test against the Python function).

An object outside the box has a `NULL` key and is not found by key ranges. The query adds `UNION ALL SELECT ... WHERE cell_key IS NULL AND exact test`; the set is empty unless an ejected star left the reach, and the plan refuses a reach the key cannot hold.

### 6.3 The sphere to runs (Python, in the web layer as well)

```text
runs = sphere_runs(center, r, budget=256)
  -> [(lo, hi), ...] ascending, merged, every sphere point's key inside one run
```

The routine of section 3 (about 1.5 ms in numpy for 256 runs; 0.9 ms at r = 1, 1.9 ms at r = 50 [M]). It uses only integers and numpy. It belongs in `planetgen/galaxy` next to `spherical.py`, with the constants read from the galaxy's plan, and it must not import numba (the web app never does). A pure-Python version should stay under 5 ms [E, not written].

**Built (GEN.276):** `planetgen.galaxy.cellkey.sphere_runs(center, r, budget=256, reach_pc=REACH_PC)`: numpy and integers only, about 2 to 15 ms for 20 to 220 runs (r = 0.3 to 50 pc at the Sun) on the build machine, runs ascending and merged, tested against `cell_keys` on random points inside spheres and on their surfaces.

A budget of 256 is the flat part of the curve: 64 gives 1.3 to 2 times the candidates, 1,024 doubles the planning for 10% fewer candidates [M, `b1024` rows in the logs]. Gap-merging (accepting runs separated by a few empty cubes) was not needed.

### 6.4 The statement

```sql
SELECT p.id, p.x, p.y, p.z,
       (p.x - :cx)*(p.x - :cx) + (p.y - :cy)*(p.y - :cy) + (p.z - :cz)*(p.z - :cz) AS d2
FROM (VALUES ROW(:lo1, :hi1), ROW(:lo2, :hi2) /* ... up to 256 */) AS rg(lo, hi)
STRAIGHT_JOIN star_systems p FORCE INDEX (idx_cell)
  ON p.cell_key BETWEEN rg.lo AND rg.hi
WHERE (p.x - :cx)*(p.x - :cx) + (p.y - :cy)*(p.y - :cy) + (p.z - :cz)*(p.z - :cz) <= :r2
ORDER BY d2 LIMIT 200;
```

`:r2` is r squared. The run bounds are integers written into the statement, not bound parameters. For a table other than `star_systems`, the same text with its name. For extended objects (nebulae, remnants) the search radius is `r + the widest placed radius` (as `near.py` already does with `_widest_placed_phenomenon_radius_ly`) and the exact test measures to the edge, not the centre.

A search over every kind is one such statement per table with the same `VALUES` list. The best shape for that was not measured: a `WITH rg(lo, hi) AS (VALUES ...)` reused by a `UNION ALL` of the per-table selects, against one `located_object` side table that every kind writes to. A side table of (key, id, x, y, z) is as fast to read as a covering index on the table itself [M: side table 36 ms against 36 at S1M, 50 pc], but costs a second row to write for every move (section 6.6); Core measures the union first.

### 6.5 Counts, nearest first, and the dense core

[M, S10M and C1M, 256 runs, covering Morton, `STRAIGHT_JOIN`; median of 12 centres]

| | S10M 5 pc | 10 | 25 | 50 | C1M 5 pc | 10 | 25 | 50 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `COUNT(*)` | 1.0 | 1.5 | 7.6 | 44 | 7.1 | 37 | 267 | 454 |
| count with `LIMIT 10001` (early stop) | 1.1 | 1.5 | 8.2 | 8.6 | 7.5 | 8.9 | 11 | 5.7 |
| top 200 by distance | 1.3 | 2.4 | 12 | 71 | 11 | 62 | 465 | 868 |

- **Counts**: the list pages show "300+" already; `SELECT COUNT(*) FROM (SELECT 1 ... LIMIT 10001) t` costs at most 11 ms even where the sphere holds a million systems [M]. Use it for any total that is shown capped.
- **Nearest first** needs every candidate's distance, so its cost follows the sphere. For "the nearest K" use a **radius ladder**: ask for 2 pc, then 4, 8, 16, 32 (stopping when K rows are found); the first sphere that holds K rows is cheap whatever the place. At the core (15 systems per pc³) 200 systems lie within 1.5 pc; at the Sun within 8 pc [C]. Do not run top-200 at 50 pc at the core: 7.9 million systems [C], about 7 s [E, from 868 ms at 1,000,000 rows].
- **A cap**: the API refuses (or counts only) a query whose early-stop count passes a limit, with the text that says to shrink the radius. The limit is Core's call; 100,000 systems is about 70 ms of scan [M].
- **Paging** keeps NAV.43's keyset cursor `(distance, kind, id)`.

### 6.6 Keeping it current

[M, 1,000,000 rows, 200,000 rows updated in 10,000-row transactions: x, y, z moved, rows per second]

| Index layout | Rows per second |
|---|---:|
| no index beyond the primary key | 607,000 |
| key index only (key unchanged) | 609,000 |
| key index, key changes for every row (worst case) | 148,000 |
| **covering (key, x, y, z), key unchanged** | **145,000** |
| covering, key changes for every row | 148,000 |
| composite (x, y, z) | 147,000 |
| three indexes x, y, z | 65,000 |
| spatial index on (x, y) | 31,000 |
| a second table (key, id, x, y, z) kept in step (my update joins the two tables, which is part of the gap) | 23,000 to 26,000 |

- **The nightly update does not move every row.** [S] `orbital-updates.md` section 1: an object updates when it has moved 0.01 mpc (1e-5 pc), and carries its own `next_update_due`. A star at 20 km/s moves 2e-5 pc a year, so it is due about every six months [C]: 100 million systems means about 550,000 updates a night, 4 seconds at 145,000 rows per second [E]. A run in "faster than real time" mode has the same 0.01 mpc threshold, so the rows moved per game year rise but each update still changes the cell with probability 1.2e-4 [M].
- **Worst case**: every one of 100 million rows moves in one run: 11.5 minutes with the covering index against 2.7 minutes without any [E].
- **Building**: the covering index on 10 million rows took 11 s, a plain key index 7.5 s, a spatial index 193 s [M]; 1 billion rows scales to at least 20 minutes and 9 hours [E, linear; real builds spill to disk and run slower]. The numpy key for 10 million rows took 3.6 s (Morton) and 12.5 s (Hilbert) [M].
- **Size**: 50 bytes per row for the covering index, 24 for a plain one [M]; at a billion systems 50 GB [E]. Tables of wide rows (a real `star_systems` row is hundreds of bytes) make this a few percent of the table.

## 7. Questions the brief asked

- **Cell size and runs.** Fine cell 0.25 pc, cover cubes a quarter to a sixth of the radius, 100 to 220 runs for r = 1 to 50 pc, 1.4 to 1.7 times the sphere in cover, 1.0 to 1.5 candidates per hit [M, section 4.7].
- **Dense core versus sparse outskirts.** Cost is the same function of the number of candidate rows. A key query at the Sun is 3 ms at 10 pc, at the core 36 to 43 ms (51,000 to 53,000 systems), and at the thin sample 2.5 ms (one system) [M].
- **Objects with no sector.** Fully supported: the key is a function of the position only, and the query never reads `sectors`.
- **Scatter placeholders** (`bright_stars`, `phenomenon_scatter`) become real points with the GEN.223 columns and get the same generated key and index. Their (shell, band, wedge) columns stay: the generator looks up "the scatter rows of this cell" with them (an equality lookup, which they serve well), the cell key serves spheres. The scatter pass can keep its bands (a shell's bands are fixed slices) and its order; the key is computed from whatever position it writes. A scatter row that is moved by the nightly update keeps its birth address and gets a new key, as a real object does.
- **Index maintenance** at nightly scale: section 6.6.
- **Cost of the recompute**: the key is recomputed by MySQL when x, y or z is written, 2.7 s for 200,000 updates including the covering index [M]; from numpy 0.36 µs per row for Morton, 1.25 µs for Hilbert [M].
- **Interaction with GEN.223 and its migration.** GEN.223's point columns are the inputs; the key and its index are two lines of DDL in the same migration (Core owns it). The migration of existing rows is the insert measured above (15 s per million rows including the covering index [M]); under the early-alpha rule it can also be a fresh galaxy.
- **The 96-bit ID.** Unchanged, unconnected. The cell key is 63 bits and lives next to the ID; GEN.230's layout is not touched. (The note `intergalactic-and-halo-contents.md` section 1 covers the layout.)
- **Why not use the sector address as the key** (Boss's second question, in full): the address is also a spatial key, and the benchmark has it as `sector_in` and `sector_ranges`. Its cells are 4 pc wedges whose neighbours are far apart in the number order (every shell is a block), so a sphere is about 660 runs at 50 pc and the near-surface cover is 4 pc coarse: 1.3 to 117 candidates per hit against 1.0 to 1.5 for the cell key, and the Python enumeration of the touched cells alone costs 46 ms at 50 pc [M]. It also changes the moment an object leaves its sector and ties queries to the grid layout (any change of the shell rules would rewrite every key), where the cell key depends on nothing but the position.

## 8. What was not measured (and how to close it)

- **Cold cache and wide rows.** Every number is warm, rows are 238 bytes, ids are in sector order. Real `star_systems` rows are wider (the index scan is unaffected, the plain-index and box designs lose more), and ids from several generation runs are less spatially ordered. Core re-times the recommended statement, the plain index and the box on a real generated galaxy once GEN.223 lands; a cold-cache run (restart with a small buffer pool and `O_DIRECT`) is part of that.
- **More than 10 million rows.** Larger tables are extrapolated: B-tree depth rises by one level at about 400 million rows, about +10% on 220 probes [E]; the box's slab grows linearly with N for a given region radius (a billion systems at the Sun's density fill a ball of radius 1,270 pc; a 50 pc box would read about 30 million rows, 10 s at the measured 0.33 µs per row [E]).
- **Several tables in one search** (section 6.4), **concurrency** (one connection), **clumpy real stellar distributions** (clusters make the cover look worse than a smooth density), and **the MySQL version on Boss's server** (8.0.46 here).
- **The 21-bit constants and the `IF` guard** (section 6.2).
- **The Hilbert key's cost in a pure-Python plan** (3 times Morton in numpy).

## 9. Proposed TODO items (working titles; the TODO thread assigns the numbers)

1. **Cell key and covering index on every positioned table** (Core). In GEN.223's migration: `cell_key` generated column, `idx_cell (cell_key, x, y, z)` on `star_systems`, the standalone phenomenon tables, standalone facilities, interstellar comets, asteroid fields and the scatter tables. Equality test: SQL key equals the Python key on a million random points including the box edges and just outside it; the `NULL` guard; tests that no position can be written without its key.
2. **`planetgen/galaxy/cellkey.py`: key and sphere runs** (Core). Pure numpy, no numba; `cell_key(points)`, `sphere_runs(center, r, budget=256)`, constants from the plan's reach (cell size `2 × reach / 2^21`); a test that every point inside a random sphere lands in a run, and a fuzz against brute force.
3. **NAV.43 on the key** (Core, then Interface). `objects_within` reads by the `VALUES`/`STRAIGHT_JOIN` statement for every kind with positions, including objects with no sector and the unbuilt scatter rows ("not generated yet" rows come from the same statement on the scatter tables); the sector-table path goes; extended objects use `r + widest radius`; counts stop early; a plan-shape test (`EXPLAIN FORMAT=JSON` names the covering index and the join order) so a MySQL upgrade that picks a full scan fails a test.
4. **Nearest-first by radius ladder, dense-core guard** (Core). NAV.12 and NAV.43's ordering use the doubling ladder; the API refuses an uncapped count above its limit with a clear message; NAV.44's page shows "300+" from the early-stop count.
5. **Scatter rows keep both keys** (Core, with the whole-sphere train). Point columns and cell key on `phenomenon_scatter` and `bright_stars`; the (shell, band, wedge) columns stay as the generator's lookup and the birth address.
6. **Re-measure on a real galaxy** (Quality, after 1 and 2). The benchmark scripts become `scripts/bench_region_queries.py` (seed, centres, the recommended statement against the box and the sector walk) and run once on a generated galaxy of a few million systems, warm and cold; the numbers go into this note.
7. **Test: the key survives motion** (Quality). The nightly update test moves objects across a cube boundary and asserts the key follows; a round trip bit for bit of the point columns and the key.
8. **Docs** (Docs). `galaxy-coordinate-system.md` and `database-schema.md` describe the cell key; this note moves to "built" when 1 to 3 land; the API reference states the cap and the ladder.

## 10. Reproducing

`/mnt/project-files/research/spherical-region-queries/scripts/`: `keys.py` (Morton, Hilbert, the cover, grid runs), `gen_load.py` (datasets and indexes, with build times), `run_bench.py` and `main_bench.py` (every design, brute-force check), `run_update.py` (section 6.6), `extra_tests.py` (section 6.5 and the generated column), `mk_tables.py` (the tables above). Raw results and logs are in `../results/`.
