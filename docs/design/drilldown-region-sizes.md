# Click-centred drill-down: step sizes and region shapes (MAP.146 study)

Research only, 2026-10-09. No code, TODO.md or PR changes. Boss decides.

Question (Boss): MAP.146 keeps the old 243 / 27 / 3 / 1 sector steps, centred on the click. Is there anything more effective or efficient, with more granular control?

## 1. Answer in brief

**Recommendation: split MAP.146 into a frame and a data window, and size the window by "cells across".**

1. **Frame (what the user sees and what a fill or backfill acts on):** a round region (the ADM.30 cylinder) centred exactly on the clicked sector. Described as a layer range, a ring range and, per ring, one slot arc (which may wrap across slot 0 or be the whole ring). Free, exact, follows the cursor for a hover preview.
2. **Data window (what the map reads and caches):** the existing nested cells (243, 27, 3, 1 sectors a side) that cover the frame. The existing per-block stage cache serves them unchanged. A window of 9 cells across touches at most 8 cached stage files.
3. **Size = number of cells across the window**, odd, 3 to 11, changed by wheel or keys; the default of 9 reproduces 243 / 27 / 3 / 1. That gives 19 distinct widths (3 to 2,673 sectors) instead of 4, with no change to the block ladder, the nested-wedge rule or the stage queries.
4. **Fill size is chosen by a system budget, not by width.** The same width holds 1,000 times more systems at the core than at the rim.
5. **Do not force nesting.** With the data on the fixed cells it costs nothing, so the open question can be answered "no".

Ranking (best first):

| # | Option | Verdict |
|---|---|---|
| 1 | Exact frame + aligned data window, size = cells across, existing 9-ary cells | Recommended. Exact centring, 729 aggregate rows per click from at most 8 warm stage files, 19 widths. |
| 2 | Same, on a 3-ary cell pyramid (1, 3, 9, 27, 81, 243, 729) | Upgrade path. 28 widths, steps never above 1.67x, but cells less square (arc 0.84 to 1.25 of an edge against 0.94 to 1.06) and three more levels to build. |
| 3 | MAP.146 default: 243 / 27 / 3 / 1, exact centre, its own child grid | Simple and correct, but every click is a unique cache key and a cold scan: about 3.2 M row reads at width 243. |
| 4 | Free zoom (any width, region follows camera distance) | Best feel, worst cache, unless it sits on option 1's data layer, where it is just option 1 with a smooth wheel. |
| 5 | Other fixed ladders (ratio 2, 3, 4, 5, 7) | Only change the step size; the cache and centring problems of option 3 stay. Ratio 3 needs 7 clicks, ratio 2 needs 10. |
| 6 | Size by a budget of systems, as the navigation rule | Good for fill and warnings, bad for navigation: the width that holds the same number of systems differs about 10x between core and rim, so a click zooms by a different factor in different places. |

What I could not measure is in section 9.

## 2. What the code does today (checked in the repo)

- `galaxy/drill.py` / `galaxyprisms.js`: nested blocks 243, 27, 3, 1 sectors a side (9, 9, 3 per step). A block is `(m, ring, wedge, slab)`. Wedge counts nest (each child ring's count is a whole multiple of its parent's), so every parent wedge line is a child wedge line.
- **There are two caches, and only one is keyed by block.** The MAP.146 text says "cached cube tiles ... keyed by fixed block". That is the **stage cache**: `tilecache.fetch_stage` stores one JSON per block key `m.ring.wedge.slab`, filled by `queryDb.galaxy_stage`. The **cube tiles** (`/api/galaxy/tiles`, the stars, clouds and bright stars) are an octree of fixed cubes chosen from the view radius (at most 3 tiles per axis, 27 per view, star budget 150 to 850 per tile, 4,000 in a finest tile). They are already free-floating and need no change for centred regions.
- `galaxy_stage` for a block reads the generated sectors in the block (member rings and layers by index, one `(ring AND slot BETWEEN)` clause per ring: 243 clauses at m = 243) and groups them into children in Python. Its cost grows with the generated sectors in the block, up to 243^3 = 14.3 M. The design doc's "scale note" already says a per-level count table is needed once generated sectors reach millions.
- Key space today (upper bounds, default galaxy): about 2,600 blocks at 243, 1.29 M at 27 and 892 M at 3. Under exact-sector centring the key is (centre sector, width): 23.9 B centres per width.
- Drawing is not the constraint: each stage draws one prism per child (729 at most) plus stars from the tile budget. None of that depends on how the region is cut.
- Cell squareness of the nested wedges (centreline arc over cell edge, innermost two block rings skipped): 9-ary 0.945 to 1.062; 7-ary 0.93 to 1.08; 5-ary 0.90 to 1.125; 3-ary 0.835 to 1.25.

## 3. What a region really holds

Model: the default galaxy (`planetgen plan` shape: disc scale length 2,600 pc, height 300 pc, bulge 1,580 pc, bar, thick disc), the real skeleton from `build_layer_extents`, 4 pc sectors. Mean systems per cell is E x relative density, with E = 6.31 systems at density 1, averaged round each ring (arm contrast is not in the figures). Totals: 3,764 rings, layers -1,020 to 1,020, 23.92 B cells (matches the 23.9 B in the sector-size study), 301 B systems.

Systems per cell: core 970, ring 100 854, ring 500 209, Sun radius (ring 2,050) 10, ring 3,000 2.3, ring 3,600 0.9. Median cell 2.7 systems; 95th percentile 47; 99th 193.

**Region contents by click position and width W** (round region, radius W/2 sectors, half-height (W-1)/2 layers; sectors / systems):

| click | W=3 | W=9 | W=27 | W=81 | W=243 | W=729 |
|---|---|---|---|---|---|---|
| core (ring 5) | 27 / 26.2 K | 621 / 602 K | 15.7 K / 15.1 M | 418 K / 392 M | 11.3 M / 8.5 B | 304 M / 75 B |
| bulge (ring 100) | 27 / 23.1 K | 621 / 531 K | 15.9 K / 13.6 M | 421 K / 354 M | 11.3 M / 7.9 B | 304 M / 73 B |
| inner disc (ring 500) | 27 / 5.6 K | 603 / 126 K | 15.6 K / 3.3 M | 420 K / 87 M | 11.3 M / 2.2 B | 304 M / 33 B |
| Sun radius (ring 2,050) | 27 / 276 | 585 / 6.0 K | 15.1 K / 154 K | 407 K / 4.1 M | 11.3 M / 96 M | 266 M / 1.3 B |
| outer disc (ring 3,000) | 27 / 63 | 603 / 1.4 K | 15.6 K / 36 K | 412 K / 943 K | 11.2 M / 22 M | 151 M / 258 M |
| edge (ring 3,600) | 30 / 28 | 630 / 582 | 16.1 K / 14.9 K | 429 K / 387 K | 7.1 M / 6.1 M | 56.8 M / 55.4 M |

Reading it:

- The sector count barely depends on where you click (about 15 K at W = 27 everywhere, less at the outline). The star count swings by about **1,000x** (W = 27: 15.1 M at the core, 14.9 K at the rim). A fixed width therefore means a fixed number of sectors to draw, but a wildly different number of systems to fill, count, colour and index.
- Rows to expect in the database are about 20 per system (generation study: about 19 non-system rows per system), so the W = 243 region at Sun radius is about 2 B rows.
- The core zone is small in volume but large in content: rings 0 to 121 hold 0.39% of the cells and **4.1% of all systems** (about 12 B); rings 0 to 243 hold 1.5% of the cells and 14.2% of the systems. Half the cells (53.6%) are inside ring 2,050; 91% of the systems are.
- Height: layers -121 to 121 (the middle slab of the 243 ladder) hold 42% of the cells and 64% of the systems.

## 4. Fill-on-demand cost (why width is the wrong dial for filling)

Rate used: 25 ms per system, one process (generation-performance-study.md, 20 to 27 ms per system, measured 2026-10-09). Eight workers divide this by roughly 8.

| click | W=3 | W=9 | W=27 | W=81 | W=243 |
|---|---|---|---|---|---|
| core (ring 5) | 11 min | 4.2 h | 4 days | 113 days | 2,457 days |
| inner disc (ring 500) | 2 min | 53 min | 22.6 h | 25 days | 622 days |
| Sun radius (ring 2,050) | 7 s | 2 min | 64 min | 28 h | 28 days |
| outer disc (ring 3,000) | 2 s | 35 s | 15 min | 6.5 h | 6 days |
| edge (ring 3,600) | 1 s | 15 s | 6 min | 2.7 h | 43 h |

**Width that holds a given number of systems** (round region; the equal-slot cube is 15 to 20% narrower):

| budget | core (r5) | ring 100 | ring 500 | Sun (r2,050) | ring 3,000 | edge (r3,600) | pole (r20, layer 400) | pole (r200, layer 900) |
|---|---|---|---|---|---|---|---|---|
| 10 K systems | 2 | 3 | 4 | 12 | 20 | 27 | 10 | 25 |
| 100 K | 6 | 6 | 9 | 27 | 44 | 59 | 22 | 55 |
| 1 M | 12 | 13 | 21 | 58 | 94 | 128 | 49 | 127 |
| 10 M | 27 | 28 | 45 | 125 | 208 | 328 | 105 | 283 |

A 100 K system budget is about an hour of one process: a 6-wide region at the core, 27 at Sun radius, 59 at the rim. The ADM.29 estimate (prefix sums, instant) is exactly the tool for this: the fill dialog should offer the widest region that fits a time or system budget and show the estimate before it runs. This is a fill and warning rule. As the navigation rule it makes each click zoom by a different factor (about 10x between core and rim for the same budget), so the breadcrumb loses its meaning (option 6 in the ranking).

## 5. Options compared

### 5.1 Fixed ladders (what sizes)

Zoom factor per step is the ratio; the volume factor per step is its cube. Clicks are from the whole galaxy (about 7,526 sectors across).

| ladder | clicks to a sector | volume per step | cell squareness | notes |
|---|---|---|---|---|
| 243 / 27 / 3 / 1 (9, 9, 3) | 4 | 729x, 729x, 27x | 0.945 to 1.062 | today. Sun radius: 96 M systems -> 154 K -> 276 |
| 625 / 125 / 25 / 5 / 1 (5) | 5 | 125x | 0.90 to 1.125 | |
| 343 / 49 / 7 / 1 (7) | 4 | 343x | 0.93 to 1.08 | |
| 729 / 243 / 81 / 27 / 9 / 3 / 1 (3) | 7 | 27x | 0.835 to 1.25 | |
| 1023 / 511 / ... / 3 / 1 (about 2, odd) | 10 | 8x | n/a (not nested) | |

Odd widths are needed so the centre cell is exact; powers of two cannot centre. A finer ladder gives finer control, but it only adds steps: it does nothing for the cache or for how centring is computed. Ratio 3 also costs squareness, because nested wedge counts must be whole multiples of the parent's (a ratio of 3 forces rounding by up to a third).

### 5.2 Region shape

The ring/slot grid is built so a slot's arc is 0.94 to 1.06 of an edge at every ring. That decides the three shapes:

| shape | physical width across the region | near the axis | verdict |
|---|---|---|---|
| **Same arc (angle)** (MAP.146 default) | grows with radius: at ring 500, W = 243 spans 0.76 to 1.24 of W; W = 729 spans 0.27 to 1.73 | inner edge reaches the axis when click ring <= half-width: the region becomes a wedge from the axis, not centred | Wedge-shaped; worst of the three away from the rim |
| **Same slot count** at every ring | constant (within 6%) | a straight strip along the click's bearing: it does not extend through the axis, so the centroid is up to 25% of W off the click. The zone is 4.1% of all systems at W = 243 and 25.8% at W = 729 | Good outside the core, wrong inside it |
| **Round region** (cylinder around the click, ADM.30's definition) | constant | centred; ring arcs come from the law of cosines; one arc per ring, or the whole ring | Best. Reuses ADM.30; needs slot ranges that can wrap across slot 0 |

Counts at W = 243 (sectors / systems): same arc 12.3 M / 9.2 B at the core, 15.1 M / 128 M at Sun radius; same slots 6.4 M / 4.8 B at the core (it misses the far side of the axis), 14.3 M / 123 M at Sun radius; round cylinder 11.3 M / 8.5 B and 11.3 M / 96 M; a ball (sphere, for nearby search) is 52% of the cube: 7.5 M / 6.1 B and 7.5 M / 68 M.

Slot ranges must be able to wrap: any region within about W/2 of the axis crosses or covers the whole of slot 0's ring, and the old ladder never did.

### 5.3 Where the cache and the query cost sit (the real difference)

| | today's fixed blocks | MAP.146 default (exact centre, own grid) | recommended: exact frame + aligned data window |
|---|---|---|---|
| cache key | block | (centre sector, width): practically unique per click | existing block stage keys; a 9-wide window touches at most 8 files (2 per axis) |
| warm for others | yes | no, except bookmarks and `?sector=` links | yes, shared by every click in the neighbourhood |
| rows read cold, W = 243 | up to the sectors in the block (14.3 M if fully generated), once | about **3.2 M** per click, every click (428 K aggregate blocks + 2.8 M single sectors, found by decomposing the 729 child boxes against the 9-ary cell pyramid at a random alignment; child boxes 27 wide leave 19.5% of their sectors as single-sector leftovers) | 729 aggregate rows from at most 8 stage files; each file is built once |
| data fetched per click | 729 | 729 | (g+1)^3 = 1,000 cells for g = 9 (x1.37), because a frame not on the lattice touches 10 cells per axis; x1.73 for g = 5, x2.37 for g = 3 |
| centring | to the block | exact sector | frame exact; data cells whole |

On the disk cache: the default cap is 200 MB (`tile_cache.max_mb`). Unique keys would turn it into a stream of cold queries and write-then-prune churn.

Other ways to get aligned data, none better: a 3-D summed-volume table in ring/slot space (the slot dimension is not rectangular, and inserts shift every later sum), or caching bricks at the old ladder's cell sizes (that is the recommendation, with no new cache).

### 5.4 Continuous or budget-sized zoom

A free wheel that sets the width continuously, or a region sized to hold N systems, are the most granular and are the worst cache keys. Both are fine as the **UI** if they sit on the aligned data layer: the wheel moves g (3, 5, 7, 9, 11) and, past an end, hops one cell level. That is the recommendation. The pure budget rule is useful as a cap on fill and as an estimate on hover (the page already evaluates the density, `meanDensity` in `galaxyprisms.js`, so a hover estimate needs no request).

### 5.5 The size menu for the recommendation (cells across g = 3 to 15)

| pyramid | distinct widths up to the galaxy | largest step |
|---|---|---|
| existing cells 1, 3, 27, 243; g 3 to 11 | 19 (3, 5, 7, 9, 11, 15, 21, 27, 33, 81, 135, 189, 243, 297, 729, 1215, 1701, 2187, 2673) | 2.45x (297 -> 729) |
| existing cells; g 3 to 15 | 26 | 1.80x (405 -> 729) |
| 3-ary cells 1, 3, 9, 27, 81, 243, 729; g 3 to 11 | 28 | 1.67x (3 -> 5) |

Prisms drawn are g^3: 27 at g = 3, 729 at g = 9, 1,331 at g = 11, 3,375 at g = 15 (a stage today has at most 729 children). Centring error of a data window on a click: at most half a cell, so 16.7% of the width at g = 3, 10% at g = 5, 5.6% at g = 9 and 4.5% at g = 11. The frame itself is exactly centred.

## 6. Edges, poles and the axis

- **Galaxy edge and top of the stack.** Monte Carlo on 150 random clicks per row (round region, half-height W/2). Share of clicks whose region has less than 98% / 75% / 50% of its nominal sectors inside the outline:

| W | clicks weighted by sectors | clicks weighted by systems |
|---|---|---|
| 9 | 2 / 1 / 0% | 1 / 0 / 0% |
| 27 | 3 / 1 / 0% | 1 / 1 / 0% |
| 81 | 10 / 7 / 1% | 1 / 0 / 0% |
| 243 | 44 / 25 / 1% | 5 / 2 / 1% |

  Sliding the centre inward so the region fits (toward the axis and toward the plane) brings all of these to 1% or less, but the click is then no longer central (by up to W/2 at the rim). Clipping at the outline is exact and costs nothing; I recommend clipping, with the label saying "reaches the galaxy edge", and no slide for W <= 81. The poles are thin: at layer 900 the outline reaches only ring 261, at layer 1,000 ring 44 (compare 3,763 on the plane).
- **Axis.** Any click within W/2 of ring 0 gives a region that covers whole rings; it is the core, 4.1% of all systems at W = 243 and 25.8% at W = 729, so it cannot be treated as a corner case. The round region handles it; the strip and the same-arc wedge do not.
- **Slot wrap.** Ranges need `(start, length)` modulo the ring's slot count, and a length of at least N means the whole ring. SQL needs two `BETWEEN` clauses where it wraps. `_block_slot_range` and the ADM.29 column range do not wrap today.

## 7. Nesting (open question 1)

Probability that a click inside a region opens a child that reaches outside it (uniform click; "2D" is a top-down click, "3D" includes the layer):

| step ratio | 2D | 3D | mean slide if clamped inside, as a share of the child's width |
|---|---|---|---|
| 9 | 21% | 30% | 2.8% (about 25% of the width when it does slide) |
| 7 | 26.5% | 37% | 3.6% |
| 5 | 36% | 49% | 5.0% |
| 3 | 56% | 70% | 8.3% |
| 2 | 75% | 88% | 12.5% |

With data on the fixed cells nothing depends on nesting, so there is no reason to clamp: the user gets what they clicked. Up should mean "same centre, one size larger", and Back the history. The breadcrumb needs a new label (centre sector and width) in place of block names. Clamping to the parent is cheap (a slide of about 3% of the width on average at ratio 9) if Boss prefers strict nesting.

## 8. How the recommendation answers MAP.146's three open questions

1. Nesting: no. See section 7.
2. First click from the whole galaxy: centred like the rest. The opening view stays the fixed pre-built view (MAP.134); the first window is cells of 243 sectors, 9 across, centred on the clicked 243-cell.
3. "Same arc or same slot count": neither. Use the round region around the click (ADM.30's definition). Same arc distorts away from the rim (0.76 to 1.24 at ring 500, W = 243); same slots misses the far side of the axis.

The three other items in MAP.146 change as follows: the cube tile cache is unaffected (octree keys, not blocks); the stage cache is kept, not rekeyed, and no new tile type is needed; the fill and backfill actions (MAP.120, ADM.29) take the frame's exact layer, ring and slot-arc ranges (the stats panel takes whole cells); the settle step takes the same ranges.

## 9. What I could not measure

- **Draw time.** No browser or GPU here. Draw cost is from the code's caps (prisms per stage, tile budgets), not a frame time.
- **Database timings.** No MySQL here. Row counts above come from the geometry and a lattice decomposition, not from queries. The 3.2 M row figure is on a cube lattice with a 9-ary pyramid; the real wedge cells differ at the edges. Stage file sizes are unknown (the 200 MB default holds an unknown number).
- **Fill rate.** 25 ms per system is from the 2026-10-09 study and the open PERF items (PERF.43 to PERF.47) may change it; the fill times scale with whatever the rate becomes. The figures are single-process.
- **Star counts.** Means from the density model, averaged round each ring (the arm contrast, +/-40% of the thin disc, is not included, and a region on an arm crest holds more than the table). Cells the model predicts below one system but the outline admits (the outline test uses the maximum over azimuth) hold 0 or 1 in practice. The totals check: 301 B in the model against about 298 B in the sector-size study.
- **Usability.** Whether 5.6% centring error of a data window, or the hover preview, feels right is untested; the frame is exact, so I expect it not to matter.
- **Wheel behaviour.** The 19-width menu is computed, not tried.
- **Monte Carlo** uses 150 clicks per row; shares of 1% to 3% are noisy.
- **Pyramid decomposition.** Does not model the real wedge geometry; the 3-ary and 2-ary rows are indicative only.

## 10. Method and files

Everything is computed from `planetgen.galaxy.geometry` (ring and slot counts), `density` (the model), `skeleton` (the outline, via `build_layer_extents`) and `drill` (wedge counts), with the pieces that import the name generator stubbed out. Scripts in `scripts/`:

- `model.py` builds the per-(ring, layer) mean-systems table. `regions.py` counts a region (same arc, same slots, cylinder, sphere).
- `t3.py`, `t11.py`: region contents. `t4b.py`: width for a system budget. `t5.py`: shares of cells and systems. `t6.py`: wedge squareness by ladder. `t7.py`: lookups to answer a box from a pyramid. `t8.py`: edge and pole fill. `t9.py`: shape widths and axis zone. `t10.py`: block key-space counts.

Run from `src` with `PYTHONPATH=.`; they need numpy and astropy (and scipy for the stubbed import chain).
