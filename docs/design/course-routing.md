# Course routing and nearby search

How the NAV page finds a route from one system to another, how routing and
the "what is within N pc" search read the galaxy database quickly, how a
route crosses unfilled space, and how a course is shown. Recorded 2026-10-02
from Boss's decision on hop length (NAV.12) and the hop-length study
(`nav-hop-length/report.md` in the project's shared files); extended
2026-10-09 with the measured routing and spatial-index research. The frames a
course is described in (bearing and mark, warp and fold speeds) are in
[navigation-frames.md](navigation-frames.md). How a hop is bent around gravity
wells and asteroid fields is in [course-avoidance.md](course-avoidance.md).
The sector grid the queries use and the order sectors are generated in are in
[galaxy-coordinate-system.md](galaxy-coordinate-system.md) and
[fill-order-curves-and-core.md](fill-order-curves-and-core.md).

Informs: NAV.9, NAV.10, NAV.11, NAV.12, NAV.34 (done), NAV.38 (done), NAV.36, NAV.39, NAV.41, NAV.42, NAV.43, NAV.44, NAV.45, NAV.47, NAV.48, UX.35 (also NAV.25 and NAV.6, which read NAV.10's corridor query)

Status: section 1 and the NAV.34, NAV.38 and TEST.79 pieces are built (checked against main at 27220d8, 2026-10-09). Everything else is planned and not built. Sections 3 to 9 are research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation.

Evidence tags: [S] seen in a search result or file, [C] computed or measured in the research (single-thread CPU seconds, the MariaDB server's own CPU and index entries read, because the machine was shared and wall-clock times were unusable), [R] recalled and unconfirmed. The research environment could only read search-result text, not papers or vendor manuals, and its search budget ran out after one query, so every outside-world fact is [R] unless marked.

## Decisions already taken

- **NAV.12 (Boss, 2026-10-02 01:53Z):** "routs will always find the nearest star they can even if it crosses sector boundaries even across multiple sectors. This is for a game mechanic I need in place and I also want a jump through unknown space is marked in red and glows to draw attention to it. Also UX change here to display the path horizontally and find a way to split it between multiple lines for mobile or limited displays." No hop limit of any kind.
- **NAV.11, NAV.41, NAV.42 (Boss, 2026-10-02 04:19Z):** the route's travel time "using that route assuming each planet gets stopped at", read as a stop at every system on the route; a readable wide course map; a course and distance on every stop.
- **NAV.43 (Boss, 2026-10-02 04:39Z):** select a place and ask what is within some distance in parsecs; the largest distance allowed is 50 pc.
- **NAV.47 (Boss, 2026-10-03 05:38Z):** system to system is by far the preference; only when a jump crosses unknown space does the route look for scattered stars, black holes, neutron stars and quasars as stops, "marked as their kind". Hypervelocity stars are included.
- **NAV.48 (Boss, 2026-10-07 11:47Z):** where unfilled sectors block a course and cannot be bypassed, offer to generate "all sectors between the two points".
- **NAV.45, issue #677 (Boss, 2026-10-08 22:24Z):** the "within N pc" action is offered for any pickable object, down to planets and moons.

## 1. Today

- `queryDb.nav_between` routes two systems in the same sector over that
  sector's own systems (the Sector Local Frame) and anything else over every
  placed system in every generated sector (`_galaxy_frame_positions`, the
  Galactic Standard Frame).
- `navGraph.build_knn_adjacency(positions, k=6)` builds a pure-Python k-d tree
  and links each system to its 6 nearest, made two-way, on every request;
  `navGraph.shortest_path` runs Dijkstra on straight-line lengths.
- Only star systems are stops. Phenomena can be the ends but never a stop;
  bright stars from the galaxy scatter are not; ungenerated sectors add
  nothing. Nothing limits or reports hop length, and warp and fold times are
  for the direct distance only.
- The NAV page lists the route as a vertical `<ol class="nav-route">` and draws
  it on one square map at one scale (`planetgen/web/maps/navmap.py`).

What goes wrong:

- **Islands.** A generated area with 7 or more systems had no link out: 2,000
  generated sectors over the disk split into 714 islands and a cross-galaxy
  course found no route. Fixed by NAV.34 (PR #427): `navGraph.join_islands`
  links each island to its `NAV_ISLAND_LINKS` (6) nearest islands, in rounds
  until one is left.
- **Hidden long hops.** A lone system always links out, however far: one 2 kpc
  above the disk was reached by a route whose last hop was 6,504 ly.
- **Speed.** The rebuild dominates: 0.5 s at 10,000 systems, 16 s at 200,000
  (Dijkstra 0.7 s). A galaxy-scope load of 1e6 systems with the project's
  dict-row pattern takes about 10 s and 830 MB before any graph exists [C].
- **A performance cliff.** The pure-Python tree takes 1.18 s against 0.03 s for
  `scipy.spatial.cKDTree` at 10,000 points [C]. `join_islands` queries every
  point of the smaller island against the larger island's tree for each close
  pair: 15 to 18 s for 10,000 points in 4 islands, 6.9 s for 53,000 points in
  394 islands, against 0.07 s and 0.23 s for a vectorised rewrite returning the
  identical edge set (compared edge by edge, 3,003 links on 106,529 points)
  [C]. The study's 0.2 s join of 714 islands was a different implementation.
  Port it before NAV.10 relies on it.

## 2. No maximum hop length (NAV.12, phase 1)

This drops the study's optional "ship range" field and its "deep-space hop
over 2 sector edges" flag.

- **A route always exists** between two placed endpoints, by joining the
  graph's islands (NAV.34). In the study, linking each island to its 6
  nearest joined all 714 in 0.2 s and gave a cross-galaxy route 1.29 times
  the direct distance.
- **Same-sector routes may leave the sector.** This is a loading rule: load
  the corridor by cells (section 3), never "that sector's systems". The
  third xfail test in `test_route_edge_cases.py`
  (`..._same_sector_route_uses_nearer_stars_next_door`) will pass because
  of graph shape (the near line stars are in the ends' 6 nearest), not the
  cost.
- **The longest hop is shown** (`longest_hop_ly`).
- **Unknown-space jumps.** A hop is flagged when its straight line crosses
  an unfilled sector other than the cells holding its two endpoints. The
  cells come from `galaxyGeometry.sectors_along_segment` (NAV.38, PR #357),
  exact at faces and edges, with browser twin
  `galaxyprisms.sectorsAlongSegment`. `/api/nav` returns per hop
  `unknown_space`, `from`, `to`, `distance_ly`, plus a `sector_id` and
  sector-local position per stop (NAV.42 needs them); these are the keys the
  xfail tests read.
- **The filled set.** `store.filled_sector_addresses(conn)` returns a Python
  set of tuples: 1.56 million addresses took 5.2 s and about 170 MB [C].
  Pack each address into one int64 (ring, layer and slot fit in 49 bits, as
  `provisional_sector_designation` packs them) in a sorted array: 13 MB, 0.5 s,
  and a membership test of 14,000 cells in 0.7 ms [C]. Cache it keyed by the
  existing `galaxy_content_state` token (`_state_token`).
- **Void is not unknown.** Cells outside the galaxy outline
  (`GalaxyBounds.contains`) are flagged "through unknown space" but never
  offered for charting (NAV.48).
- **Cheap tests first.** When both endpoints share a cell or sit in
  face-adjacent filled cells, skip the ray cast; call the exact function
  only for hops longer than about one edge.
- **Cost stays plain length** (section 3.3). Route edge cases are tests
  first (TEST.79, PR #427).
- **Built (NAV.12).** `route.hops` carries `{from, to, distance_ly,
  unknown_space}` per hop and `route.longest_hop_ly` the longest
  (`query.nav_between`, `corridor.unknown_space_flags`: a hop is flagged when
  a cell on its line, other than the two holding its ends, has no `sectors`
  row). A same-sector pair with a galaxy placement is routed in the galaxy
  corridor and shifted back into the sector's frame. The NAV page states the
  longest hop and the number of unknown-space jumps; the red glow is NAV.36
  and the horizontal strip UX.35. Not yet built: the per-stop `sector_id` and
  local position, the packed filled-set cache and the cheap adjacent-cell
  shortcut (the check reads the cells a line crosses with one query per ring
  and layer).

- **Built (UX.35).** The route is an ordered list (`role="list"`) laid out left to right with
  flex wrapping, each stop a link followed by the hop to the next (distance, and "unknown space"
  in italics for a hop through ungenerated sectors); a stop and its hop stay on one line and the
  list wraps as the panel narrows, with no device check. A route of more than nine stops shows a
  short list instead (the first stop, the last three, both ends of the longest hop and of every
  unknown-space hop, with an ellipsis where stops are left out) and keeps the whole route in an
  "All N stops" `<details>`. Checked at 320 px and 900 px.

## 3. Routing that scales (NAV.10, built)

- **Corridor search.** `query.nav_between` at galaxy scope no longer loads
  every placed system. `corridor.positions_near_segment` reads the sectors
  near the straight line between the ends from `sectors`' center index, cut
  into pieces so the read grows with the line's length, then the systems of
  those sectors by their indexed sector id, and keeps those within a
  half-width of the line. The half-width starts at a quarter of the direct
  distance, between 25 and 100 ly, and doubles (to at most 6,400 ly) while
  the route found has a hop longer than half of it: a hop that long may be
  hugging the corridor's edge with stepping stones just outside it. A
  corridor of more than 300,000 systems is not widened. The joined graph
  (NAV.34) still guarantees a route among whatever the corridor holds.
  Because the graph is built from the corridor's systems, a route can differ
  from the one a whole-galaxy graph would give where a nearest neighbour lies
  outside the corridor; the benchmark below shows the length stays within a
  few percent.
- **A\*.** `nav_graph.shortest_path` takes the node positions and uses the
  straight-line distance to the goal as its heuristic. The length is the
  same as Dijkstra's (a test checks it); fewer nodes are visited.
- **Objects near a line.** `corridor.objects_near_segment` returns every
  generated system, star and phenomenon within a distance of a segment, in
  order along it (NAV.6's steering will use it).
- **No new indexes.** The existing center indexes (`idx_sectors_center` and
  each phenomenon table's own) and `idx_star_systems_sector_id` serve the
  boxes, so NAV.10 needs no schema change.
- **Measured.** `scripts/bench_nav.py` builds a synthetic galaxy of any size
  with bulk SQL and times a corner-to-corner route. On 100,000 sectors and
  400,000 systems (a flat grid, 4 pc apart, four systems each), a 4,519 ly
  route took 3.3 s through the corridor and found a 5,162 ly path of 561
  stops; loading every system and rebuilding the graph, as before, took 88 s
  (2.9 s to load, 84.6 s to build). On 10,000 sectors and 40,000 systems: 0.9
  s against 4.9 s, a 2,098 ly path against 2,101 ly with the whole graph.
  Measured on the build server with a test suite running beside it.
- **Research recommendations not taken.** The research proposed building the
  corridor from the sector cells along the segment through the address index
  (1,000 to 52,000 rows in 15 to 60 ms and a route in 2 to 34 ms on a
  1e7-system, 1.56 million-sector test database [C], section 6), a cKDTree
  for the k-nearest graph and the island join, and noted that
  `idx_sectors_center` prunes only on X because it is `(center_x_pc,
  center_y_pc, center_z_pc)`. The built version reads the centre index in
  pieces along the line; revisit the cell-based read if the benchmark above
  grows slow on a full galaxy. The stored `nearest_systems` table cannot back
  routing (3 neighbours per object, within 4 pc, no island data).

### 3.1 Corridor rule and widening

kNN-6 plus island links rebuilt on the corridor's points only, against the
whole-graph optimum, on 52,983 systems in 394 areas, 40 random routes
between areas (median 2,545 pc, 50 hops, 1.21 times the straight line) [C]:

| Half-width | Corridor points (median) | Length over whole-graph optimum (median, p90) | CPU, pure-Python graph |
|---|---|---|---|
| max(20 pc, 2% of L) | 545 | 0.844, 0.909 | 0.05 s |
| max(50 pc, 5% of L) | 1,533 | 0.864, 0.918 | 0.15 s |
| max(100 pc, 10% of L) | 2,895 | 0.896, 0.961 | 0.60 s |
| max(200 pc, 25% of L) | 11,913 | 0.973, 1.000 | 3.15 s |

- **A route always exists whatever the corridor**, because the island join
  runs on the corridor's own graph. Widening is for quality only; this
  replaces NAV.10's "widened if no route is found".
- Corridor routes are shorter than the whole-graph "optimum", which links far
  islands from the global picture. So the rule fixes the result: default
  half-width `max(50 pc, 5% of L)`, one doubling if the route exceeds 1.3
  times the straight line.
- **Density makes the corridor explode:** a 1,000 pc route at 50 pc half-width
  is 7.9e6 pc^3, about 785,000 systems at solar density and 30 million near
  the core. Choose the width by a point budget (for example 50,000): sum
  `sector_stats.actual_systems` along the dilated cells and shrink the width
  until it fits.
- Dense regions give many tiny hops: 58 stops for a 100 pc route in the
  solar patch (2.1 pc per hop) [C]. UX.35 must cope (section 10).

### 3.2 Search algorithm

Same archipelago, 211,138 edges, same 40 routes, median CPU and nodes
expanded [C]: Dijkstra in Python 0.207 s, 28,612 nodes; A* with the
straight-line heuristic 0.041 s, 4,994 nodes (p90 27,009 against 50,525);
bidirectional Dijkstra 0.133 s, 24,128 nodes; weighted A* (heuristic x 1.5)
518 nodes, route 1.023 times optimal; `scipy.sparse.csgraph.dijkstra` over
the whole graph from one source 0.016 s.

A* is exact: every edge costs its Euclidean length and the heuristic is the
Euclidean distance (admissible and consistent, island edges included). Use
scipy's Dijkstra on numpy arrays for plain routes and Python A* when edges
need callbacks (unknown-space penalty, NAV.47). Weighted A* is a fallback for
very long routes at a 2% median loss.

### 3.3 Connectivity and cost

Connectivity graphs on the same 52,983 systems; stretch is route length over
straight line for 120 random pairs between areas [C]:

| Graph | Edges per node | Build CPU | Stretch (median / p90) | Longest hop vs minimax |
|---|---|---|---|---|
| kNN-6 + join to 6 nearest islands (today) | 4.0 | 0.1 s + join | 1.21 / 1.35 | 1.37x |
| Delaunay | 7.4 | 2.1 s | 1.06 / 1.11 | 1.64x |
| Gabriel | 3.0 | 2.9 s | 1.19 / 1.27 | 1.35x |
| Euclidean MST | 1.0 | 2.2 s | 3.47 / 5.25 | 1.00x |

Relative-neighbourhood and kNN-plus-MST graphs reach the minimax hop but
zigzag (stretch 1.51 and 2.38). Keep kNN-6 plus nearest-island links: 1.2
stretch, the cheapest build, a local update. Add MST edges only if a
"smallest possible longest hop" feature (the dropped ship range) returns. For
NAV.12: the longest hop on a cross-gap route is 1.37 times the unavoidable
gap, because Euclidean cost picks long hops.

Cost, 60 pairs 80 to 150 pc apart on a 6-nearest graph [C]: plain length d
gives 58 hops, mean 2.09 pc, longest 3.44 pc, stretch 1.20 in the solar patch
(sparse patch: 28, 5.44, 9.52, 1.23); d squared gives 104 hops, longest 2.96,
stretch 1.41; d + 3 pc per hop gives 52 hops, stretch 1.25. The graph already
caps hop length, so the cost mostly changes the hop count. Use length; any
cost c(d) >= d keeps the Euclidean heuristic admissible, which is why an
unknown-space multiplier is allowed and d squared is not.

### 3.4 Voxel traversal (Amanatides-Woo) is not needed

On this grid the axes are the layer planes (linear in t), the ring cylinders (a
quadratic, not monotone) and the slot half-planes (an angle, with a slot count
that changes by ring), so A-W needs two phases split at the closest approach to
the axis. `sectors_along_segment` already collects every cut, sorts, and looks
each piece up with `sector_address_at` (about 27 microseconds per cell; 12.6k
cells for a 40 kpc segment in 0.32 s [C]). A streaming version would only buy
early exit at the first unfilled cell; add a generator later only if the flag
on very long hops needs it.

## 4. The spatial index: use the sector grid

**No new index is needed for NAV.43 or NAV.10.** `sectors(ring_index,
layer_index, ring_slot_index)` has a unique index (`uq_sectors_address`),
`star_systems.sector_id` is indexed, and every phenomenon table has
`sector_id` and a center index. Reading the cells a ball touches, then the
rows in them, reads 1.3 to 1.7 times the rows in the answer at 50 pc, against
7 to 45 times for a composite (x, y, z) index, which prunes on its first
column only [C].

### 4.1 SPATIAL indexes cannot do 3D

Both servers index only 2D shapes (an R-tree over the bounding rectangle) and
their `ST_` functions work in X and Y [R for MySQL 8.4 and MariaDB 11.4].
Measured on MariaDB 10.11: `POINT Z (1 2 3)` does not parse (NULL) and
`ST_Distance` of two `POINT Z` values is NULL [C]. A 2D index on (x, y) plus a
z filter read 6 times the answer for a 100 pc slab (126,150 entries for 19,873
rows at 50 pc, solar patch) and took 25.4 s to build against 3.8 to 6.6 s for
each B-tree at 1e6 rows [C]. It would also duplicate X and Y in galaxy-frame
points that must follow `advance_galactic_positions` and sector refiling. Do
not use it.

### 4.2 Alternatives, measured

"Everything within R pc", MariaDB 10.11, server CPU and index entries read, at
1e7 systems (1.56M sectors), R = 50 [C]. A is the composite (x, y, z) bounding
box; B address ranges OR'd; C the address join through `sectors` then
`sector_id`; D2 Morton (Z-order) ranges at a 2 pc leaf in a temporary table.

| Case (rows in answer) | A bbox | B OR ranges | C sectors join | D2 Morton temp join |
|---|---|---|---|---|
| core (1,345,778) | 7.7 s, 2,498,698 | 17.6 s, 1,542,978 | 8.2 s, 1,559,464 | 7.8 s, 1,437,298 |
| solar (55,406) | 1.61 s, 1,150,105 | 0.77 s, 67,618 | 0.43 s, 87,266 | 0.34 s, 66,484 |
| sparse (1,830) | 160 ms, 160,881 | 30 ms, 2,901 | 30 ms, 6,225 | 30 ms, 8,250 |

At R = 10, 1e7, solar (518 rows): A read 230,321 entries, B 1,276. Morton at an
8 pc leaf OR'd was no better than B.

- The composite index is the worst at scale (45 times the answer at R = 10,
  sparse) and the only option needing a new index.
- The cell approach stays within 1.3 to 1.7 times the answer and is flat as
  the table grows. Morton needs about 3,000 ranges at 2 pc (830 at 4 pc, 235
  at 8 pc) and a new column; its one gain, fewer rows at small R (1.5x
  against 2.4x volume overscan at R = 10), does not pay for that. A Hilbert
  decomposition was not benchmarked and would gain little [R].
- Each row costs about 6 microseconds of server CPU for the exact-distance
  test, so 1.3 million rows cost 8 s whatever the index. That is why NAV.43
  caps totals and takes page 1 from shells (section 5).
- **Do not send hundreds of OR'd ranges.** 3,000 Morton ranges as OR cost
  11.5 s of server CPU against 0.95 s joined from a temporary table; 560
  address ranges as OR cost 2.0 s against 1.1 s for a join through `sectors`
  [C]. MySQL falls back to a scan past `range_optimizer_max_mem_size` (8 MiB
  default) [R, untested]. Keep range lists under about 200 terms or use a
  temporary-table join, which both servers support.
- **Optional covering index** `star_systems(sector_id, position_x_mpc,
  position_y_mpc, position_z_mpc)` made the same query 1.6 to 2.1 times
  cheaper at about 35% of the table's other indexes' size [C]. Add it only if
  profiling on the production database asks for it.
### 4.3 The padding rule

`galaxyGeometry.enumerate_sectors_within_radius` returns the cells whose
**center** is within the radius, not the cells the sphere touches. For a 10 pc
sphere, 21.5% of the points inside it sit in cells the function does not list;
at 50 pc it is 3.1% [C]. The farthest corner of a cell from its center is
4.0 pc in ring 0 and 3.5 pc elsewhere at a 4 pc edge [C]. The rule for every
caller: **ask for `radius + one edge` (4 pc) and filter by exact distance.**
At R = 50 the unpadded call lists 8,060 cells in 14 to 25 ms; the padded one
about 10,000 cells in 31 ms [C]. `store.sectors_reached_by` already pads, by
the nominal half diagonal (`_sector_half_diagonal_pc`, 3.46 pc), a thin margin
against 3.5 pc in the slotted rings; pad by one edge there too. The
enumerator's docstring should say it lists centers, and it yields in no useful
order, so sort before applying any limit
(see fill-order-curves-and-core.md, section 3.3).

### 4.4 In-memory caches per request

cKDTree, one thread [C]: build 0.05 s (1e5 points), 0.78 s (1e6, 21 MB extra),
12.5 s (1e7, 256 MB); k = 6 for all points 0.43 s, 8.6 s and 111.7 s; a ball
of R = 50 returns in 0.9 to 5.3 ms in the solar patch and 2 to 127 ms in the
core. A whole-galaxy graph per request at 1e7 needs about 1 GB for neighbour
arrays alone, so a corridor is the way. `pykdtree` 1.4.3 gave identical
neighbours and was no faster (10.3 s against 7.6 s at 1e6) and is
LGPL-3.0-or-later [S: pypi.org/pypi/pykdtree/json]; `hnswlib` 0.8.0 is
approximate with no radius query [S]. Stay with cKDTree (`workers=` and
`return_length=` in the pinned scipy 1.13 [R]). Avoid `np.unique` on large
integer arrays: 9.4 s on 6e6 int64 on numpy 2.5.3 against 0.11 s for
sort-plus-diff [C].

## 5. Everything within N pc of a place (NAV.43, NAV.44, NAV.45)

`queryDb.objects_within` and `GET /api/near`, the "What's nearby" page and the
map actions.

**As built (NAV.43 and NAV.44, PR #804, `planetgen/db/near.py`).** The search
is `near.objects_within(conn, place, distance_pc, kinds, limit, offset)` behind
`GET /api/near`; a place is an object reference or a point in parsecs, the
limit is 50 pc, pages default to 50 rows (largest 200) and use `limit` and
`offset`. It reads the sectors the sphere can reach with one bounding-box query
on `idx_sectors_center`, measures exact distances in the galaxy frame, loads
only per-system body counts for the whole sphere and the bodies of the systems
on the requested page, and counts ungenerated sectors without filling them.
The shell-by-shell page 1, keyset cursor and "300+" capped totals described
below are not built; use them if a measurement at the core of a full galaxy
shows the bounding-box read too slow (this section's 1.3 million systems
figure is the case to test). The design below is the research behind the
choice.

1. **Cells.** `enumerate_sectors_within_radius(center, R + 4.0, edge)`, about
   10,000 cells at R = 50. Drop cells not in the cached generated set
   (section 2); their count is "sectors not generated yet".
2. **Rows per kind** by the cell list in batches of about 48 cells, ordered by
   center distance, through `sectors` (systems, rogue objects, facilities) or
   by `sector_id` (phenomena; extended ones also need the bounding-box test
   `_placed_phenomenon_rows` already does), using the address ranges or the
   temporary-table join of section 4.2, never a long id list.
3. **Page 1 from shells.** Read cells in batches by center distance; stop
   when 50 rows are held and the 50th best distance is under the next batch's
   center distance minus 4 pc. Page 1 cost 10 to 30 ms of server CPU where
   reading the whole ball cost 30 to 740 ms at 1e6 systems [C]. Later pages use
   a keyset cursor `(distance, kind, id)`, not OFFSET. Bodies (stars, planets,
   moons, belts, comets) take their system's distance, so rank systems and
   expand.
4. **Capped totals.** A 50 pc sphere at the core holds 1.3 million systems
   (8 s of server CPU just to count). Totals are exact under a cap, else
   "300+", as `SEARCH_COUNT_CAP` (300, `db/query.py`) does. The 50 pc limit
   stays; the row count, not the radius, needs the cap.
5. **Not generated yet.** The address indexes give for free the count (and an
   optional list) of pre-placed bright stars and scattered black holes and
   neutron stars in unfilled cells inside the sphere
   (`bright_stars.star_system_id IS NULL`, `phenomenon_scatter.built_at IS
   NULL`).

NAV.44 and NAV.45 add no routing content; they show "300+" past the cap, the
"not generated yet" counts and the kind filter.

**NAV.9** (references from search and locate): `galaxy_locate` returns sectors
and systems only and has no `ref` or parents. Build `ref` with
`object_ref.format(kind, id)` and take the parent chain from the JOIN each
panel already runs (planets and moons: `star_systems` then `sectors`), not
`resolve_object` per row, which makes about four queries per object. Names are
FULLTEXT-indexed on sectors, systems, stars, planets and moons; phenomena have
only a plain `KEY (name)`, so prefix matching works and substring matching
scans those small tables.

## 6. Near a line segment (NAV.10, NAV.25, NAV.6)

- **In memory:** a chain of spheres of radius sqrt(w^2 + (step/2)^2) centered
  every `step = sqrt(2) w` covers the cylinder of half-width w with a fetched
  to cylinder volume ratio of 1.73, the minimum; then filter by exact
  point-to-segment distance. cKDTree at 1e6, L = 1000 pc, w = 20 pc: 1.4 ms
  (1,491 points, solar) and 26 ms (52,152, core); a full-array bounding-box
  scan costs 33 to 70 ms at 1e6 and grows with N [C].
- **In the database:** `sectors_along_segment` gives the exact cells (300
  random segments checked against dense sampling: none missing, 43
  corner-touching extras in 19,619 cells [C]). Dilate with `neighbor_addresses`
  for w up to one edge; for larger w drop cells whose center is further than
  `w + 4 pc` from the segment (the one-ring dilation over-fetches 4.4 times at
  w = 4 pc). Intersect with the cached generated set, then join through a
  temporary table.
- **End to end, 1e7 systems** [C]: solar, L = 300 pc, w = 12: 3,253 cells,
  1,183 rows in 20 ms server CPU plus 15 ms client, route 3 ms; inner disk,
  L = 100: 10,319 rows, 60 + 44 ms, route 6 ms; core, L = 60, w = 4: 52,513
  rows, 340 + 222 ms, route 34 ms. This is the "measured time" NAV.10 asks for
  (on 1.56 million sectors rather than the TODO's 100,000).
- NAV.25 and NAV.6 use this same query. [course-avoidance.md](course-avoidance.md)
  bends each hop one at a time and does not change which systems are the stops.

## 7. Stops in unknown space (NAV.47)

Candidates are `bright_stars` (address index; absolute `position_*_mpc`;
`star_system_id` set once the sector is filled) and `phenomenon_scatter`
(address index only; `kind`, absolute position, `velocity_*_kms` for
hypervelocity stars, `built_at` set once filled). Planetary nebulae and
remnants are in the table but not on Boss's list. The scatter table holds about
5e8 rows (tens of GB; the repo's own comment puts runaway stars at 1.3 billion),
1,000 to 10,000 of them hypervelocity stars [C]. Stone density in the solar disk
is about 5e-4 per pc^3 for neutron stars and 1e-4 for black holes, plus bright
stars (about 6e7 galaxy-wide at 500 solar luminosities [R]); the simulation
uses 6e-4.

**Algorithm (refine flagged hops, lazily):**

1. Route as in section 3 on generated systems only; a flagged hop is a gap.
2. For each gap (u, v) of length L, take the unfilled cells in the cylinder of
   half-width w = max(8 pc, 0.1 L) (section 6). Read `bright_stars` rows with
   `star_system_id IS NULL` and `phenomenon_scatter` rows with `built_at IS
   NULL` and `kind IN ('black-hole', 'neutron-star', 'quasar',
   'hypervelocity-star')` in those cells, via a temporary table of cells.
3. Build a reach graph over {u, v} plus the stones, linking pairs closer than
   `r = 2 rho^(-1/3)` (rho the stone density measured in the corridor); double
   `r` until u and v connect or `r` reaches L. A* with cost = length; keep the
   direct edge in reserve at cost L (1 + unknown_penalty).
4. Replace the gap by the stone chain if one exists, marking each stone with its
   kind; otherwise keep the single unknown jump.

Simulation (Poisson stones, 30 trials, cost = length) [C]: at 6e-4 per pc^3,
D = 100 pc gives 5 hops (longest 23.5 pc, route 1.04 of D) and D = 500 pc 28
hops; outer disk 1e-4, D = 500 pc: 16 hops, longest 43 pc; halo-like 1e-5:
D = 100 pc has no usable chain (the single jump stays), D = 500 pc gives 8
hops, longest 90 pc. A d squared cost gave more, shorter hops (75, 16.6 pc, 1.31
of D at 6e-4), so cost stays length and the reach sets the hop size.

**Rules.** Stones apply to unfilled cells only: once a sector is filled its
scatter rows are built into `black_holes`, `neutron_stars` and `star_systems`,
so a filled cell's black hole is not a stop (it is still an obstacle for
[course-avoidance.md](course-avoidance.md)). A stone in an unfilled sector has
no sector row, so its hops are described in the Galactic Standard Frame
([navigation-frames.md](navigation-frames.md)).

**The `kind` lookup.** `phenomenon_scatter` has no `kind` index, so `WHERE kind
= 'hypervelocity-star'` scans about 5e8 rows. Either add `KEY (kind, ring_index,
layer_index, ring_slot_index)` (10 to 15 GB and a long online `ALTER`) or write
the hypervelocity and nucleus rows (a few thousand, inserted last by
`special_rows`) to a small `phenomenon_scatter_special` table read whole into
memory. Recommended: the side table. The corridor read of the big table is a
few hundred rows by the address index.

**Hypervelocity stars "at the time".** The sim clock is real elapsed time
(`orbit_simulation_state.last_updated_at`; `planetgen orbits` runs about
monthly). At 500 to 1,000 km/s a star moves 0.51 to 1.02 pc per 1,000 years
(1 km/s = 1.0227 pc per million years), 2.6e-3 pc a month, about 200 AU a year
[C]. A scatter row holds the plan-time position, so the position at time t is
`p0 + v (t - t_plan)`; `t_plan` is not stored, so add
`galaxy_shape.phenomenon_scatter_at` when the scatter is written. The effect is
under 0.01 pc a year, so a stone stays in its cell: the formula is for
correctness, not the index. Trap: `store.advance_galactic_positions` moves a
stored system by a rotation about the galactic axis (its galactic orbit) and
rotates its velocity vector, but never adds velocity times elapsed time, and
skips scatter rows, so a hypervelocity star in a filled sector orbits instead of
travelling outward. NAV.47 should define "position at the time" once and use it
for both.

## 8. Charting the sectors that block a course (NAV.48)

- **Bypass test.** Run the corridor search with flagged edges at infinite cost
  (flag evaluated lazily, memoised). If a route exists through known space,
  report its length against the best route and do not offer charting (or offer
  it as optional). If none exists, the flagged hops are the blockers: take
  `sectors_along_segment` of each, drop filled cells and cells outside
  `GalaxyBounds`; that is the list to chart.
- **Which cells.** Boss's wording is "all sectors between the two points".
  Charting every cell on the A-B line is far more than needed when the route
  already hops through known areas, so the default is the cells of the unknown
  hops only, then re-plot (new systems change the graph). A one-cell border
  adds about 3 to 9 times the cells; make it a checkbox.
- **Cost, before offering.** `generation.stats.estimate(sectors, stats,
  database, workers)` takes `(density, expected_systems)` per sector and returns
  seconds, bytes and a disk refusal; `run_galaxy.generate_sector_neighborhood`
  already does this for a sphere (`_neighborhood_candidates`,
  `_neighborhood_batch`, `estimate_only`), and NAV.48 needs the same with an
  explicit address list. A line in the disk plane crosses 0.26 to 0.32 cells
  per parsec, a random 3D direction about 0.4 [C]. At 6.3 systems per
  solar-density cell that is about 6 CPU-minutes per 1,000 pc at the PERF.3
  default (`DEFAULT_SECONDS_PER_SYSTEM` = 0.2 s) and about 30 at one second per
  system, the one-off system page's figure (`architecture.md`, flow 4); a core
  cell holds up to 241 systems; a cross-galaxy line (about 7,500 cells) is 14
  hours at solar density and one second per system, far worse through the
  bulge. The existing "ask again past about 5,000 sectors" confirmation is the
  right gate, with the estimate (time, disk, refusal) shown first.
- **Also needed:** the bright-star backfill (`backfill_bright_stars`, GEN.30)
  runs once around a center for neighborhoods; a corridor needs it per block or
  the dimmer tiers are missing along the line. The job must call it.

## 9. Travel times for the route (NAV.11, phase 1)

The route gets the same warp and fold tables as the direct distance (those in
[navigation-frames.md](navigation-frames.md), which agree with the project's
`travel-speeds/speeds.md`), per hop (unknown-space jumps included) and in
total, the total assuming a stop at every system on the route (Boss,
2026-10-02 04:19Z).

| Option | Meaning | Comment |
|---|---|---|
| A (default) | A stop is each system on the route; each hop timed rest to rest at the chosen factor; no stay | Sum of `direct.distance_ly / speed`; matches the constant speeds and no acceleration model |
| B | A plus an optional per-stop "stay" (default 0) added once per intermediate stop | One number on the page and API; models refuelling without a rule change |
| C | Visit every planet of every system on the route | Not worth modelling: 1 AU is 499 s of light-travel; 30 AU is 4.2 h at warp 1 and 1.2 s at warp 9.995 [C] |

Ship A with B as an optional parameter; "each planet" reads as "each stop on
the route". The total says how many hops are unknown-space jumps and the total
including stays. Worked example [C]: a 100 pc route in the solar disk has 58
hops of 2.09 pc mean, 121 pc (395 ly) in all; at warp 9.995 (33.41 ly a day)
that is 11.8 days, at warp 1 about 395 years, against 9.8 days for the direct
326 ly. Per-hop times go in a collapsible detail, the total in the summary.

## 10. Showing the route

- **A readable course map (NAV.41, phase 0, bug).** A wide panel across the
  usable width, labels at least body-text size, the route list below (Boss,
  2026-10-02 04:19Z). One 6,504 ly hop shrinks the other nine to a dot at one
  scale; a flagged hop could be drawn broken and dense parts inset (NAV.5,
  NAV.20).
- **Course and distance per stop (NAV.42, phase 1).** Each stop shows the
  course to the next in the existing notation, for example `045 mark 012, 3.2
  ly`, worked out in that hop's frame; a stone stop in an unfilled sector uses
  the Galactic frame.
- **Horizontal, wrapping (UX.35, phase 1, alongside NAV.12).** Keep
  `<ol class="nav-route">` and add `role="list"` (Safari drops list semantics
  from `list-style: none` [R]). Each `<li>` is one unit, the stop link plus the
  hop that follows (distance, then NAV.42's course text), `inline-flex` and
  `white-space: nowrap` so a stop is never split, `overflow-wrap: anywhere` for a
  very long name. The container is a wrapping flex row (`gap: 0.25rem 0.75rem`)
  in a wrapper with `container-type: inline-size; container-name: nav-route`,
  following the project's named container queries (`nav-row` in
  `html/static/style.css`): below about 30rem the hop line goes under the stop;
  below about 20rem fall back to the vertical list. Subgrid is not needed.
- **Long routes.** 50 to 100 or more stops are normal [C]. Collapse the middle
  in a `<details>` ("show 92 more stops"), always showing the first and last
  three stops, the longest hop and every unknown-space hop. Test
  (`test_web_a11y.py`) at 320 px with a 60-stop route and a 40-character name.
- **Red and glowing (NAV.36, phase 2).** A flagged hop is red with a glow in the
  route strip, on the NAV page's map and on the Galaxy Map course, with a legend
  entry; under `prefers-reduced-motion` red without the pulse; readable in both
  themes; the list also labels it "unknown space" in text.
- **Saved courses (NAV.39, phase 2).** A saved course (NAV.17) keeps which hops
  were unknown-space jumps; opening it checks again, since sectors may have
  been filled, and says what changed.
- **Waypoints (planned 2026-10-07, NAV.49):** objects picked in Star select mode
  become waypoints marked at every map level; two or more plot a course, which
  stays drawn until cleared.

## 11. Order

| ID | Piece | Phase | Needs |
|---|---|---|---|
| TEST.79 | Route edge cases, written first | 0 | **Built, PR #427** |
| NAV.34 | Join the route graph's islands (bug) | 0 | **Built, PR #427** |
| NAV.38 | Every sector a straight line passes through | 0 | **Built, PR #357** |
| NAV.10 | Corridor search, A*, the segment query | 1 | **Built** |
| NAV.12 | No hop limit, longest hop shown, unknown-space flag | 1 | NAV.34, NAV.38, TEST.79, NAV.10 |
| UX.35 | Horizontal, wrapping route strip | 1 | alongside NAV.12 |
| NAV.11 | Travel times per hop and for the route | 1 | NAV.10, NAV.12 |
| NAV.42 | Course and distance per stop | 1 | UX.35 |
| NAV.43 to 45 | Nearby search: query, page, map actions | 1 | NAV.7; NAV.43 then 44 then 45 |
| NAV.47 | Stones in unknown-space hops | 1 | NAV.12 |
| NAV.48 | Chart the cells that block a course | 1 | NAV.12 |
| NAV.36 | Unknown-space jumps red and glowing | 2 | NAV.12, UX.35, NAV.20 |
| NAV.39 | Saved courses re-check their unknown-space jumps | 2 | NAV.12, NAV.17 |

## 12. Why it works this way

- **No cap.** With most of the galaxy ungenerated, every capped test in the
  study that crossed ungenerated space found no route, at every cap from 5 ly
  to 5,000 ly, and the answer would change as sectors are generated. At top
  speed even the longest hop takes years (97,236 ly in 7.97 years at warp
  9.995, 4.66 at fold 8.5), so no hop is impossible under the current rules.
- **Why the route has stops at all.** With straight-line costs the direct edge
  is always shortest; stops exist only because the graph leaves long edges
  out. The 6-nearest rule is today's implicit limit (about 8 ly at the Sun's
  radius, 34 ly 1 kpc above the disk).
- **Flag instead of forbid.** A long hop across ungenerated space hides a gap
  behind a list that looks like a chain of nearby stars, and makes NAV.6's
  steering costlier. Flagging tells the user without refusing the course.

## 13. What carries over from the uploaded documents, and what does not

- **`Orbital Position and Vector Update Algorithms.md`, "Final Design of Vector
  Generation" (`resolve_neighbor_sectors`):** the idea carries over. A point
  becomes a ring, azimuth and column address, neighbouring rings' azimuth
  windows come from an angular buffer, and the records of the 25 to 33
  overlapping cells are merged. The padded enumerator of section 4.3 is the
  same dilation (one cell edge in place of the document's `r_search / r`).
  Superseded: its constants (11.5 ly thickness and arc, slot counts in
  multiples of 6, "10 ly cubic sector" early on). The shipped grid has a 4 pc
  edge (13.05 ly) and the master-wedge slot rule of `ring_sector_count` (3, 9,
  15, 21, 27, 36, ...; Boss, 2026-09-30); lookups use `sector_address_at` and
  `neighbor_addresses` (6 to 8 face neighbours per cell, measured in
  fill-order-curves-and-core.md, section 2.1).
- **Spatial hash and Morton keys:** that document has no Morton text. The
  proposal (4 pc cubic cells, a 64-bit Morton code, a 27-cell stencil,
  "replacing cylindrical geometries") is in `Computational Astrodynamics.md`,
  "Synthesis and Structural Recommendations". It was turned down for the
  shipped grid (`galaxy-drilldown-navigation.md`; `orbital-updates.md` section
  11), and fill-order-curves-and-core.md shows a cubic index cannot address the
  ring, layer, slot cells. Section 4.2 measures Morton ranges as a nearby-search
  key and finds no gain over the address index. Superseded; add no Morton
  column. The octree depth and point-mass table material belongs to
  `orbital-updates.md`.
- **The study's range cap and "deep-space hop over 2 sector edges"** are
  superseded by Boss's 2026-10-02 01:53Z decision. Its `lib/navmap.py` is now
  `planetgen/web/maps/navmap.py`. Its travel-time table agrees with
  `travel-speeds/speeds.md` and navigation-frames.md (6.7 ly takes 4.8 hours at
  warp 9.995).

## Evidence notes

[R] items to verify when paper and vendor access is allowed (the research
could only read search-result text):

- MySQL 8.4: spatial index needs NOT NULL and an SRID; WKT with Z rejected;
  `ST_` functions 2D. MariaDB 11.4: same 2D limit (measured only on 10.11).
- MySQL `range_optimizer_max_mem_size` behaviour with long OR/IN lists
  (MySQL 8.4 and MariaDB 11.4 were not installable; only 10.11 was measured).
- Generated-column indexes in both engines; Hilbert versus Morton range
  counts; scipy 1.13 `workers=` and `return_length=`; Safari `role="list"`.
- The bright-star total (60 million at 500 solar luminosities) is from
  `galaxy-coordinate-system.md`, not recomputed.

[C] caveats: server CPU resolution is 10 ms; the 1e6 and 1e7 databases are
synthetic (four density patches, half uniform and half Gaussian clumps of
sigma 4 pc; 164,064 sectors at 1e6, 1,563,383 at 1e7), not generated systems,
and the 2D-slab amplification used +-150 pc patches (a real disk is thicker);
MariaDB 10.11.14 only; pure-Python routing timings are the repo's own code
under load; A* and corridor results are one archipelago layout (seed 99).

## Sources

- Repository (checkout at c30d74a, rechecked at 27220d8): `docs/TODO.md`,
  `navigation-frames.md`, `orbital-updates.md`, `galaxy-coordinate-system.md`,
  `architecture.md`; `galaxy/nav_graph.py`, `galaxy/geometry.py`,
  `db/query.py`, `db/store.py`, `db/schema.sql`, `generation/stats.py`,
  `run_galaxy.py`, `phenomenon_scatter.py`, `tuning.py`,
  `src/tests/test_route_edge_cases.py`.
- Project files: `nav-hop-length/report.md`, `travel-speeds/speeds.md`,
  `nav-map-preplan/plan.md`; uploaded `Orbital Position and Vector Update
  Algorithms.md` and `Computational Astrodynamics.md`.
- PyPI JSON: pypi.org/pypi/pykdtree/json (1.4.3), pypi.org/pypi/hnswlib/json
  (0.8.0, 2023-12-03), pypi.org/pypi/scipy/json (lock file pins 1.13.1).
- One web search (MySQL 8.4 spatial, 2D index, Z): oneuptime.com and
  PlanetScale posts and the dev.mysql.com 8.4 spatial functions page described
  a 2D R-tree but not Z.
- Benchmark scripts (research scratchpad, not in the repository): `gen_data.py`,
  `gen_arch.py`, `bench_db.py`, `routing_bench2.py`, `routing_bench3.py`,
  `join_scipy.py`, `stones_sim2.py`.
