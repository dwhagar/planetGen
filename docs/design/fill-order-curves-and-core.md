# Fill order, region selection and the galactic core

Which order should the generator fill sectors in when it fills a ball, a
disc, a span or the galactic core; whether the "pruned Hilbert curve with
face-adjacent unit steps" wording of GEN.101 can be met; how to count and
enumerate regions cheaply (ADM.29, ADM.30, GEN.97, NAV.43); and what the
GEN.24 core fill costs and what it contains. The grid it works on is
described in [galaxy-coordinate-system.md](galaxy-coordinate-system.md) and
the density it fills with in [galaxy-disk-density.md](galaxy-disk-density.md).

Informs: GEN.101, ADM.29, ADM.30, GEN.24, GEN.23, GEN.97, ADM.31, PERF.18, NAV.43, PERF.3

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [C] computed in this research and re-runnable (scripts named
under Sources), [S] read directly from a package page or API, [R] recalled
from memory and unconfirmed. The research environment could not read papers
and its web search budget was already spent, so every astronomy number and
every paper citation is [R] and is listed under Evidence notes. Everything
about the repository and every experiment is [C].

## Decisions already taken

- **GEN.101 fill order (Boss, 2026-10-07 11:47Z):** fill sectors nearest the
  original sector first, circling outward, using a pruned Hilbert octree
  over the integer lattice, deterministic and O(N), "face-adjacent unit
  steps ... visiting every valid voxel centroid exactly once with zero
  backtracking".
- **GEN.101 (Boss, 2026-10-07 17:11Z):** keep the Hilbert order, with a
  logged jump where the ball cuts the curve.
- **GEN.24 core (Boss, 2026-10-01 14:55Z):** the admin chooses the size, as a
  radius or a number of rings, on the Generate page and the CLI; suggested
  default the bulge scale radius (200 pc, about 50 rings, about 7,900
  sectors); layer 0 only; the generate-around sphere and bright-star
  backfill (GEN.23) run around the core as around any generated sector; fill
  from ring 0 outward so a run that stops early leaves a solid disc.
- **ADM.30 (Boss, 2026-10-07 11:47Z):** "x sectors radius (1 = minimum for
  contiguous orthogonal connection between each sector and its adjacent
  sectors)".

## Findings in short

1. **The unit-step, visit-once, zero-backtrack guarantee cannot hold for a
   pruned ball on any curve, Hilbert included.** It is a parity obstruction
   (section 1.1): no integer radius from 1 to 40 is walkable, and radius 1
   needs at least 4 jumps. Boss's recorded decision (keep Hilbert, log the
   jump) is the only wording that can be true, and it means "the order is
   face-adjacent for about 99% of steps and every other step is logged with
   its length", not a guarantee. The NP-completeness remark in the GEN.101
   quote is correct in general but is not the reason; balls fail for a
   simpler one.
2. **Plain pruned Hilbert also fails "nearest first"**: it starts at a corner
   of the enclosing cube, 0.99 R from the centre (section 1.2). "From the
   origin sector outward" and "face-adjacent" cannot both come from the
   Hilbert index.
3. **The grid is not a voxel lattice** (6 to 8 face neighbours per cell), so
   a Hilbert walk over voxels cannot address sectors (section 2).
4. **Recommended order:** ring serpentine for the core and any axis-centred
   disc or cylinder (zero jumps, nearest first); deterministic greedy
   nearest-first walk elsewhere (96% to 99.5% face-adjacent, nearest-first
   quality about 0.9); Hilbert only as a sort key or block tie-break
   (section 2.4).
5. **The GEN.24 default of 200 pc is 7,857 sectors and about 7.45 million
   systems**: 17 days on one worker, about 400 GB, refused by PERF.3 on any
   disk under about 1.6 TB. Recommend 50 pc (section 4.3). The "bulge scale
   radius 200 pc" it cites is a revision-2 value; the shipped default is
   1,580 pc.
6. **Two traps for ADM.30** (section 3.3): a Euclidean disc of radius 1 edge
   returns 3 cells, not the point and its face neighbours, and the
   enumerator's order is ring-major, so "nearest first as enumerated" in
   `_neighborhood_batch` is false.

## 1. Unit-step fills on a ball of lattice points

### 1.1 The parity obstruction

The cubic lattice graph is bipartite. Colour a point by (x + y + z) mod 2;
every unit step changes colour, so a path alternates colours and its two
colour counts differ by at most 1. Take a ball with e even points and o odd
points, imbalance D = |e - o|. Any cover of the ball by p vertex-disjoint
unit-step paths has D <= p, so a single walk with jumps in between needs
at least D - 1 jumps. A second necessary condition: a point with one
neighbour in the set must be a path end, and a path has two ends.

Balls centred on a lattice point, centre counted as even [C, `exp1.py`,
`exp1c.py`]:

| R | points N | even | odd | imbalance D | jumps forced (D - 1 or more) |
|---|---|---|---|---|---|
| 1 | 7 | 1 | 6 | 5 | 4 |
| 2 | 33 | 19 | 14 | 5 | 4 |
| 4 | 257 | 141 | 116 | 25 | 24 |
| 6 | 925 | 459 | 466 | 7 | 6 |
| 13 | 9,171 | 4,585 | 4,586 | 1 | 0 by parity, but 6 degree-1 tips, so at least 2 |
| 20 | 33,401 | 16,757 | 16,644 | 113 | 112 |
| 40 | 267,761 | | | 257 | 256 |
| 60 | 904,089 | 452,261 | 451,828 | 433 | 432 |

The rows for R = 1, 2, 4, 6, 13 and 20 were recomputed while writing this
document and match.

Among all 1,335 distinct balls with R squared at most 1,600, 16 pass the
parity test (R squared = 3, 30, 48, 169, 177, 204, 214, 217, 334, 432, 745,
802, 1144, 1150, 1203, 1218). Only R squared = 3 (27 points) was confirmed
walkable, by a depth-first search; the others pass a necessary condition
and were not searched. No integer R from 1 to 40 is walkable (R = 13 fails
on the degree-1 tips). So walkable pruned balls are the exception, and for
the integer radius in sectors that an admin would type, essentially never.

Theory, all [R]: Hamiltonian path in general grid graphs is NP-complete
(Itai, Papadimitriou, Szwarcfiter, SIAM J. Comput. 1982; Papadimitriou and
Vazirani 1984), solid grid graphs have a polynomial cycle algorithm (Umans
and Lenhart, FOCS 1997), and the 3D case is believed no easier. This does
not matter for the fix: parity already decides it.

### 1.2 The pruned Hilbert curve

Implementation: Skilling's transpose algorithm (Skilling, "Programming the
Hilbert curve", AIP Conf. Proc. 707, 2004 [R]; the alternative is Butz 1971
[R]), vectorised in numpy, about 35 lines (`hil.py`). Verified on full cubes
of side 2 to 32: bijective, every step has L1 length 1, start (0,0,0), end
on the same edge [C, `t1.py`]. A pure-Python version is about 25 lines. PyPI
has `hilbertcurve` 2.0.5 (MIT, numpy only, last release 2021-03-29) and
`numpy-hilbert-curve` 1.0.1 (MIT, 2020-11-07) [S]; both are old and small,
so vendoring the 25 lines is preferable to a dependency.

Setup: a ball of radius R centred on a lattice point, in a cube of side
2^b centred on it, b minimal with 2^(b-1) > R. The depth formula in the
GEN.101 quote, D = ceil(log2(2R)), is one level short when R is a power of
two (R = 1, 2, 4, 8, 16, 32): the ball is then 2R + 1 points across, which
does not fit in 2R. Order is the Hilbert index; points outside the ball are
skipped.

Results for R = 60 (N = 904,089); other radii behave alike [C, `exp1.py`,
`exp1b.py`, `exp1d.py`]:

| Order | non-unit steps | mean / median / max jump (voxels) | path per point | nearest-first Q at 1% / 10% / 50% | prefix reach, max radius over ideal at 1% / 10% |
|---|---|---|---|---|---|
| Hilbert, pruned | 8,432 (0.93%) | 2.3 / 1.4 / 32.8 | 1.01 | 0.00 / 0.12 / 0.49 | 4.66 / 2.16 |
| Morton (Z-order), pruned | 50.0% | 1.9 / 1.4 / 120 | 1.47 | 0.00 / 0.01 / 0.49 | 4.66 / 2.16 |
| z-major serpentine raster | 0.7% | 1.9 / 1.4 / 16 | 1.01 | 0.00 / 0.00 / 0.50 | 4.64 / 2.15 |
| Sorted by distance (ties by Hilbert) | 100% | 10.4 / 7.9 / 102 | 10.4 | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 |
| Shells 1 voxel wide, Hilbert inside | 37.0% | 2.3 / 1.4 / 95 | 1.49 | 0.99 / 0.99 / 0.99 | 1.01 / 1.01 |
| Shells 4 voxels wide, Hilbert inside | 9.3% | 2.4 / - / 98 | 1.13 | 0.83 / 0.98 / - | 1.24 / 1.01 |
| 4x4x4 blocks by block distance, Hilbert inside and between | 1.6% | 4.7 / - / 91 | 1.06 | 0.88 / 0.93 / - | 1.46 / 1.12 |
| Greedy nearest-first walk (window 256) | 1.0% | 16.0 / - / 109 | 1.15 | 0.90 / 0.91 / - | - / 1.27 |
| Lower bound on jumps from parity | at least 432 (0.05%) | | | | |

Q(n) is the fraction of the first n visited points that lie inside the
sphere containing the n nearest points. The first voxel of the pruned
Hilbert curve is at distance 59.6 of 60; the sorted orders start at the
centre.

Reading the table:

- Pruned Hilbert's jump count is about 2.35 r squared (240 at r = 10, 2,125
  at r = 30, 8,432 at r = 60), about a fifth of a jump per surface voxel;
  mean jump about 2.3 voxels at any size, with occasional long ones (up to
  54 at r = 40) where a whole sub-cube is cut away. Shifting the ball by a
  few voxels inside the cube changes the jump fraction by under 10% (0.009 at
  r = 60 in all six offsets tried), so alignment is not a lever.
- Morton and Gray-code-style curves (Faloutsos 1986 [R]) are never unit-step
  (half of all steps jump) and no better at nearest-first. Use them only as
  fast sort keys.
- The Moore curve (a closed Hilbert variant, 3D forms in Haverkort's work
  [R]) makes start and end adjacent, which is irrelevant to a ball. 3D Peano
  (3x3x3 serpentine recursion [R]) is unit-step on a cube of side 3^k but
  behaves like Hilbert once pruned. Gosper-style curves are 2D hexagonal;
  there is no standard 3D one.
- Locality theory for Hilbert (Moon, Jagadish, Faloutsos, Saltz, IEEE TKDE
  2001; Niedermeier, Reinhardt, Sanders 2002; Haverkort and van Walderveen
  2010 [R]) is about the cost of range queries, which is what a database
  wants when it stores sectors in curve order. It says nothing about
  visiting a ball nearest first. The "onion curve" (Xu, Nguyen, Tirthapura,
  ICDE 2018 [R]) is the published relative of the shell ordering tested here.
- Best-first search by minimum cell distance over an octree (Hjaltason and
  Samet's incremental nearest neighbour [R]) reproduces exact distance order:
  perfect Q, but 100% non-unit steps and a path of 10 steps per point at
  r = 60. Right for "what do we generate first if the run stops early",
  wrong for locality.

### 1.3 What this means for the GEN.101 wording

| Demand in the quote | Verdict |
|---|---|
| Deterministic and O(N) | Met by every order above (sorting is O(N log N); the greedy walk is O(N x window)). |
| Nearest first | Not met by plain pruned Hilbert; met by greedy, shells plus Hilbert, blocks, serpentine. |
| Locality | Met by greedy, shells plus Hilbert, blocks. |
| Unit steps, every point once, zero backtracking | Cannot be met for a ball, with or without Hilbert. The achievable statement is "about 99% of steps are face-adjacent; each other step is logged with its length". Parity says no algorithm can beat 0.05% jumps at radius 60, and at radius 1 four of the six steps must jump. |
| Begin at the centre and expand | Not possible for one Hilbert curve: its prefixes are unions of dyadic cubes. Placing the origin at the cube corner and filling one octant per curve starts at the origin, but the eight octants then run one after another, which is not nearest first either. |

## 2. The real grid: ring, layer, slot

Code read: `src/planetgen/galaxy/geometry.py` (the TODO items call it
`galaxyGeometry.py`), `src/planetgen/generation/run_galaxy.py`
(`_neighborhood_candidates`, `_neighborhood_batch`, `_generate_addresses`,
`run_shell`).

### 2.1 How irregular the neighbours are

Rings 0 to 59, layer 0 [C, `exp2b.py`]; the neighbour counts were recomputed
while writing this document and match.

| Property | Measured |
|---|---|
| Face neighbours per cell (`neighbor_addresses`, includes up and down) | 6 (4,272 cells), 7 (6,408), 8 (633) |
| In-plane face neighbours | 4, 5 or 6 |
| Neighbours of a ring-0 cell in ring 1 | 3 (3 slots against 9), not 1 or 2 |
| Neighbour relation symmetric and connected | yes, 0 asymmetric pairs |
| Distance to same-ring neighbour, in edges | 0.87 (ring 0, 3 slots) to 1.045; mean 1.00 |
| Distance to vertical neighbour | exactly 1.00 |
| Distance to neighbour in an adjacent ring | 1.00 to 1.38, mean 1.07 |
| Largest in-plane face-neighbour distance, rings below 2000 | 1.379 edges (ring 30 slot 174 to ring 29 slot 164) |
| In-plane face neighbours farther than 1 edge (rings below 200) | 33% |

"Face-adjacent unit step" therefore has no Euclidean meaning here: it means
"is in `neighbor_addresses`". Where the slot count changes between rings a
cell touches one or two cells in the next ring (three from ring 0 into ring
1), the brick-course offset the docstring describes. The slot-0 seam is not
special in the adjacency graph (slot N - 1 and slot 0 are neighbours modulo
N) and does not exist in Cartesian coordinates. It matters only to orders
that walk (ring, slot) as coordinates: a ring is a cycle, not a line. The
serpentine in section 2.4 turns where rings align.

### 2.2 Mapping the grid onto a cubic voxel lattice does not work

Putting each layer-0 cell centre of a 200 pc disc into a 4 pc cube: with
cubes cornered on the origin, 6,623 cubes get one cell and 617 get two; with
cubes centred on the origin, 6,168 get one, 840 get two and 3 get three.
Inside r <= 190 pc, 584 to 754 of about 7,090 cubes get no cell [C,
`exp2d.py`]. "One voxel = one sector" is false. A Hilbert index can only be
a sort key on quantised centre coordinates, never a lattice walk. Defining
the curve on (ring, layer, slot) index space instead (a compact Hilbert index
for unequal side lengths, Hamilton and Rau-Chaplin, Inf. Process. Lett. 2008
[R]) would put slot 0 and slot N - 1 at opposite ends and treat ring 3 slot 7
and ring 4 slot 7 as neighbours although they sit at different angles, so it
is not recommended.

### 2.3 Orders compared on real cells

"Face-adj" is the share of consecutive steps that are in
`neighbor_addresses`; path is the sum of centre distances in edges per cell;
Q as in section 1.2 [C, `exp2c.py`, `exp2f.py`].

| Region (cells) | Order | Face-adj | Mean / max jump (edges) | Path per cell | Q at 10% | First cell distance |
|---|---|---|---|---|---|---|
| 50 pc sphere on the axis (8,169) | Enumeration today | 97.0% | 2.2 / 24 | 1.03 | 0.36 | 48 pc |
| | Distance sort | 77.8% | 5.4 / 24 | 1.97 | 1.00 | 2 |
| | Hilbert key (centres quantised to 0.5 edge) | 90.1% | 1.9 / 10 | 1.12 | 0.10 | 49.5 |
| | Shells of 1 edge + Hilbert | 66.9% | 2.5 / 16 | 1.50 | 1.00 | 2 |
| | Greedy walk, window 256 | 99.5% | 5.5 / 19 | 1.03 | 0.59 | 2 |
| 50 pc sphere off the axis, ring 700 (8,121) | Enumeration today | 89.4% | 10.0 / 25 | 1.96 | 0.00 | 49.7 |
| | Distance sort | 0.1% | 12.2 / 25 | 12.2 | 1.00 | 0 |
| | Hilbert key | 94.9% | 2.1 / 9.8 | 1.06 | 0.11 | 49.0 |
| | Shells of 1 edge + Hilbert | 61.0% | 2.3 / 16 | 1.49 | 0.94 | 4 |
| | Greedy walk, window 256 | 96.1% | 3.2 / 11.8 | 1.09 | 0.90 | 0 |
| Core disc, layer 0, R 200 (7,857) | Enumeration today | 99.4% | 1.4 / 1.5 | 1.00 | 1.00 | 2 |
| | Hilbert key | 80.3% | 2.1 / 51 | 1.24 | 0.00 | 198 |
| | Greedy walk | 100% | 0 | 1.00 | 1.00 | 2 |
| Cylinder, 30 rings x 11 layers (30,987) | Ring serpentine | 100% | 0 | 1.00 | by cylindrical radius | ring 0 |

Plain distance sort is the worst for locality; off the axis its steps are
essentially never adjacent, because it hops between the two sides of the
centre as distances tie. The greedy walk is the only order that is both
near-perfect on adjacency and near-perfect on nearest-first. The Q = 0.59 at
10% on the axis comes from a 3D spiral that fills layers in a leapfrog way;
Q at 50% is 0.88. Cost: an O(N log N) sort plus an O(N) pass, about 10
microseconds per cell in Python for the neighbour lookups [C].

### 2.4 Recommended orders

**Ring serpentine (core, axis-centred discs and cylinders).** Walk ring 0
slots upward, ring 1 slots downward, ring 2 upward, and so on. Every ring's
slot 0 starts at angle 0 and its last slot ends at 360 degrees, so
consecutive rings always touch at the turn. Built for a 30-ring, 11-layer
cylinder (30,987 cells): 0 non-adjacent steps, every cell once [C,
`exp2d.py`]. On layer 0 the current ring-major enumeration is already
nearest-first with 49 near-jumps (the ring-to-ring wrap); the serpentine
removes them. For several layers, fill one ring through all its layers
(serpentine in the layer index) before the next ring. This is exactly the
GEN.24 rule "fill from ring 0 outward" and needs no Hilbert code.

**Greedy nearest-first walk (any other centre).** Deterministic; the order
depends only on the centre, the radius and the grid; no random numbers.

1. Take the region's cells (the enumerator, trimmed to the galaxy outline).
2. Key of a cell: (distance to the centre rounded to 1e-9 pc, ring, layer,
   slot). Sort the cells by key; start at the first.
3. At each step move to the unvisited member of `neighbor_addresses(current)`
   with the smallest key.
4. When no such neighbour exists, choose, among the next 256 unvisited cells
   in key order, the one nearest the current cell, ties by (distance, key).

On a 50 pc sphere (8,169 cells) this gives 99.5% face-adjacent steps on the
axis and 96% off it, nearest-first quality about 0.9, mean jump 3 to 5
edges, against Hilbert-key ordering at 90 to 95% adjacent but Q about 0.1.

**If Hilbert is kept (Boss's 2026-10-07 17:11Z decision).** Use it as a
tie-break key, not as the walk: world-aligned 4-voxel blocks ordered by
distance, Hilbert inside the blocks, gives 1.6% jumps, Q 0.93 at 10% and a
path of 1.06 steps per sector (table in section 1.2); the first sector is
the centre. Plain pruned Hilbert starts at the far edge. Log each non-adjacent
step with its length, and report the count and the mean at the end of the
run.

### 2.5 Fill order does not change results, only what exists when a run stops

A sector's contents depend only on the galaxy seed and its address
(`reproducible-galaxies.md`, section 3: "never on the order sectors run in
or the worker count"). `_submit_batch` queues all sectors and waits; with
more than one worker, completion order is not submission order. So the order
matters for: what `--limit N` keeps (ring, column, block and shell modes cut
the batch at N before submitting: `run_ring_batch`, `_generate_addresses`),
the state left by a cancelled or failed job, the progress display, and
database locality if sector ids follow submission order. It does not matter
for correctness. This is the reason to prefer "cheap, nearest first, logs
its jumps" over a curve with a theorem attached: with parallel workers "zero
backtracking" has no observable meaning.

## 3. Regions: spans, cylinders and sphere enumeration

### 3.1 Counting without enumerating

- Sectors per layer in rings 0 to K - 1 is a prefix sum of
  `ring_sector_count`; `GalaxyBounds` already keeps it (`_cumulative`). A
  ring span [a, b] on one layer costs `cum[b + 1] - cum[a]`; a layer span
  costs O(layers); a column span (ring, slot range) is the sum over its rings
  of the layers that reach each ring (`galaxy_column`). Prefix sums for 4,000
  rings take 0.02 s [C, `exp5.py`]. ADM.29's "estimate shown first" is
  therefore free; only PERF.3's per-sector density needs sampling (3.4).
- On-axis disc of radius rho on one layer: exactly the sum of N_i over rings
  with (i + 1/2) x e <= rho, about pi x floor(rho / e + 1/2) squared [C,
  `exp2e.py`, recomputed with `ring_sector_count`]:

  | rho | Rings | Sectors | Sectors if 11 layers all exist |
  |---|---|---|---|
  | 20 pc | 0..4 | 75 | 825 |
  | 50 pc | 0..12 | 531 | 5,841 |
  | 100 pc | 0..24 | 1,965 | 21,615 |
  | 200 pc | 0..49 | 7,857 | 86,427 |
  | 400 pc | 0..99 | 31,425 | 345,675 |

- Off-axis disc: about pi x (rho / e) squared. Measured actual over formula
  for 40 random centres per radius: rho = 8 pc 0.88 to 0.96, 20 pc 0.90 to
  0.99, 50 pc 0.99 to 1.08, 100 pc 0.97 to 1.00, 200 pc 0.99 to 1.00. Good
  for an estimate; enumerate when the exact list is needed.
- Ball of radius R: about (4/3) pi (R / e) cubed cells: 8,169 at 50 pc
  (formula 8,181), 65,181 at 100 pc, 524,835 at 200 pc [C].

### 3.2 Enumerating a cylinder

A cylinder enumerator is the sphere enumerator with the vertical offset
fixed. For each ring from ceil((r_c - rho) / e - 1/2) to floor((r_c + rho) /
e - 1/2), take one slot window from acos((r_m squared + r_c squared - rho
squared) / (2 r_m r_c)) (all slots if the disc contains the axis), tested
against the cell centre. A prototype (`disc_cells` in `exp2e.py`) matched
brute force on 60 random discs with 0 differences [C]. The existing
`enumerate_sectors_within_radius` already does this per (ring, layer) pair;
a cylinder is its loop with `planar` held at rho for every layer, a ten-line
variant. Each (ring, layer) meets the disc in one contiguous slot window
(wrapping at most once), so the unique index `uq_sectors_address (ring_index,
layer_index, ring_slot_index)` answers "which of these are generated" with
one range predicate per (ring, layer) instead of a list of ids.

### 3.3 ADM.30: what "radius x sectors" means

Boss: "1 = minimum for contiguous orthogonal connection between each sector
and its adjacent sectors". Two readings differ on this grid [C, `exp2b.py`,
`exp2g.py`]:

- **Euclidean disc of radius x edges** around the centre cell. At x = 1 it
  returns the centre plus 2 cells (3 in total), because most face neighbours
  sit 1.0 to 1.38 edges away. It fails Boss's definition. The face-neighbour
  set is inside the disc only from radius 1.38 edges up.
- **Face-adjacency steps** (breadth-first over in-plane `neighbor_addresses`):
  x = 1 gives 5 or 6 cells per layer (4 to 6 neighbours plus the centre);
  x = 2 gives 13 to 17; x = 3 gives 25 to 34, against pi x squared = 3.1,
  12.6, 28.3. A diamond-like shape, not a disc.

Recommendation: a Euclidean disc of radius (x + 0.385) edges. At x = 1 that
is 1.385 edges, which returned exactly the centre plus its in-plane face
neighbours in all 1,331 trials (rings 0..59, 100, 300, 1000; every neighbour
is within 1.379). Larger x gives a round shape connected by face adjacency
from the centre. Boss gave no height: default to x layers either side (a
can of 2x + 1 layers), so x = 1 includes the cells above and below, which
are face neighbours too, and make the half-height a separate field.

Second trap, in the existing code: `enumerate_sectors_within_radius` yields in
no particular order (ring-major in practice), so the `_neighborhood_batch`
docstring "nearest first as enumerated" is false. In an off-axis test the
first 10% of the list contained none of the nearest 10% [C, `exp2a.py`].
Sort before applying any limit.

### 3.4 Cost: time, memory, streaming, estimates

- Speed on this machine (Python 3.13): about 1.3 microseconds per cell for
  the enumerator, 14 for `relative_density`, 15 for `sector_address_at` [C].
  Sphere of 50 pc: 0.02 s; 100 pc: 0.07 to 0.15 s; 200 pc: 0.7 s. Python 3.9
  on users' machines may be 1.5 to 2 times slower [R]. Enumeration is never
  the bottleneck; generating one sector takes seconds.
- Memory if the list is materialised (`_neighborhood_candidates` builds
  one): about 170 bytes per cell [C: 1.4 MB for 8,169 cells], so 90 MB at
  200 pc and about 1.4 GB for the 8 million cells of a 500 pc sphere, before
  `sector_args` per sector, which is heavier. For runs above about 50,000
  cells, stream the enumerator in chunks (for example 1,000 cells) into the
  work queue, keeping only the key list (three small ints per cell) for the
  nearest-first sort, and build `sector_args` lazily. Sorting (key, address)
  pairs costs a few tens of MB at 500,000 cells.
- PERF.3 estimate: `relative_density` on every cell costs 14 microseconds x N,
  0.11 s for the 7,857-cell core but 14 s for a million cells. Sampling 200
  random cells estimated the 7.455 million-system core to within 1% (7.37 to
  7.51 million over 20 trials) [C, `exp5.py`]. Sample above about 20,000
  cells and label the figure an estimate.
- Warn and confirm: the repository uses `LARGE_RING_WARNING_THRESHOLD = 2000`
  sectors for ring, block and shell modes. Use one rule for all region modes
  (spans, cylinders, core, neighbourhoods): always show the PERF.3 estimate;
  require `--yes` or the page's confirm above 2,000 sectors; refuse to
  enumerate in a web request above about 2 million cells (a 200 pc sphere is
  0.5 million; a 400 pc sphere is 4.2 million and 6 s, fine in a worker, not
  in a request).

### 3.5 GEN.97, NAV.43, GEN.23, PERF.18, ADM.31

- **GEN.97 (N random neighbourhoods).** Expected cell count is N x (4/3) pi
  (r / e) cubed minus overlaps; 100 neighbourhoods of 12 pc is about 11,000
  cells. Deduplicate with a set of addresses while streaming; fill each
  neighbourhood in the greedy order and the neighbourhoods in draw order.
  `GalaxyBounds.random_address` already draws centres uniformly by volume.
- **NAV.43 (`objects_within`).** The enumerator names sectors whose centre is
  within r. A sector whose centre is outside but whose cell reaches the
  sphere would be missed; a cell extends at most about half a diagonal
  (about 3.7 pc at 4 pc: 1.06 x sqrt(3) / 2 x 4) past its centre. Query with
  r + 3.7 pc (or r + one edge for margin), then filter on exact object
  distance. For r = 50 pc that is about 10,000 cells; use one `(ring, layer,
  slot BETWEEN)` predicate per (ring, layer) pair (about 600) rather than a
  10,000-id IN list, or the `idx_sectors_center (center_x_pc, center_y_pc,
  center_z_pc)` bounding box with an exact distance test. Cells that were
  never generated simply return nothing.
- **GEN.23 and PERF.18 (backfill around the core).**
  `backfill_bright_stars_around` calls the enumerator once per generated
  centre over the largest tier radius (100 ly = 30.7 pc). For a 50 pc core
  that is 531 calls and 1.0 million enumerated cells to find 16,527 unique
  cells (5.2 s); for 100 pc, 1,965 calls, 3.7 million cells, 44,247 unique,
  16.8 s; extrapolated to 200 pc about 14.8 million enumerated cells (about
  70 s) for roughly 100,000 unique ones [C, `exp4.py`]. A 60-fold
  redundancy. For a convex run region it is cheaper to enumerate once around
  the region's bounding disc and keep cells whose distance to the region
  (here max(0, distance from axis - R_core) on layer 0) is within the tier
  radii, or to enumerate only the margin ring.
- **ADM.31 (show on the Galaxy Map).** The fill list is also the "what was
  made" list. Return a 500,000-cell result as ring, layer and slot-window
  ranges (one per (ring, layer)), not as a list of addresses, to keep the map
  request small.

## 4. The galactic core (GEN.24)

### 4.1 What the code does at the centre today [C, `exp3.py`]

- Defaults (`planetgen plan` in `src/planetgen/cli/generate.py`; the Generate
  page `generate_page.py`): thin disk scale length 2,600 pc, height 300 pc,
  **bulge scale radius 1,580 pc** (axes 1,580 x 620 x 430 pc), bulge
  amplitude 3.11, 2 arms, `k_norm` 36.6.
  [galaxy-disk-density.md](galaxy-disk-density.md) section 5 records that
  revision 2 had a 200 pc sphere and amplitude 1 and that GEN.118 replaced
  them with the COBE/DIRBE bar. The "bulge scale radius, 200 pc" in the GEN.24
  text (and in `orbital-updates.md` section 11, "bulge scale radius 200 pc")
  therefore describes a model that no longer exists. No other 200 pc constant
  remains in the density module; the nuclear stellar disc "inside 200 pc" is
  explicitly listed there as not modelled.
- Density is nearly flat through the core: relative density 154 at the
  centre, 140 to 169 by azimuth within 25 pc, 163 at 200 pc along the bar and
  132 across it. With `expected_system_count_at_density_1` = 6.31 systems per
  4 pc sector, a core sector holds about **880 to 1,060 systems** (14 to 17
  systems per cubic parsec, against 0.1 in the solar neighbourhood [R]). The
  model has no central cusp, no nuclear cluster and no nuclear disc.
- There is no per-sector star cap: `MAX_SYSTEMS_PER_SECTOR = 1e9` in
  `generation/stats.py` only limits the estimate text. Generation draws a
  Poisson count around the expected value (`SpaceSector.grow_from_seed`) and
  places each system by Poisson-disk growth with Hill-sphere separation.
- Hill spheres are not the constraint at that density.
  `calculate_hill_sphere(d, m, M_MW)` with the whole Milky Way mass (1.15e12
  Msun) gives about 27 AU for a 1 Msun star 2 pc from the centre and 1,360 AU
  at 100 pc [C], against a mean spacing of 0.25 pc (about 52,000 AU) at 16
  systems per cubic parsec. The formula uses the total galaxy mass, not the
  mass enclosed at the radius (about 1e9 Msun inside 200 pc [R]), so it
  understates the tidal radius by roughly a factor 10 near the core. Not a
  blocker; worth a comment in GEN.24 or the physics review.

### 4.2 Real numbers, all [R] and to be verified

| Structure | Size | Mass | Notes |
|---|---|---|---|
| Sgr A* | point | 4.3e6 Msun | distance 8.2 to 8.3 kpc |
| Nuclear star cluster | half-light radius about 4 pc (Schodel et al. 2014, A&A 566, A47; Feldmeier-Oelkers et al. 2014) | 2.0 to 2.5e7 Msun | mean density inside 4 pc about 4.5e4 Msun/pc^3 (derived: 1.2e7 Msun in 268 pc^3) |
| Nuclear stellar disc | radius about 100 to 230 pc, scale height about 30 to 45 pc (Launhardt et al. 2002; Sormani et al. 2022) | about 1e9 Msun | mean density inside 200 pc about 50 to 150 Msun/pc^3 (derived, assuming 0.7e9 in pi x 200^2 x 90 pc^3) |
| Central molecular zone | about 200 pc radius (about 1.5 degrees of longitude) | gas 3 to 5e7 Msun | gas, not stars; a ring or twisted ring near 100 pc plus clouds |
| Bar and bulge (boxy) | half-length about 2 kpc; long bar to about 4 to 5 kpc | about 1.5e10 Msun | the model's 31% bulge share |

Star densities: the nuclear stellar disc volume holds roughly 100 to 300
stars per cubic parsec (about 0.5 Msun per star averaged over a population
with many low-mass stars [R]), so about **6,000 to 20,000 stars per 4 pc
sector**. The generator's 880 to 1,060 systems (about 1.3 stars each by the
PERF.3 default `_stars_per_system` = 1.3, so 1,100 to 1,400 stars) is 5 to
15 times lower. In the nuclear cluster the density rises further: about
1e5 Msun/pc^3 near 1 pc and 1e6 to 1e7 Msun/pc^3 within 0.1 pc (derived
from an r^-1.4 slope normalised to the mean above; low confidence, check
Genzel et al. 2010 and Schodel et al. 2018 [R]). The sector holding Sgr A*
(a 4 pc cube, r of about 2 to 2.5 pc) contains of order 1e6 to 1e7 stars, a
thousand times anything generated here, and cannot be filled star by star
with planetary systems (a star every 0.005 pc, a few hundred AU, with Hill
spheres of 27 AU or less). If the owner wants the real centre it should be
a special-case object (a nuclear cluster phenomenon with a statistical star
population and Sgr A* as a black hole), not a sector fill. Ring 0 slot 0
already holds the black-hole nucleus (`add_galactic_nucleus`).

### 4.3 Time and storage for the core fill

Speeds are the PERF.3 defaults (`DEFAULT_SECONDS_PER_SYSTEM` = 0.2 s,
`DEFAULT_BYTES_PER_SYSTEM` = 48 KiB), not measured on a server [C,
`exp3.py`]:

| Core radius (layer 0) | Rings | Sectors | Expected systems | 1 worker | 8 workers | Storage (+10%) |
|---|---|---|---|---|---|---|
| 20 pc | 0..4 | 75 | 73 thousand | 0.2 days | 0.03 days | 4 GB |
| 50 pc | 0..12 | 531 | 0.51 million | 1.2 days | 3 hours | 28 GB |
| 100 pc | 0..24 | 1,965 | 1.89 million | 4.4 days | 13 hours | 102 GB |
| **200 pc (GEN.24 suggestion)** | 0..49 | **7,857** | **7.45 million** | **17.3 days** | **2.2 days** | **403 GB** |
| 400 pc | 0..99 | 31,425 | 28.5 million | 66 days | 8.2 days | 1.5 TB |

The time is per-system generation, so a faster fill order cannot help: GEN.24
is limited by density, not traversal. PERF.3 refuses a run that would take
more than `MAX_DISK_SHARE` = a quarter of the database disk or leave under
`MIN_FREE_BYTES` = 5 GB, so the 200 pc default is refused on any disk below
about 1.6 TB. Reducing the cost needs a smaller radius, a lower per-sector
density in the core (a cap), or an accepted cost.

### 4.4 What the bulge scale radius controls

`bulge_scale_radius_pc` is x0 of the Dwek G2 profile, the half-axis along the
bar; y0 and z0 are fixed fractions of it (0.39 and 0.27). The bulge's
contribution at the centre is `bulge_amplitude x k_norm` whatever x0 is.
From the real model [C]: centre density 154.2 for any x0 at amplitude 3.11;
density at 200 pc from the centre is 118.9 (x0 = 200), 160.6 (x0 = 800) and
163.2 (x0 = 1,580); with amplitude 1 the centre is 77. A toy galaxy with
scale length 500 and bulge radius 300 still has 154 at the centre. So the
systems per core sector come from the bulge amplitude and the disc, not from
the scale radius; only the extent of the dense region follows x0; and a
fixed 200 pc core is deep inside the real bulge at the Milky Way default but
covers most of a small galaxy's bulge.

### 4.5 Recommended default and thresholds

- Core radius default: **50 pc** (531 sectors, rings 0..12) at the Milky Way
  shape, scaled for other shapes as min(50 pc, 0.25 x
  `bulge_scale_radius_pc`), so a toy galaxy with a 100 pc bulge gets 25 pc.
  200 pc stays a typed value ("the nuclear stellar disc") with the estimate
  in front of the confirm.
- Warn at 2,000 sectors (the existing threshold, reached near 100 pc) and
  name the estimate; the PERF.3 refusal then does what it already does.
- Order: ring serpentine from ring 0, layer 0 only (Boss's answer). Ring 0
  slot 0 holds the nucleus and is filled first, so `--limit N` leaves a solid
  disc.
- Backfill (GEN.23, PERF.18): see 3.5; around a 200 pc core compute the
  margin once.
- Owner's call: cap systems per sector in the core, or mark the very centre
  as a special cluster object, since the model's 1,000 systems per sector is
  already 5 to 15 times below the real nuclear disc and the cluster sector is
  off the scale.

## 5. Statements in older documents that this research supersedes

- **Bulge scale radius 200 pc** (TODO GEN.24; `orbital-updates.md` section
  11): stale; shipped default 1,580 pc (section 4.1).
- **`Computational Astrodynamics.md`, "Synthesis and Structural
  Recommendations" (spatial partitioning paragraph and its table row):** 4 pc
  uniform cubic cells indexed by a 64-bit Morton code with a 27-cell stencil,
  "replacing cylindrical geometries". Not what is built: the shipped grid is
  the cylindrical ring, layer, slot grid, `galaxy-drilldown-navigation.md`
  turned Morton keys down, and `orbital-updates.md` section 11 keeps the
  current sectors ("this sector and every sector touching it"). The 4 pc edge
  matches. Section 2 shows why a cubic index cannot address these cells.
- **`Orbital Position and Vector Update Algorithms.md`, "Database Storage
  Methods", worst-case setup:** a 10 ly core sector with 1 black hole, 50
  stars and 500 planets (551 rows). The model puts 880 to 1,060 systems in
  each 4 pc core sector (section 4.1), so that table's 70 KB point-mass table
  and 32 MB octree figures are lower bounds on a generated core sector, and
  the real centre is denser still (section 4.2). The black-hole part (an
  analytic point mass for Sgr A*) stands.
- **GEN.101 text:** "D = ceil(log2(2R))" is one level short for R a power
  of two, and "face-adjacent unit steps ... zero backtracking" cannot be met
  (section 1).

## Evidence notes

Everything marked [C] is reproducible from the scripts. Everything outside
the repository is [R] and unverified; the research environment could read
only search-result text, not papers, and web search was unavailable for this
topic. Verify when paper access is allowed:

- All astronomy figures in 4.2: nuclear star cluster mass (2.0 to 2.5e7
  Msun) and half-light radius (about 4 pc), Schodel et al. 2014 and
  Feldmeier-Oelkers et al. 2014 (medium confidence); nuclear stellar disc
  mass (about 1e9 Msun), radius and scale height, Launhardt et al. 2002 and
  Sormani et al. 2022 (medium); CMZ gas mass and radius (medium); Sgr A* mass
  4.3e6 Msun (high); bulge mass 1.5e10 Msun (matches the repository's BHG16
  citation). The inner-parsec densities are derived from an assumed slope
  and are low confidence.
- Stars per solar mass in the dense core (about 2) and the 0.1 systems per
  cubic parsec solar-neighbourhood figure are common but unchecked. The
  enclosed mass inside 200 pc (about 1e9 Msun) is recalled.
- Papers: Skilling 2004 (AIP Conf. Proc. 707); Butz 1971 (IEEE Trans.
  Computers C-20); Hamilton and Rau-Chaplin 2008 (Inf. Process. Lett. 105);
  Moon et al. 2001 (IEEE TKDE 13); Niedermeier, Reinhardt, Sanders 2002;
  Haverkort and van Walderveen 2010; Haverkort 2017 (J. Comput. Geom.); Xu,
  Nguyen, Tirthapura, ICDE 2018 (onion curve); Faloutsos 1986; Itai,
  Papadimitriou, Szwarcfiter 1982 (SIAM J. Comput.); Papadimitriou and
  Vazirani 1984; Umans and Lenhart 1997; Hjaltason and Samet. These are
  believed to exist; titles and years need checking. That 3D grid-graph
  Hamiltonian path is NP-complete is believed, not confirmed. It does not
  affect any conclusion: the parity proof is self-contained.
- The 1.5 to 2 times slowdown of Python 3.9 against the 3.13 timed here is a
  rule of thumb.
- The PERF.3 defaults (0.2 s and 48 KiB per system) are constants in
  `generation/stats.py` and have not been measured on a real server; the
  core time and storage rows scale with them.

## Sources

- Repository files read: `docs/TODO.md` (GEN.101, ADM.29, ADM.30, GEN.24,
  GEN.97, ADM.31, PERF.18, NAV.43), `galaxy-coordinate-system.md`,
  `galaxy-disk-density.md`, `reproducible-galaxies.md`,
  `galaxy-drilldown-navigation.md`, `orbital-updates.md`,
  `Computational Astrodynamics.md`, `Orbital Position and Vector Update
  Algorithms.md`; `src/planetgen/galaxy/geometry.py`, `density.py`,
  `skeleton.py`, `viewport.py`, `src/planetgen/generation/run_galaxy.py`,
  `run_common.py`, `stats.py`, `src/planetgen/physics/orbits.py`,
  `src/planetgen/tuning.py`, `src/planetgen/db/models.py` (sector indexes).
- PyPI JSON API: https://pypi.org/pypi/hilbertcurve/json ,
  https://pypi.org/pypi/numpy-hilbert-curve/json ,
  https://pypi.org/pypi/pyzorder/json
- Raw source (module docstring only):
  https://raw.githubusercontent.com/galtay/hilbertcurve/main/hilbertcurve/hilbertcurve.py
- Experiment scripts (research scratchpad): `hil.py` (Hilbert and Morton,
  metrics), `t1.py` (Hilbert verification), `exp1.py`, `exp1b.py`,
  `exp1c.py`, `exp1d.py` (ball experiments, parity), `exp2a.py` to `exp2g.py`
  (real grid: counts, neighbours, orders, serpentine, discs), `exp3.py` (core
  density and cost), `exp4.py` (backfill redundancy), `exp5.py` (estimation
  cost, sampling); raw results in `exp1.json`. The scripts are not in the
  repository; a developer implementing GEN.101 should port the greedy walk
  and serpentine into `src/tests/` with the 50 pc sphere and 30-ring cylinder as
  fixtures.
