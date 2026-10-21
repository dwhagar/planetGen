# Spherical-shell sectors with a round core: what it would mean

Research only, 2026-10-15. No code, TODO.md or schema changes. Boss decides.

**Update, 2026-10-15 06:39Z:** Boss decided to go ahead with spherical shells, with pole caps and wedges starting at a 4 pc gap. The plan and the final rules (which differ from S1 below in the caps and the multiples-of-6 wedge rule) are in [spherical-sector-coordinates.md](spherical-sector-coordinates.md). The recommendation in this note was to keep the cylinder; that is superseded by Boss's decision.

Question (Boss, 2026-10-15 05:51Z): what would happen if sectors were spherical shells cut into arc segments, with the core a spherical sector 4 pc across, so that at sector level it is still basically a cube for most of the galaxy?

Informs: MAP.151 (the 3-ary block ladder, MAP.151b, is the next piece), MAP.152, GEN.101, GEN.24, ADM.30.

Evidence tags: [C] computed in this research and re-runnable (scripts and raw output under `/mnt/project-files/research/scripts/spherical-sectors/`, `results.txt` beside them), [S] read in the repository. Nothing here is [R] recalled; every number below is from the scripts or the code. The density is the real default model (`planetgen plan` shape, 4 pc edge, 6.31 systems per cell at density 1, the skeleton's own qualification test).

## 1. Answer in brief

**It works, it is not worth doing, and nothing in the project asks for it yet.** My view, plainly: keep the cylindrical grid.

1. **The grid is buildable and tiles exactly.** Shells 4 pc thick, an odd number of latitude bands per shell (the plane stays in the middle of a band), slots per band proportional to the band's circumference. A point maps to its cell in closed form (20,000 random points, 0 misses) and the cell volumes sum to the ball's volume [C]. It does not bring back the problems of the archived design (section 3).
2. **For the disk it changes almost nothing.** Candidate cells: 23,917,889,703 against 23,920,622,421 on the cylinder (-0.011%), both 3,764 rings or shells [C]. Cell volume is within 0.999 to 1.001 of an edge cubed for 99.8% of cells (the cylinder's master-wedge rule gives 0.945 to 1.062). 99.99999% of cells have all three edges within 10% of each other. "Basically a cube for most of the galaxy" is true, and was already true.
3. **The core does not get lighter.** The 4 pc ball holds about 509 systems, but its neighbours are denser: the densest cell is 1,187 systems against 1,114 on the cylinder (section 4.3). The model has no central cusp, so a smaller centre cell does not remove the heavy cells next to it. The core-readability lever stays presentation, as the 2026-10-08 study concluded.
4. **The cost is the cylinder's best property.** Aligned columns and flat layers go away. A block or region query that today needs one range per ring needs one per (shell, band), 27 times more for a 27-wide block; face neighbours double (6 to 7 become 10 to 15); the block ladder (243/81/27/9/3/1), `drill.py`, `region.py`, `window.py` and the three map scripts are rebuilt. 46 source files (882 lines naming ring, layer or slot) and 79 of 373 test files change; the galaxy and every sector, object ID and seed-derived sector content are regenerated (section 5).
5. **If the wish is a round, special centre**, the cheaper way is the one already proposed in `fill-order-curves-and-core.md` section 4.2: the nucleus as a special object in ring 0 slot 0, not a different grid (section 7).

If Boss wants it anyway, decide **before** MAP.151b and the rest of the Maps ladder merge: every ladder merge on the cylinder raises the price (section 7).

## 2. What is proposed, made precise

Edge `e` = 4 pc. This is my reading of Boss's sentence; the two variants are mine.

- **Core (shell 0):** a ball of radius `e/2`, diameter 4 pc, volume 0.524 e^3 (33.5 pc^3).
- **Shell `k` >= 1:** `r` in `[(k - 1/2) e, (k + 1/2) e)`, centreline radius `k e`. Same thickness as a ring; the core ball is half a shell, so shell `k` sits where ring `k - 1/2` does today.
- **Latitude bands:** `n(k)` = the odd integer nearest `pi k`, equal steps of colatitude. Odd puts the galactic plane in the middle of a band, as layer 0 straddles it today. A band's arc at the centreline is within 0.1% of `e` at large `k`.
- **Slots (S1, the variant studied):** band `b` of shell `k` has `N = max(1, round(2 pi k sin(colatitude of its centre)))` equal slots, so every cell is close to a cube at every latitude, up to the axis. Near the axis the first bands hold 3 slots, a polar cap of three pie wedges like ring 0 today.
- **S2 (aligned meridians):** one slot count per shell, `round(2 pi k)`, the same in every band, which keeps every slot boundary a vertical half-plane through the axis like the cylinder. Studied as the "keep the columns" alternative and rejected (section 4.2).

Address: `(shell, band, slot)`, unique, closed form both ways. Local axes, positions inside a cell and the stored centre work as today; sampling inside a cell is separable (cube root for `r`, uniform in `cos` of colatitude, uniform longitude).

## 3. What changed since the archived design

The archived "Track C" design (`archive/galaxy-coordinate-system-shells.md`, built 5.3.6, replaced in 6.0.0) put **cube** sectors on spherical shells at Fibonacci-sphere points. Its problems and what happens to each:

| Problem then | Now |
|---|---|
| Cubes cannot tile a sphere, so gaps and overlaps, then a Voronoi prism per sector (`sector_vertices`) and a Fibonacci-number neighbour search | Gone. Arc segments are the cells; they tile by construction and every lookup is closed form [C]. |
| Most shell slots sit far off the disk | Still true of the address space, but the stored outline (one band range per shell, about 3,764 rows, like `galaxy_layer`) bounds it, and the outline holds the same 23.9 billion cells as the cylinder [C]. |
| A per-shell band skeleton (`galaxy_shell_band`) | Replaced by the same bounds-before-generating rule; the table would be a shell-keyed twin of `galaxy_layer` and `galaxy_column`. |

What the archive could not know, because it came before them: the 243/81/27/9/3/1 ladder and its nested wedge rule, `drill.py` and the stage cache, region and window frames (MAP.151a/c), the 80-bit object IDs (the first layout; 96 bits is the default since Boss, 2026-10-15) with a 40-bit birth-sector field, per-sector seeds derived from an address string, bright-star and phenomenon scatter tables keyed by ring, layer and slot, the Redis scatter queue, and the fill-order work. All of them are written against ring, layer and slot, and the archive's reasons for leaving shells (a cube cannot tile a sphere, off-disk addresses) do not apply to arc segments, so **the old reasons no longer settle the question; the cost of the new things does.**

## 4. What I measured

All for the default Milky Way shape, outline built with the skeleton's own test (bound density at the cell centre at least 1/6.31).

### 4.1 Counts and cube-ness

| | Cylinder (today) | S1 (varying slots per band) | S2 (one slot count per shell) |
|---|---|---|---|
| Rings or shells | 3,764 | 3,764 | 3,764 |
| Candidate cells | 23,920,622,421 | 23,917,889,703 | 26,138,221,065 |
| Volume per cell, share of cells within | 0.945 to 1.062 (all) | 0.999 to 1.001 (99.8%), extremes 0.85 to 1.15 near the axis and centre | median 1.00; 6.7% of cells below 0.5, 1% below 0.08 |
| Cells with the three edges within 10% | 100% (arc 0.94 to 1.06, others exactly 1) | 99.99999% | 82% |
| Cells more than 25 degrees off the plane | n/a | 11% (4.3% above 45) | 18% (10% above 45) |

S1 is the only variant that meets "a cube almost everywhere". S2 pinches cells in the bulge and over the axis, where a fifth of its cells sit above 25 degrees.

### 4.2 Why S2 does not rescue the columns

S2 keeps vertical slot planes, but a longitude slot at latitude `L` is `cos L` as wide, so cells get narrower with height. In the default outline 7% of cells would be under half an edge cubed, and the cell count grows 9%. The cylinder is the only grid here that is a cube at every height, because a cylinder is the same at every height. That is the property being traded away.

### 4.3 The core [C]

| | Cylinder | S1 |
|---|---|---|
| Centre cell | ring 0 layer 0 slot 0, a pie wedge, volume 1.05 e^3, about 920 systems at its centre density | ball, volume 0.52 e^3, about 509 systems |
| Densest cell | 1,114 systems (rings 1 and 2, any layer 0 to 3) | 1,187 systems (shells 3, 5 and 7, polar-cap cells of volume 1.15 e^3) |
| Shell 1 | n/a | 12 cells of 1.13 e^3 whose outer face is 9 times the inner |
| Cells whose outer face is within 10% of the inner | rings 10 and out (40 pc) | shells 21 and out (84 pc), 41,453 cells inside |

Density at the centre is 154 times the local density and almost flat across the core. A round centre cell makes that one cell lighter and leaves 12 neighbours at about 1.1 e^3 and the same density. Only 0.0002% of cells are inside shell 21, so Boss's "most of the galaxy" is safe; the core is the one place where the shape is not a cube, on both grids.

The nucleus itself is not a sector-sized problem on either grid: `fill-order-curves-and-core.md` section 4.2 puts the real centre at 10^6 to 10^7 stars in a 4 pc cell, which is a statistical object, not a fill.

### 4.4 Neighbours, addressing, outline [C]

- **Face neighbours:** cylinder 6 to 7 (mean 6.2); S1 10 to 15 (mean 12.5 near the plane, 13.6 at 60 degrees); S2 10 to 14 (mean 12). A neighbouring shell's bands and slots are offset in two directions instead of one, so a cell touches up to four cells above and four below instead of one or two. Anything that walks neighbours (the greedy fill, the nav graph, `neighbor_addresses`) does about twice the work and the "brick course" offset becomes a two-way offset.
- **Closed form:** shell from `|p|`, band from `acos(z/|p|)`, slot from `atan2`. No search anywhere.
- **Outline and ID bits:** the outline reaches the axis up to shell 927 (about 3,700 pc above the plane); the largest band offset from the plane is 1,456 and the largest slot count is 23,650. The 40-bit birth-sector field of the 80-bit object ID (the first layout; 96 bits is the default since Boss, 2026-10-15; 12 + 12 + 16 bits, `object-id-options.md`) still fits: shell 3,764 < 4,096, band offset 1,456 < 2,048, slot 23,650 < 65,536. The field widths survive; every value in it changes.

## 5. What it would change in the project

**Gone with the cylinder:**

- **Layers.** There is no global height index. A band is a cone: band `j` of shell `k` sits at `z = k e sin(j pi / n)`, about `j e - j^3 e^3 / (6 r^2)`, so it is flat to a fraction of a parsec between neighbouring shells but drifts down across the galaxy (about 3 pc low at 1,000 pc above the plane at 8 kpc, 18 pc at 3 kpc), and no slab is flat. Everything that reads "layer" means a flat slab today: the slab slider, `galaxy_layer` and `galaxy_column`, ADM.30's decision "cylinder height = X layers", GEN.24's "layer 0 only" (a plane band instead), the Zone and slab wording on the Galaxy pages.
- **Aligned columns and one arc per ring.** `region.py` lists a round region as one slot arc per ring, the same for every layer; `window.py`, `span.py` and the stage queries follow. With bands, each (shell, band) has its own arc: a 27-wide block is 27 SQL range clauses today and 729 on S1 (243 wide: 243 against 59,049) [C]. The stage queries already pay this pattern (`drilldown-region-sizes.md` section 5.3 counts 243 clauses at m = 243); it would need a different key (a sortable curve key per cell) to get back to a few ranges, and `sectors`, `bright_stars`, `sector_stats`, `phenomenon_scatter`, `galaxy_layer`, `galaxy_column` and `galaxy_block_stats` are all keyed by `(ring, layer, slot)` [S].
- **The block ladder and master wedges.** 3-ary blocks cut on shared lines need a latitude master rule (bands a multiple of a master count) on top of the longitude one. A crude try with odd 3^z masters gives band arcs of 0.89 to 1.11 of an edge, against 0.94 to 1.06 for wedges today [C]. Longitude masters can continue to within 4 to 7 degrees of the axis at shells 100, 1,000 and 3,000, past which blocks over the axis need a special case [C]. A 27-block is also a frustum in both angles instead of one, so blocks, tiles, stage keys, `drill.py`, `galaxyprisms.js` (1,047 lines), `galaxystageview.js` (2,561) and `galaxymap3d.js` (3,306) are rebuilt rather than adapted. This is design work, not a rename. I did not design the S1 ladder; the numbers above only say it is possible.

**Rebuilt or renamed:**

- 46 non-test source files name ring, layer or slot (882 lines); the heavy ones are `db/query.py` (124 mentions), `db/store.py` (107), `galaxy/geometry.py` (112), `generation/run_plan.py` (63), `generation/bright_stars.py` (62), `galaxy/skeleton.py` (50), `generation/phenomenon_scatter.py` (40), `generation/run_galaxy.py` (39) [S]. 79 of the 373 test files use them (707 lines), including `test_galaxy_gen.py`, `test_galaxy_geometry.py`, the fuzz and boundary tests and the map tests.
- That is the Database (migration, seven tables, block stats), Maps (ladder, stage view, prisms), Generation (fill order, scatter, bright stars, region fills) and UX (Galaxy and Sector pages, slab wording) lanes at once, plus the Test Suite thread.
- The in-flight Maps ladder (MAP.151b with its data-only migration, then MAP.152) is built for the cylinder and would be set aside, and MAP.151a/c (`region.py`, `window.py`, merged) redone. The 8.1 fly-through headline rests on that ladder.

**Reset:** like v32, v33 and v35, the migration deletes every galaxy-placed sector, system and phenomenon and rebuilds the outline. Because a sector's seed is the SHA-256 of the galaxy seed and its address string, **the same galaxy seed fills every sector differently**, and because object IDs carry the birth sector, **every object ID and sector ID changes**. The no-backward-compatibility rule makes this acceptable; it is the same kind of reset Boss accepted on 2026-09-24, on a codebase many times larger.

**Unchanged:** the density model and the outline test, E = 6.31 systems per cell, the cell count and so the generation and storage cost (+/-0.01%), seeds and reproducibility rules, orbits and the galactic potential (they are continuous), Hill-sphere and radius handling (`galactic_radius_pc` is stored already).

## 6. What it would give

Real, but small:

1. **A centred nucleus.** The quasar sits mid-cell in a ball instead of on the vertex shared by three wedges at ring 0, layer 0, slot 0, and "the core out to a radius" becomes whole shells (GEN.24's 50 pc is shells 0 to 12, a ball, not a one-layer disc).
2. **An address that is a navigation coordinate.** Shell is distance from the core; band and slot are mark and bearing in the galactic frame, the same course the NAV page and the object-ID bearing and mark use (`navigation-frames.md`). A sector designation would read as a quantised course.
3. **Even cells.** Expected systems per cell vary by +/-0.1% instead of +/-6% from rounding (back to about +/-11% in latitude if the ladder's master rule is added).
4. **Natural for spherical things.** Globular-cluster tidal spheres, a halo, a spheroidal galaxy (`multiple-galaxies.md` has `axis_ratio_thin` 1.0 for ellipticals). The cylinder holds these too, as cubes at every height, so this is tidiness, not capability.

Not given: fewer cells, lighter cells, a cheaper fill, a faster map, a better disk. The disk is where almost all the volume is, and there the two grids are the same cells to a fraction of a percent.

## 7. Options

| # | Option | Verdict |
|---|---|---|
| A | **Keep the cylinder** and, if a round centre is the point, make the nucleus a special object in ring 0 slot 0 (a statistical nuclear cluster with the black hole, `fill-order-curves-and-core.md` section 4.2). | **Recommended.** Zero grid change. The cluster's half-light radius is about 4 pc [R, in that note], twice the ball's radius, so neither grid makes it one cell. |
| B | **S1 spherical shells** as a full reset. | Buildable, a different galaxy afterwards, a rebuild of the Maps, Database and Generation work done and in flight. Worth it only if Boss wants spherical symmetry as a feature (nav-native addresses, halo and cluster work, spheroidal galaxies) and accepts the cost. If yes, **decide before MAP.151b** and the rest of the ladder merge: it is the point where the cost is smallest. |
| C | S2, one slot count per shell. | No: it keeps the columns but loses the cube (section 4.2). |
| D | A ball only at the core, cylinder elsewhere. | Not a grid: removing a ball from the three ring-0 wedges leaves cells that are neither wedges nor shells, and every layer above and below still needs ring 0. Only as option A's special object. |

The core-readability question was asked and closed on 2026-10-08 ("leave it all be", 4 pc stays). Nothing measured here is new information against it: the heaviest cells are as heavy on shells as on rings.

## 8. What I could not measure

- **Draw time and database timings.** No browser or MySQL here. The 27 times more range clauses is a count of clauses, not a query time; cell prisms with curved faces were not drawn.
- **The S1 block ladder.** Only the arc ranges of a latitude master rule and the longitude run-out near the axis were computed; the nested ladder, stage keys and the wire format are not designed.
- **Neighbour counts** come from random samples (400 per row, 2,000 for the cylinder), not all cells.
- **Arm contrast** is not in the counts; the qualification test uses the same upper bound as the skeleton, so both grids see the same outline.
- **Effort.** I give the size of the surface (files and lines), not a time; the cost of the reset on a server (regeneration time) is the existing per-system rate times the sectors Boss chooses to fill and is the same on both grids.
- **Shape of the core on S1** assumes the 3-slot polar caps and the rounding rule above; a different cap rule moves the shell 1 to 5 figures, not the rest.

## 9. Method and files

`sphgrid.py` (the grid, outline and statistics) and four runners, in `/mnt/project-files/research/scripts/spherical-sectors/`: `01_counts_and_cubeness.py` (sections 4.1, 4.2), `02_core_cells.py` (4.3), `03_neighbours_and_tiling.py` (4.4, the SQL clause counts and the 20,000-point round trip), `04_outline_bits_masters.py` (outline bits and master rules); `results.txt` is their output. Run from the repository's `src` directory with `python3 -I` and numpy and scipy installed; they import `planetgen.galaxy.density`, `skeleton` and `geometry` and use the real default shape and outline test. Code surface counts are `grep` of `ring_index`, `layer_index`, `ring_slot_index`, `ring_sector_count`, `ring_master_count`, `sector_address_at`, `neighbor_addresses` and friends over `src`, 2026-10-15, main at 8.0.1248.
