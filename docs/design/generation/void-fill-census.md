# Filling the voids on the spherical grid, again (GEN.102)

Boss (2026-10-07 11:47Z): "Investigate filling all void space (space that plots to have a 0 or nearly 0 star density) at
once, since those sectors have the least in them, we can afford to fill them." The first answer
([sampling-backfill-and-resume.md](sampling-backfill-and-resume.md) section 5, 2026-10-09) counted the old ring grid and
recommended against a pre-fill. The grid has since become spherical shells, bands and wedges (GEN.235), the reach is 262 kpc
with a halo (GEN.234, GEN.301), and a cell over the star ceiling is only partly filled (GEN.281). This note does the count
again on today's grid and says what changes.

Informs: GEN.102, GEN.276, GEN.301.
Status: research, 2026-10-16. Evidence tags: [C] computed here and re-runnable, [S] read from the repository, [R] recalled and
unconfirmed.

## Short answer

The recommendation stands: do not pre-fill voids; keep the fill on visit (`ensure_sector_generated`). The count barely
moved, and a new fact makes the question larger than it was. Inside the outline a "void" costs the same fixed time per
cell as any other, so it is the worst place to spend a pre-fill. Beyond the outline there are 1.2 x 10^15 cells and no
cell holds even one expected system, so filling "all the space" cell by cell cannot be done at any price and the halo has to come from the
region and scatter paths (GEN.301, GEN.276). Nothing below is a fixed number: it is the census of the galaxy the
generator makes from today's shape defaults, and it is re-run from the shape.

## 1. How it was counted [C]

The shape is the generator's present default (`planetgen generate` shape options: disk scale length 2600 pc, height 300
pc, bulge 1580 pc and amplitude 3.11, two arms, 15 degree pitch, arm and inter-arm density 1.4 and 0.6, nuclear amplitude 1)
with no components, since the globular clusters and satellites are drawn from the galaxy seed and add systems only in
their own cubes. Those defaults are Milky Way scale numbers; when the shape is derived from the seed
(the no-fixed-content rule) this census is re-run on the derived shape and the table below changes with it.

- The outline is `skeleton.build_shell_extents`: 3,763 shells (the last is 3,762, 15 kpc out), 2.418 x 10^10 candidate cells.
- `E` (systems per cell at density 1) is 6.306. A cell's expected systems are `E * relative_density * (cell volume / edge^3)`.
- 3,000 cells were drawn uniformly from the outline (`GalaxyBounds.random_address`, seed 102), and the density was read at each
  cell's centre. Sampling error on the shares is about one percentage point.
- Beyond the outline: 24 radial shells from 15 kpc to the 262 kpc reach, 400 random directions each.

Re-run with `/mnt/project-files/galaxy-studies/gen102-void-census/voids.py` and `halo.py` from a checkout (a few seconds).

## 2. Inside the outline

Density quantiles (1.0 is the calibration point): minimum 0.072, 1% 0.085, 5% 0.119, median 0.42, 95% 7.0, 99% 35. As
before, nothing inside the outline is near zero: a cell is kept only if its best longitude reaches 1/E = 0.159.

| "Void" definition | Share of cells | Cells | Systems per cell | Share of all systems |
|---|---|---|---|---|
| density below 0.1 | 2.4% | 5.9 x 10^8 | 0.55 | 0.10% |
| expected systems below 1 | 13.3% | 3.2 x 10^9 | 0.78 | 0.80% |
| expected systems below 2 | 40.3% | 9.8 x 10^9 | 1.21 | 3.8% |
| the whole outline | 100% | 2.4 x 10^10 | 12.9 | 100% (3.1 x 10^11 systems) |

These agree with the first count (12.4% of cells below one expected system, 0.77% of systems), so the change of grid did
not change the answer. The sparse cells sit at a median radius of 12 to 14 kpc and a height of 0.4 to 0.6 kpc: the outer disk
between the arms. About 14% of all cells draw no system at all (the mean of exp(-expected)).

## 3. Time and storage

The cost model is the timing report's (PR #657): about 40 ms fixed per cell and about 16.4 ms per system on one worker, 1.87
times faster on four [S, `/mnt/project-files/bright-star-timing/report.md`]. Storage is about 40 KB per system
(not re-measured). The Generate page's own estimate (`generation/stats.estimate`) uses the server's recorded speeds and,
before any run, 0.2 s and 48 KB per system, so a fresh server would quote the last two columns' larger figures.

| Fill | Cells | Systems | One worker | Four workers | Storage | At the page's default speed |
|---|---|---|---|---|---|---|
| density below 0.1 | 5.9 x 10^8 | 3.2 x 10^8 | 333 days | 0.5 years | 13 TB | about 2 years |
| expected systems below 1 | 3.2 x 10^9 | 2.5 x 10^9 | 5.4 years | 2.9 years | 100 TB | about 16 years |
| expected systems below 2 | 9.8 x 10^9 | 1.2 x 10^10 | 18.5 years | 9.9 years | 470 TB | about 75 years |
| the whole outline | 2.4 x 10^10 | 3.1 x 10^11 | 193 years | 103 years | 12.5 PB | not computed |

Time is a fixed cost per cell, so a sparse cell costs 70 to 90 ms per system against 16 for an ordinary one. Boss's premise,
that the sparsest cells are the cheapest to fill, does not hold on the cost model.

## 4. What is new: the space beyond the outline

The halo is in the density now. Its relative density is 4 x 10^-4 at 15 kpc, 5 x 10^-5 at 30 kpc, 2 x 10^-7 at 100 kpc and
3 x 10^-9 at 262 kpc, all under the outline's threshold of 0.159, so none of it is a candidate cell [C].

| Beyond 15 kpc, to the 262 kpc reach | |
|---|---|
| Cells | 1.18 x 10^15 |
| Systems in all of it | 5.2 x 10^9 (1.6 x 10^9 of them the halo, the rest the disk's outer edge) |
| Systems per cell | 4.4 x 10^-6 on average, and under 1 in every one (a cell with one would be inside the outline) |
| A cell-by-cell fill at 40 ms a cell | about 1.5 million years on one worker |

So there are two different questions. The sparse rim inside the outline is a question of how much cost to pay
(section 3). The space outside it is not a fill at all: nearly every cell is empty, and the 5 billion systems it holds have to be placed
by drawing the count for a region and then the positions, which is what the block-first draw and the region query do
(sampling-backfill-and-resume.md sections 2 to 4; GEN.276, GEN.301). That path never visits an empty cell. A pre-fill of
the outer space is not on the table, and the outline keeps it out of the candidate cells, so nothing needs to change.

## 5. Recommendation

1. **Do not pre-fill the sparse rim.** Keep the fill on visit. It costs nothing until someone looks, and a visited sparse cell is
   saved and marked generated even when it draws no system (GEN.76, `run_galaxy.ensure_sector_generated`).
2. **The old open question is answered.** "Unfilled" and "filled and empty" can already be told apart: a filled cell has a
   `sectors` row (unique on shell, band and wedge) whether or not it holds a system, so "empty needs no row" is not needed
   to save writes and would lose that.
3. **Define a void as expected systems below 1**, which matches the existing rule that a cell needs one expected system to
   qualify, and not as "near zero". It is 13% of the outline's cells and 0.8% of its systems.
4. **If the rim is wanted anyway,** the only lever is the fixed cost per cell (the timing report's three levers could roughly
   halve the 40 ms), and the first target is density below 0.1 (333 days on one worker, 13 TB), not the whole rim.
5. **The halo comes from regions,** not cells (section 4). That is GEN.301's work and needs no void rule.

## 6. What this does not settle

- The census uses today's shape defaults. When the shape is derived from the seed, the shares move; the shape of the
  answer (a fixed cost per cell, a thin sparse rim, a huge empty outside) does not.
- The per-system cost and the 40 KB are the timing report's, not re-measured since the cell star ceiling and the orbit
  passes landed. A fresh timing run on the current build would tighten the third table. This is a test candidate.
- Boss's question (default no): do you still want the sparse rim filled?
