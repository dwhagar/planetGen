# Spherical sector coordinates: the plan

Status: plan, 2026-10-15. Boss decided to move the galaxy's sector grid from the cylindrical ring/layer/slot grid to spherical shells (Boss, 2026-10-15 06:39Z). Nothing is built. This file is the brief for the dedicated thread "Spherical Sector Coordinates" and the source for the TODO items in section 6. The study behind it is [spherical-shell-sectors.md](spherical-shell-sectors.md); its scripts and raw output are in `/mnt/project-files/research/scripts/spherical-sectors/` (`05_caps_and_wedge_schemes.py` and `results-caps.txt` are the ones for the rules below).

Evidence tags: [C] computed with those scripts on the real default density model, [B] Boss's words, [S] read in the repository.

## 1. What Boss decided

- [B] Use math terms: **shell** (radius), **band** (horizontal), **wedge** (vertical: two meridian lines from pole to pole). Shell, band and wedge together name one 3D segment of the sphere.
- [B] The core is a ball 4 pc across, shell 0 (the old "layer 0 is the sphere").
- [B] The depth of a shell is always 4 pc. Everything else may be adjusted.
- [B] Keep bands and wedges as close to 4 pc arcs as possible, especially at the poles.
- [B] Each shell has a circular **cap** at each pole: a circular segment, nominally 4 pc radius and 4 pc deep. Wedges start at a 4 pc gap on the cap's rim, so if the cap must be a little other than 4 pc to make that work, balance it.
- [B] This goes to TODO at top priority. All other work stops when each lane's current task finishes (no chaining), and a new thread starts for this.

## 2. The grid, made precise

Edge `e` = 4 pc. All of this is closed form in both directions.

- **Core (shell 0):** ball of radius `e/2` = 2 pc.
- **Shell `k` >= 1:** `r` in `[(k - 1/2) e, (k + 1/2) e)`, centreline radius `R = k e`.
- **Cap:** the part of the shell within colatitude `alpha` of a pole, one cell at each pole. Its rim is cut by `W1 = 6` wedges whose rim gap is exactly `e`, so the rim circumference is `6 e` and the cap's chord radius is `rho = 6 e / (2 pi) = 3.82 pc` (nominal 4 pc is 6.28 wedges, so 6 is the nearest whole number). `alpha = asin(rho / R)`. Volume about 2.9 edge cubes [C].
- **Bands:** between the caps, `n` bands of equal colatitude step, `n` the odd integer nearest `R (pi - 2 alpha) / e`, so the galactic plane runs through the middle of the equatorial band (band 0, signed bands above and below).
- **Wedges in a band (Boss, 2026-10-15 07:17Z: one rule everywhere, no special regions):** every band of every shell, in the galaxy or in any other shape of sphere, has `W(b) = 6 * m` wedges, where `m` is a number whose prime factors are only 2, 3, 5 and 7 and which is nearest (in ratio) to `2 pi R sin(phi_b) / (6 e)`, `phi_b` being the band's centre colatitude; `m` is at least 1. The six pole meridians therefore bound every band, and the rule never looks at the outline, the window or the neighbouring shell, so `address_at` stays closed form. Caps keep 6 wedges at the rim. Measured on shell 3,000 [C, `10_uniform_smooth_rule.py`]: arcs 0.981 to 1.019 of 4 pc, 95% of neighbouring band steps nested, 444 count changes per shell; 240 distinct equator counts over shells 4 to 3,764 (about 16 shells per count). Near the poles, small counts let arcs stray further (about 0.76 to 1.06 at shell 10). This replaces the earlier "fine in the window, nested above it" rule (06:57Z), which depended on the outline. Wedge 0 starts at the +X meridian.
- **No undefined space (Boss, 2026-10-15 06:59Z):** the core ball, every shell, every cap and every band and wedge together tile all of space: each point has exactly one address, with no gap and no overlap, whether or not that sector has been generated. The outline (item 2) only says which cells may be generated; it never decides which cells exist, and `address_at` answers for any point, inside or outside the outline, as far as the 40-bit birth-sector field reaches (shell below 4,096, about 53,000 ly, beyond the galaxy's 3,764 outline shells). Shells 1 to 3 (hand-tabled), the caps, and the primes-up-to-7 wedge rule all have to meet this.
- **Address:** `(shell, band, wedge)`; a cap is band `+/-` its end index with wedge 0. Closed form: shell from `|p|`, band from `acos(z / |p|)`, wedge from `atan2(y, x)`.
- **Local frame and sampling:** unchanged in meaning (+X radial, +Y toward increasing longitude, +Z galactic north). Sampling inside a cell is separable: cube root for `r`, uniform in `cos(colatitude)`, uniform in longitude; a cap samples a disc.

Measured on the default galaxy, shells 4 and out, outline built with the skeleton's own bound test [C]:

| | Cylinder today | This rule (multiples of 6) |
|---|---|---|
| Candidate cells | 23,920,622,421 | about 23.933 billion (+0.05%) |
| Wedge and band arcs | 0.945 to 1.062 of 4 pc | 99.8% of cells within 0.5% of 4 pc |
| Cell volume (edge cubes) | 0.945 to 1.062 | 0.997 to 1.003 for 99.8%; 0.74 to 1.44 in the few cells near a pole |
| Cells with edges within 10% | 100% | 100% |
| Polar caps | none | 2 per shell up to shell 927: about 2,000 cells of 2.9 edge cubes |

Other rules tried: any wedge count per band (volume 0.94 to 1.09, but boundaries shift between bands so blocks cannot cut on shared meridians) and doubling 6, 12, 24 (only 23% of cells within 10% of a cube). Cap size by rim wedges: 5 gives radius 3.18 pc (2.0 cubes), 6 gives 3.82 (2.9), 7 gives 4.46 (3.9).

Known costs, from the study: face neighbours go from 6 to 7 up to roughly 10 to 15; flat height slabs and the aligned vertical columns are gone, so a region is an arc per (shell, band), not per ring (27 times the range clauses for a 27-wide block); a cap near the dense core holds about 2,800 systems against 1,114 for today's densest cell; the 80-bit object ID's 40-bit birth-sector field still fits (shell below 4,096, band offset up to 1,456 below 2,048, wedge up to about 23,650 below 65,536).

## 3. Decisions still open

The plan assumes the recommended answer for any row not marked Decided.

| Question | Assumed | Alternatives |
|---|---|---|
| Wedge count per band | Decided (Boss 07:17Z): 6 * m, m with primes 2, 3, 5, 7 only, nearest the 4 pc arc, same rule in every band and shell | Blocks cut on cells: decided no, blocks hold the cells whose centre they hold (Boss 06:58Z), so the ladder item picks a block count per level |
| Wedges meeting a cap | 6 (cap radius 3.82 pc) | 5 (3.18) or 7 (4.46) |
| The plane | Middle of the equatorial band (odd bands) | On a band boundary (even bands) |
| Galaxy Map ladder | Decided (Boss 07:22Z, "proceed with your plan"): built after coordinates and migration | In parallel (not adopted) |
| Object IDs in the conversion | Decided (Boss 07:32Z): renumber; every object gets a new ID from its new sector in the DB.32 conversion, positions unchanged | Keep the old IDs (not adopted) |
| Equal-volume bands | Decided (Boss 07:22Z): not adopted; the wedge rule stays primes up to 7 | Bands of equal volume (not adopted) |

Settled in the first item (GEN.224): shells 1 to 3 are too small for the cap rule to look right (shell 1's cap arc radius is 3.82 pc on a 4 pc sphere), so they are a hand-checked table in `spherical.py` (`SMALL_SHELL_TABLE`): shell 1 has caps cut at 37.5 degrees and one band of 12 wedges (cell volumes 0.9 to 1.4 edge cubes, against 4.8 for the caps and 0.67 for the cells under the rule), shells 2 and 3 keep the rule's own result. Still open: whether the migration's data delete needs a fresh word from Boss when it merges (it is destructive; Boss accepted the same for the v32 and v33 grids).

## 4. What is already built that this replaces

[S] `galaxy/geometry.py` (821 lines), `galaxy/skeleton.py`, `galaxy/region.py`, `span.py`, `window.py`, `drill.py`, `galaxy/uid.py` and the sector parts of `object_uid.py`; `static/galaxyprisms.js`, `galaxystageview.js`, `galaxymap3d.js`, `galaxyblocks.js`, `galaxytilecut.js`, `galaxystages.js`; `db/store.py`, `db/query.py`, `db/models.py`, `db/block_stats.py`; `generation/run_plan.py`, `run_galaxy.py`, `bright_stars.py`, `phenomenon_scatter.py`; the seven address-keyed tables (`sectors`, `bright_stars`, `sector_stats`, `phenomenon_scatter`, `galaxy_layer`, `galaxy_column`, `galaxy_block_stats`); 46 non-test source files naming ring, layer or slot (882 lines) and 79 of 373 test files (707 lines).

TODO items that assume the cylinder and need rewording, holding or retiring (the TODO thread sweeps `ring`, `layer`, `slot`, `slab` and `wedge` for the rest): MAP.151 (151a and 151c are built on rings; 151b is unmerged), MAP.152, GEN.101, GEN.24, ADM.29, ADM.30, and anything keyed to the 243/81/27/9/3/1 ladder.

## 5. How the work should run

Switch structure (2026-10-15): DB.32, GEN.227 and UX.111 are PRs into a shared integration branch `claude/spherical-switch`, owned by the Spherical Sector Coordinates thread; one final PR goes to main after Boss's word. MAP.205 and MAP.206 are built in parallel on their own branch and merge after that.


One dedicated thread, "Spherical Sector Coordinates", owns the grid (items 1 to 6). It is the only migration lane while this runs, so the Database lane's DB.2x items queue behind it unless they touch other tables. Maps, Generation and UX lanes stay stopped until their items below are open. The Test Suite thread owns the TEST items. Each item merges on a green `scripts/test.py changed`; the full suite runs when item 6 finishes.

Rules that carry over: no backward compatibility, no shims; PRs add a `changes/` note and never bump the version; no `.major.md` note without Boss's word; update.sh refuses while jobs are running.

## 6. Work items (TODO assigns the numbers; order is the build order)

| # | Working title | Lane | Needs | Done when |
|---|---|---|---|---|
| 1 | Spherical grid module: shells, bands, wedges, caps and the core ball | Generation | none | `planetgen/galaxy/spherical.py` (pure geometry) with the rules in section 2: layout per shell, cell bounds, volume, `address_at` and `cell_center`, face neighbours, local frame, a cap-aware `SectorCell` sampler; shells 1 to 3 tabled; tests for exact tiling (cell volumes sum to the ball, and random points inside and outside the outline each land in exactly one cell, with no gap or overlap at shell, band, wedge, cap or core boundaries), round trip, counts near 23.93 billion, and the measured arcs |
| 2 | Outline and bounds on shells | Generation | 1 | `skeleton.py` builds the outline as a band range per shell (and the column view); `GalaxyBounds.contains` and `random_address` use it; cells counted against the cylinder's |
| 3 | Addresses: designation, sector ID, seeds, object ID | Generation | 1 | `provisional_sector_designation`, `sector_uid`, the 40-bit birth-sector field and the per-sector seed string use shell, band, wedge; width checks added to the plan step; GEN.69/170 rules unchanged otherwise |
| 4 | Schema: rename and migrate | Database | 1, 2, 3 | A new Alembic revision renames `ring_index`, `layer_index`, `ring_slot_index` to `shell_index`, `band_index`, `wedge_index` in the seven tables, deletes galaxy-placed sectors, systems and phenomena, rebuilds the outline, adds the old-schema fixture; merged only with Boss's go-ahead for the delete |
| 5 | Generation paths on shells | Generation | 1 to 4 | `run_plan`, `run_galaxy`, `bright_stars`, `phenomenon_scatter`, regions, spans and fills (GEN.24 core by whole shells, ADM.29/30 as shell and band ranges, GEN.101 fill order) use the new grid |
| 6 | Reads and web wording | Database, UX | 4, 5 | `query.py`, `store.py`, the API, the Galaxy and Sector pages, URLs and bookmarks say shell, band, wedge; "layer" and "slab" mean nothing any more; pages and API keep working on a regenerated small galaxy |
| 7 | Galaxy Map ladder for shells | Maps | 6 and the ladder card | A design note for nested blocks on shells (master rules in longitude and latitude, caps, stage keys, `galaxy_block_stats`), then `drill.py`, `region.py`, `window.py`, the block stats table and the migration |
| 8 | Galaxy Map drawing | Maps | 7 | `galaxyprisms.js`, `galaxystageview.js`, `galaxymap3d.js`, tile cuts and picking draw and pick shell cells and caps; slab slider becomes a band control |
| 9 | Tests | Test Suite | each of 1 to 8 | The 79 affected test files rewritten as each item lands; a MISSED line is a selection bug |
| 10 | Docs | Docs | 6 | `galaxy-coordinate-system.md` rewritten, the cylinder design archived, README, INSTALL and the docs sweep updated |

After item 6 Boss can regenerate a galaxy and look at it before Maps work begins.

## 7. Risks

- **Heavy polar cells.** A cap near the core holds about 2.5 times today's heaviest cell. If generation time for those cells hurts, the cap wedge count (card 2) or a cap fill cap is the lever.
- **Block nesting.** Item 7 needs a latitude master rule on top of the longitude one. A first try with odd 3^z masters gave band arcs of 0.89 to 1.11 of 4 pc [C]; longitude masters reach to within 4 to 7 degrees of the poles. Not designed.
- **Wide change set.** Items 1 to 6 touch four lanes' code; the order above keeps the tree buildable by keeping the pure geometry (item 1) first and the delete last.
- **Silent meaning change.** "Layer" and "ring" must not be reused for the new indexes; the rename in item 4 is what makes a missed spot fail loudly.
