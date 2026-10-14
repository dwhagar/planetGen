# Math and library efficiency: where the heavy math runs, what already uses optimized libraries, and what to change

A survey of the generator's numerical code against the libraries it could use: which hot spots are still
plain Python, which already run on NumPy or SciPy, what each candidate fix measured, and what each new
library would cost to install and keep. Measured on a real two-sector fill and a real plan, with the fixes
prototyped outside the repository and run end to end. No code in the repository was changed.

Informs: GEN.105, GEN.106, GEN.115, PERF.31, PERF.46, PERF.72, GEN.198, GEN.199, OPS.20

Status: research, 2026-10-14; no generator code was changed. Everything marked "recommend" is a proposal for Boss.

Evidence tags: [C] measured in this research (scripts under `/mnt/project-files/research/scripts/mathperf/`);
[S] read in the repository or on a package page; [R] recalled and unconfirmed.

## Method and limits

- One 4-core Linux container (Ubuntu 24.04, the VPS's size and OS), MariaDB 10.11 and Redis on the same box,
  Python 3.13.16, NumPy 2.5.4, SciPy 1.18.1, astropy 8.0.1, numba 0.68.0 [C]. The project's minimum is
  Python 3.9 and apt's NumPy on 24.04 is 1.26 [S][C]; speeds on the VPS's Python 3.12 will differ a little.
- The workload is the PERF.31 benchmark's: seed `...BEEF`, plan to ring 12 (2,041 layers), then `galaxy
  --ring 12 --limit 2 --yes` with one worker (1,916 systems, 17,284 bodies that get a sector path) [C].
  Profiles were taken on 8.0.9xx (`2c135ba`); the end-to-end comparisons on 8.0.986 (PERF.73 merged).
- Stage seconds come from the run's own stage records, not from cProfile, which inflates code made of many
  small calls about 2.5 times; profile seconds are used for shares only. Timings are one run each and
  vary by about 10% between identical runs (the baseline "Generate the sectors" stage took 46.3, 48.1 and
  52.5 s on three runs), so only differences well above that are claimed.
- The prototypes live in `/mnt/project-files/research/scripts/mathperf/` and are run with `patched_run.py`,
  which patches them into a normal `planetgen` run. Nothing here was run on the VPS, at the default
  galaxy's scale, on MySQL 8.4, or with more than one worker.

## The answer in short

1. **The numeric core is already on NumPy and SciPy where it was looked at.** The phenomenon scatter is
   vectorised (PERF.63), routing uses `cKDTree` (NAV.52) and A* on a k-nearest graph, Kepler's equation has a
   vectorised Halley solver, nebula meshes use scikit-image's marching cubes. No library the project lacks
   would speed those up.
2. **Three places are still plain Python and cost real time, and one of them needs a new library.**
   - The sector-path integrator (`physics/sector_path.integrate_path`) is 1.0 ms a path; 17,284 paths are
     the "Save the sector paths" stage, 34 s of a 92 s fill, and the orbital update (the 8.0 headline) runs it
     again for every body that moves. A numba port runs it in 34 to 46 us (22 to 31 times); with it the stage
     fell from 34.2 s to 13.4 s, the rest being database writes. A NumPy batch version managed only 1.8 times.
   - `SpaceSector.nearest_neighbors` sorts every other system for each system. `scipy.spatial.cKDTree`, already
     required, answers it 190 to 590 times faster with identical neighbours; the generate stage lost about 10%.
   - The bright and massive star scatters evaluate the galaxy density one point at a time (556,000 calls, 7.9 us
     each); the NumPy form the phenomenon scatter already uses is 47 times faster per point.
3. **numba is the one library worth adding, and only for the integrator and later kernels like it.** It costs 196
   MB installed, Python 3.10 or later, and a compile cache. Everything else recommended below needs no new
   dependency.
4. **Two libraries looked at are not worth adding:** `orjson` (JSON encoding is 1 to 4% of a tile request) and
   `mysqlclient` (no faster on inserts, GPL).
5. **One library could come out:** astropy is imported for about twenty constants; reading them from a table
   saves 280 ms and 29 MB per process start and 53 MB of install.

## What already uses optimized libraries (checked, no change needed)

| Area | How it is done | Evidence |
|---|---|---|
| Phenomenon scatter | NumPy arrays per layer; density, angle bins, positions, masses vectorised | `generation/phenomenon_scatter.py` `_vector_rows` [S]; 0.32 s of a 16 s layer pass is the random draws [C] |
| Routing graph | `cKDTree` k-nearest, A* with `heapq` | `galaxy/nav_graph.py` [S]; one A* query visits few nodes, `csgraph` would not beat it [R] |
| Kepler's equation | vectorised Halley (`solve_eccentric_anomalies`), `brentq` fallback | 1.36 us/orbit against 20.9 us scalar [C]; not yet called by production code, only by tests [S] |
| Nebula shape | scikit-image marching cubes, imported lazily | `galaxy/nebula_shape.py` [S]; `import skimage.measure` costs 12 ms [C] |
| Black hole spin | `scipy.special.betaincinv` | `physics/spin.py` [S] |
| Constants | astropy, read once at import | see "Constants" below |

## Findings

### 1. Sector paths: 1.0 ms a path, 22 to 31 times faster in numba

`integrate_path` is a leapfrog integrator written with tuples: per step it loops over up to 32 point masses
twice (step size, acceleration), tests the cell with `sector_address_at`, and a bisection of 40 more tests finds
the exit. Every one of the 17,284 paths came out as two knots (a straight line within the 0.2% tolerance), at
about 90 steps each [C]. Where the 1.0 ms goes: in a profile of 300 paths `_step_size` and `_acceleration` (31 masses by
90 steps) took 57%, and `relevant_masses` measured alone is 180 us (18%); the cell test and bisection are the rest [C].

| Variant (8,678 paths of one sector) | ms/path | Against shipped |
|---|---|---|
| V0 shipped | 1.02 to 1.08 | 1 |
| V1 memoise `ring_sector_count` and friends (`functools.lru_cache`) | 0.92 to 1.00 | 1.1 |
| V3 NumPy batch: all bodies advance together, 192 lock-step iterations | 0.53 | 1.8 to 2.8 |
| V2 numba integration, NumPy mass selection | 0.21 | 5 |
| V2b numba mass selection, integration and spline check in one call | 0.034 to 0.046 | 22 to 31 |

- The NumPy batch is the best that needs no new library and it is not worth its complexity: the bodies
  need different step counts, so the batch is mostly masked work.
- A pure-Python shortcut exists: 69% of bodies have a summed deflection bound under the tolerance, so a
  straight line would do for them [C]. That caps the saving at about 3 times and changes how the exit knot is
  found; numba keeps the algorithm as it is.
- **End to end** (`patched_run.py --kdtree --numba`): "Save the sector paths" 34.2 s to 13.4 s. What is left is
  the database: for each path the delete-by-key, the insert and the knot insert are separate statements [S].
- **Agreement.** Without `fastmath`, numba uses the same arithmetic in the same order. 34,566 of the 34,568 saved
  knots are bit-identical to the shipped run; the other two differ by at most 7e-15 pc in position and 4e-16 in
  relative time. An exit time agreed to 0.0 in the micro-benchmark. The prototype falls back to the shipped code
  for any path needing more than two knots (none did here). Tests that compare saved paths should use a
  tolerance, not equality.
- **Scale.** One million bodies is about 1,000 s of single-core integration in Python and about 40 s in numba.
  The orbital update recomputes the paths of every sector holding a body that moved
  (`cli/orbits.py`, `db/sector_paths.compute_sector_paths`) [S], so this is the cost that grows with the 8.0
  headline, not a one-time fill cost.

### 2. Neighbours: sort-everything against a k-d tree

`SpaceSector.nearest_neighbors(entry, count)` sorts all other entries by `entry.distance_to` for every entry
saved: N^2 log N, each `distance_to` going through two `position` properties (`SpatialPosition3D.get_coordinates`,
0.43 us). In the fill profile that chain was 8.3 million `get_coordinates` calls, 11% of profiled time [C].

| Systems in the sector | Shipped | `cKDTree` | Speed-up |
|---|---|---|---|
| 300 | 0.22 s | 3.2 ms | 68 |
| 1,000 | 2.28 s | 3.9 ms | 587 |
| 1,900 | 9.22 s | 47.9 ms | 192 |

The benchmark sectors hold about 960 systems; denser inner rings hold more, and the shipped cost is quadratic.
End to end the "Generate the sectors" stage went from 46 to 52 s (three baseline runs) to 42.7 s with the tree,
about 10%. The neighbour lists written to the database were identical in all 51,861 rows. The placement loop
(`sector.py` `_massive_neighbors`, 919,000 `distance_between` calls) has the same shape and can use the same tree.

### 3. Galaxy density, one point at a time

`bright_stars._ring_bins` and `layer_expected_stars` call `density.population_densities` per angle bin: 556,000
calls in the plan, 7.9 us each (4.4 s) [C]. `density.total_density_array`, built for the phenomenon scatter,
gives the same sum at 0.17 us a point (47 times) [C]. A per-population array version is a copy of that function
returning the four terms. Share of the plan: about 4 s of the 15 to 20 s the massive and bright star stages take;
the rest of those stages is waiting on the database (20.7 s of the 82 s profile sat in socket reads) [C].

### 4. Random draws can be bulk with no change to any galaxy

`draw.Stream` builds everything on `random.Random.random()`, so every galaxy depends on that exact sequence.
CPython's `random` and NumPy's `MT19937` are the same generator with the same 624-word state and the same
53-bit float, so the state can be handed to NumPy, a batch drawn in C, and the state handed back. The floats are
identical, and so is the generator's state afterwards; `getrandbits(63)` (two draws, shifted) also matches [C].

| Draws in a batch | Python loop | Bulk | Speed-up |
|---|---|---|---|
| 6 | 1.7 us | 102 us | 0.02 (slower) |
| 600 | 99 us | 113 us | 0.9 |
| 6,000 | 581 us | 141 us | 4 |
| 60,000 | 5.7 ms | 0.46 ms | 12.5 |

The 102 us is the fixed cost of moving the state. It pays off above about 1,000 draws a call, which only the
scatter layers reach (`_uniforms`: 6 per object; `rng.getrandbits(63)` per object). At the benchmark's scale the
draws are 0.3 s of 16 s, so this is a default-scale concern (hundreds of millions of rows), not a today concern.
No reseed is needed, unlike the move to a different generator.

### 5. The orbital update's per-row Python

`advance_galactic_positions` moves each due system with a chain of small functions per row (turn, rotate,
`local_to_galaxy_pc`, `sector_address_at`, `galaxy_to_local_pc`, `classify_octant`). The chain costs 25.7 us a
row; the same maths on NumPy arrays 3.8 us a row (7 times), and the refiled sector address agreed on all
100,000 rows [C]. The SQL side (`UPDATE ... executemany`) costs about as much per row again, so the whole update
gains about 3 to 4 times, not 7. Planet and moon motion is already a set-based SQL `UPDATE` [S].

### 6. Smaller costs in the same code

| Item | Cost | Fix | Gain |
|---|---|---|---|
| `geometry.ring_sector_count` and `_ring_centerline_edges` (3.2 and 6.5 million calls in the fill profile) | 1.15 us | `functools.lru_cache` | 0.064 us; 8% of the integrator, a few % of the fill |
| `util.checks.finite_domain` wrapper (786,000 calls) | 0.96 us against 0.08 us bare | PERF.46: Boss keeps it on | none proposed |
| pymysql parameter escaping on the sector save (5.7 million `escape` calls; pymysql as a whole is 8% of profiled time) | 8.6 us a 13-column row | numeric-only formatter as PERF.73 did for scatters | 2.9 us a row, 3 times |
| Kepler, scalar `brentq` | 20.9 us | vector solver where a batch exists | 1.36 us (15 times); a numba Halley loop 0.15 us (138 times) but not needed yet |

### 7. Libraries looked at and not recommended

- **orjson.** Encoding a Galaxy Map tile of 120 to 700 kB takes 0.2 to 2.8 ms against 1.7 to 13 ms for `json`
  (4 to 10 times faster), but the tile query itself is 37 to 660 ms [C]; encoding is 1 to 4% of the request.
- **mysqlclient (the C driver).** `executemany` of 200,000 13-column rows took 5.0 s against 4.4 s for pymysql
  (the server, not the driver, is the limit), `fetchall` 1.6 s against 3.0 s [C]. It is GPL-2.0-or-later, which
  the CC0 project has not accepted (orbital-solvers-and-integrators.md section 4.3), and it needs a compiler and
  `libmariadb-dev` to build [C]. `LOAD DATA LOCAL INFILE` was 2.7 times faster on inserts but needs
  `local_infile` on both ends, which PERF.73 already declined.
- **Cython or Pythran for the integrator.** Likely as fast as numba [R], but needs a compiler and Python headers
  on the VPS and a build in `update.sh`. Not measured.
- **A NumPy batch integrator.** 1.8 times; see finding 1.

### 8. Start-up cost of the imports

Measured one fresh interpreter per import, interpreter alone 11 MB [C]:

| Import | Time | Resident memory |
|---|---|---|
| `numpy` | +82 ms | 25 MB |
| `scipy.optimize` (for one `brentq` and `newton`) | +471 ms | 69 MB |
| `scipy.spatial` | +362 ms | 59 MB |
| `astropy.constants` (for about 20 constants) | +279 ms | 40 MB |
| `numba` | +307 ms | 86 MB |
| `nltk` (for one word list) | +1,102 ms | 108 MB |
| `planetgen.generation.run_plan` | +1,363 ms | 126 MB |
| `planetgen.web.app` | +2,347 ms | 143 MB |

The scatter workers now pre-import before they fork, so this matters for each command line start, each
mod_wsgi process and each worker restart, not for each job. astropy on disk is 41 MB, with pyerfa 3 MB and the
IERS data 8.6 MB; `physics/constants.py` is the only user [S]. The constants are fixed numbers (CODATA and IAU),
so a table in the repository (generated once from astropy, checked against it by the existing
`test_astropy_constants.py` with astropy kept in the `test` extra) removes the dependency and the import.
`nltk` is only used for the `words` corpus; reading the same file directly would drop 1.1 s and 97 MB, but
that is a names question, not a math one, and is listed only so it is not lost.

## Dependency cost of each candidate

| Library | Role | Installed size | Licence | Fit with install and update |
|---|---|---|---|---|
| **numba 0.68.0 with llvmlite 0.50.0** (the one to add) | integrator kernels, later vector kernels | numba 24 MB and llvmlite 172 MB in a venv; wheels 3.8 and 59.7 MB [C][S]; +75 MB resident per process that imports it, +307 ms | BSD; llvmlite BSD-2 with Apache-2.0 and the LLVM exception [S] | No `python3-numba` on Ubuntu 24.04 (apt has only `python3-llvmlite` 0.42) [C], so pip, through the NumPy stack group in `scripts/install-python-deps.sh`. Needs Python 3.10 or later and NumPy below 2.6 [S]; the project's floor is 3.9. First compile of the integrator kernels 5.9 s, from the disk cache 0.27 s [C]. Apache's code tree is read-only to the web user, so set `NUMBA_CACHE_DIR` (say `/var/cache/planetgen/numba`, made by `create-cache-dir.sh`) and warm it in `update.sh` as `warm_map` is. |
| scipy.spatial.cKDTree | neighbours, placement | none: scipy is already required | BSD | none |
| orjson 3.13.0 | JSON | wheel about 0.1 MB [S] | MPL-2.0 and Apache-2.0 or MIT [S] | pip wheel; not recommended |
| mysqlclient 2.3.0 | MySQL driver | needs `libmariadb-dev`, `pkg-config`, a compiler [C] | GPL-2.0-or-later [S] | apt has `python3-mysqldb` 1.4.6 [C]; not recommended |
| Cython 3.3.0 | compiled kernels | wheel 3.4 MB [S] | Apache-2.0 | compiler and headers on the VPS; not measured |

Boss's rule of no compatibility shims argues against a pure-Python fallback beside the numba kernel: make
numba a hard requirement and raise the minimum Python to 3.10 (the VPS runs 3.12). Memory is the VPS's tighter
limit with MySQL on the same box: importing numba adds 75 MB to each process that does, so only the processes
that integrate paths (the settle step, the orbit update, workers) should import it, not the web app.

## Ranked recommendations

| # | Change | Gain | New dependency | Risk |
|---|---|---|---|---|
| 1 | `nearest_neighbors` and the placement distance loop on `cKDTree` | 190 to 590 times on the call, about 10% of the generate stage, removes the quadratic term in dense sectors | none | Low: identical neighbours in all 51,861 rows |
| 2 | Sector-path integrator on numba (selection, integration, spline check in one call) | 22 to 31 times; settle 34.2 s to 13.4 s; the orbital update's main cost | numba, llvmlite, Python 3.10 | Medium: new dependency and cache; paths agree to 4e-16 |
| 3 | Per-population density as NumPy arrays in the bright and massive star scatters | 47 times on 556,000 calls, about 4 s of a 15 to 20 s pair of stages | none | Low: same formula as `total_density_array` |
| 4 | `lru_cache` on the ring and slot count functions | 18 times on the call | none | None |
| 5 | Sector save: numeric-only insert formatter (PERF.73's approach) for the sector fill's tables | 3 times on the formatting, which is about half of pymysql's 8% of the fill | none | Low if limited to numeric columns |
| 6 | Galactic move as NumPy batches in the orbital update | 7 times on the maths, 3 to 4 on the whole update | none | Low: refiled addresses agreed on 100,000 of 100,000 rows |
| 7 | Constants table instead of astropy | 280 ms and 29 MB a process, 53 MB of install | removes one | Low: values identical, checked by a test |
| 8 | Bulk draws for the scatter layers (state hand-over) | 4 to 12 times on the draws; matters at default scale | none | Low: draws and state identical, no reseed |
| 9 | Use the vector Kepler solver wherever updates are batched | 15 times against scalar | none | None |

Items 2 and 6 belong with the orbital update (GEN.105 chain) and should be built together, since both rewrite
the per-body loop; 1, 3, 4, 5 and 7 are independent small changes.

## Open questions for Boss (default stated)

1. Accept numba as a hard requirement and raise the minimum Python to 3.10? Default: yes; the VPS is on 3.12 and
   the alternative is a Python fallback, which the no-shims rule rules out.
2. Where may the numba cache live? Default: `/var/cache/planetgen/numba`, owned by the web user, warmed by `update.sh`.
3. Drop astropy for a constants table? Default: yes.

## What this does not predict

- **Default scale.** All timings are from a ring-12 galaxy of 1,916 systems. The quadratic term in finding 2 and
  the draw volume in finding 4 grow with scale; the numba gain is per body and should hold.
- **The database.** After the fixes in findings 1 and 2 the sector-path stage is 13.4 s of which nearly all is
  writes; the next gain there is batching the statements, not faster maths.
- **Workers.** Everything ran with one worker. With two workers the fill speeds up and the settle step does not
  (it is one process) [S]; numba matters more as the fill scales out.
- **mod_wsgi.** Importing numba inside Apache processes was not tried.
- **Python 3.9 to 3.12.** Speeds will differ; the ratios between variants should not.

## Proposed items (the TODO thread assigns IDs)

See `/mnt/project-files/research/handoff/math-and-library-efficiency.md`.

## Evidence notes

- The [C] numbers come from `paths_baseline.py`, `paths_variants.py`, `paths_batch_numpy.py`, `neighbours.py`,
  `density_bench.py`, `rng_bulk.py`, `galactic_move.py`, `kepler_bench.py`, `tile_json.py`, `db_driver.py`,
  `import_cost.py` and `patched_run.py` in `/mnt/project-files/research/scripts/mathperf/`.
- End-to-end stage seconds: baseline Generate 52.5, 48.1 s, Settle 34.2, 35.6 s; with `cKDTree` Generate 42.7 s,
  Settle 35.0 s; with `cKDTree` and numba Generate 45.1 s, Settle 13.4 s. The Generate difference between the
  patched runs is within the run-to-run spread; the Settle difference is not.
- Output comparison used the saved `nearest_systems` and `sector_path_knots` rows of the baseline and patched
  databases. `planetgen fingerprint` (about 62 s a database here) gave the same plan and region digests for the
  baseline and the `cKDTree` plus numba run; whether the digest covers sector paths was not checked.
- [R] items: that a Cython build matches numba's speed; that `csgraph` does not beat A* for one query.
