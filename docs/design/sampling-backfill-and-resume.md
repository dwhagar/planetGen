# Sampling, backfill, void filling, resumable runs and directives

How to draw stars for a large volume without visiting every cell (the exact
block-first draw that replaces GEN.42's "drop sectors by probability"), the
tests that prove a draw exact and reject the drop-and-boost form, what the
backfill costs, the arithmetic behind GEN.102 (fill the voids) and its
recommendation, a resumable-run design (PERF.29, PERF.30), and the sampling
rules for generation directives (GEN.96) and random neighbourhood centres
(GEN.97). The grid is in [galaxy-coordinate-system.md](galaxy-coordinate-system.md),
the density in [galaxy-disk-density.md](galaxy-disk-density.md), the per-unit
seeds in [reproducible-galaxies.md](reproducible-galaxies.md) section 3 and
[generation-determinism.md](generation-determinism.md), and fill order and
region enumeration in [fill-order-curves-and-core.md](fill-order-curves-and-core.md).

Informs: GEN.40, GEN.41, GEN.42, GEN.43, GEN.96, GEN.97, GEN.99, GEN.102, PERF.18, PERF.29, PERF.30 (and the GEN.39 rule that one seed gives the same stars wherever a run lands)

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [C] computed in this research and re-runnable (scripts named
under Sources), [S] read from the repository or a downloaded package, [R]
recalled and unconfirmed. The research environment could read only
search-result text, not papers, and its web-search budget was spent, so every
claim about the outside world is [R] and listed under Evidence notes. Another
thread is tuning the black hole and neutron star rates; nothing here touches
those.

## Decisions already taken

- **GEN.40 (Boss, 2026-10-01 22:22Z):** use the star density calculation to
  weed out sectors so "the back fill would literally have less work to do",
  without over-filtering, "because bright stars in odd places is really neat and
  realistic". **PERF.18 (2026-10-01 22:13Z):** parallelise the backfill on the
  work queue.
- **GEN.44 (Boss, 2026-10-02 03:25Z):** every sector stores its backfill level
  (-1 untouched, the lowest L_sun reached, 0 generated); built, PR #425.
- **GEN.96, GEN.97, GEN.99 (Boss, 2026-10-03 05:38Z):** an Override button with
  directives (density, at least x stars of a type, at least x habitable worlds);
  "x number of random neighborhoods"; a nebula that needs certain stars triggers
  a volume backfill down to 750 L_sun with skewed probabilities, and existing
  backfilled stars are regenerated in place with a similar luminosity if
  possible, otherwise left alone.
- **GEN.102, PERF.29, PERF.30 (Boss, 2026-10-07 11:47Z):** "Investigate filling
  all void space (space that plots to have a 0 or nearly 0 star density) at once,
  since those sectors have the least in them, we can afford to fill them"; record
  which sectors "might be only partially filled by identifying what runs haven't
  been made yet"; finish a block or sector "on next start".
- **GEN.39 (Boss, 2026-10-02):** a unit's numbers depend only on the galaxy seed
  and its address, never on run order or worker count.

## Findings in short

1. **GEN.42's drop-and-boost pass is not needed and the literal version is
   wrong.** Dropping 90% of cells and boosting survivors by 10 keeps the mean
   but gives dispersion 1.136 (theory 1.144) and a multiplicity chi-square of
   743,000 on 14 degrees of freedom. Exact methods pass. The saving Boss wants
   exists exactly as a block-first draw (section 2); GEN.43 then holds by
   construction.
2. **At today's run sizes speed is not the problem.** A 100 ly backfill (1,800
   sectors) costs 2 ms per sector warm, 0.06 ms of it for draws. The first
   backfill in a process costs 5 s (40 s with no ceiling) to build
   `_bright_table`, and a forked RQ work horse loses that cache every job
   (PERF.18).
3. **GEN.102: nothing inside the outline is near zero density.** The default
   galaxy has 23,920,622,421 sectors and about 3.0 x 10^11 systems; the lowest
   density is 0.069. The 12.4% of sectors with under one expected system hold
   0.77% of systems and would take about 5 years of one worker. Do not
   pre-fill (section 5).
4. **PERF.29 and PERF.30:** paths and the backfill around a run take "which
   sectors" from the dead process's `started_at`. Add a transactional pending
   mask and resume by re-running the unfinished `generation_runs` command
   (section 7).
5. **GEN.96:** 24.5% of systems have a habitable world. Star-count directives are
   exact by truncated Poisson inversion; habitable-world directives need bounded
   whole-sector rejection with the attempt in the seed key (section 8).
6. **GEN.97:** plain rejection on density is unbiased (about 77 proposals per
   centre); density-weighted centres cost 14 times more to fill (section 9).
7. Earlier figures this supersedes (500 L_sun, 61 million stars, 88 billion
   systems, the timing report's "26B" column, "cost per cell is unmeasured") are
   in section 10.

## 1. What is exact and what is not

Poisson facts [R, Kingman, *Poisson Processes*, 1993; Devroye, *Non-Uniform
Random Variate Generation*, 1986]:

- Counts in disjoint regions are independent and Poisson with mean equal to the
  integral of the intensity: cell i holds Poisson(lambda_i) stars, lambda_i = E x
  rho_i x (share of stars in the band), E = 6.306 systems per sector at density 1.
- **Marking:** labelling each point by type with probabilities p_t gives
  independent Poisson processes with intensities p_t x lambda. This makes "at
  least N stars of a type" exact (section 8) and GEN.99 safe (section 6).
- **Given a count n, the points are i.i.d. with density proportional to the
  intensity,** so a block's count can be drawn first.
- **Thinning** (Lewis and Shedler, 1979, *Naval Research Logistics Quarterly* 26)
  [R]: simulate a homogeneous process at an upper bound and keep a point at x
  with probability lambda(x) / bound. Exact; cost is proportional to bound x
  volume, so the bound must be tight. **Skip-ahead** is the same idea with
  exponential gaps between candidates.
- **Bernoulli is not Poisson.** Bernoulli(lambda) never gives two stars in a cell;
  Poisson gives P(2) = lambda^2 / 2. `backfill_cells` draws Poisson counts, so
  Poisson is what a faster method must match.
- **Bounds, not sums, make it cheap.** The project already has one:
  `skeleton.bound_relative_density_at(shape, r, z)` is the maximum over angle of
  the density at a ring and layer and cuts the galaxy's outline
  ([galaxy-disk-density.md](galaxy-disk-density.md) section 3). Inversion of
  cumulative counts and Walker/Vose alias tables (Walker 1977; Vose 1991) [R]
  were considered; the first ties results to a region's total, the second needs
  a set-up pass over every cell, so neither suits a sparse region.

## 2. The block-first draw

### 2.1 Procedure

For one band (`canonical_bands`, eight to a decade) and one block, a (ring,
layer, 256-slot chunk):

1. `bound = E x bound_relative_density_at(shape, ring radius, layer z) x s_max`,
   `s_max` the largest population share of the band, taken as at least the
   `MIN_RELATIVE_DENSITY` floor (GEN.78). Valid because every cell of the block
   sits at the ring radius and layer height the bound is evaluated at.
2. Draw the candidate count n ~ Poisson(bound x cells in the chunk) from one
   keyed uniform by inversion.
3. Each candidate gets a uniform slot in the chunk and an acceptance uniform u;
   keep it when `u x bound` is below the cell's true rate (`E x sum over
   populations of density_p x share_p`). A block with no candidates costs one hash.
4. Pick a kept candidate's population last, in proportion to the per-population
   rates, and draw its star from a stream keyed by (block key, candidate index).

The pre-placed scatter works the same way per ring and angle bin. **Choose per
block:** use the candidate path when expected candidates are under about 5% of
the block's cells, the per-cell path otherwise. Both are exact, and in the bulge,
where a quarter of cells hold a star, per-cell wins (the prototype was 10 times
slower there, section 4.2). One process per population with its own bound would
cut candidates 5 to 10 times.

### 2.2 What the draws must be keyed by

Block-first sampling meets GEN.39 if and only if:

1. The block grid is fixed globally, independent of the run's region; a region
   only chooses which blocks to read and filters candidates by cell and tier.
2. Nothing sequential crosses block boundaries.
3. Star parameters come from a stream keyed by (block key, candidate index).
4. The bound depends on the block, not the region, and float sums are not
   reduced in a partition-dependent order (`planetgen/util/detmath.py`,
   [generation-determinism.md](generation-determinism.md) section 3).

Tested in C: both block methods give bit-identical star lists processed in
index order, in random order, split among four workers, or restricted to a
sub-region (5 seeds each) [C, `determ.py`]. TEST.77's golden galaxy should
include one block-first backfill once it lands.

### 2.3 Cost of keyed randomness

CPU time per keyed uniform [C, `micro.py`, `micro2.py`]: `draw.Stream(str)`
(`random.Random(str)`: SHA-512 plus Mersenne Twister seeding, the current
per-cell, per-band cost) 9.5 us; `random.Random(int)` 7.8 us; SHA-256 to a float
1.3 us; numpy `Philox` constructor 14 us, then 38 ns per number; a splitmix64
finaliser in C 1.9 ns. A Mersenne Twister is the wrong tool for one tiny
decision per cell. In Python use SHA-256 for the one or two decisions per block
and a `random.Random` only for blocks with candidates. Philox and Threefry
(Salmon et al., "Parallel random numbers: as easy as 1, 2, 3", SC11, 2011) [R]
suit C or numpy. `exp` and `log` can differ in the last bit between machines, a
probability of about 10^-16 per decision of flipping an acceptance test, which
TEST.77 watches for.

### 2.4 Code sketch

Distilled from `proto_bf.py` (112 lines, on the project's real geometry); not
in the repository.

```python
def block_candidates(seed, band, ring, layer, chunk, cells, bound):
    key = sha256(f"{seed}:bf:{band}:{ring}:{layer}:{chunk}")   # one hash per block
    n = poisson_inv(bound * cells, uniform(key, 0))            # split means above 30
    for i in range(n):
        slot = chunk * 256 + min(int(uniform(key, 1 + 3 * i) * cells), cells - 1)
        yield slot, uniform(key, 2 + 3 * i), uniform(key, 3 + 3 * i)
# caller keeps (slot, u_accept, u_pop) when u_accept * bound < rate(slot, band);
# the population is the categorical on the per-population rates, using u_pop
```

`poisson_inv` must split a mean above about 30. The first C version capped
inversion at 60 events while block means reach 350; the harness showed a mean z
of -290 at once, which is why section 3's tests are the done-test.

## 3. Tests that reject drop-and-boost

### 3.1 Why the literal pass is detectably wrong

Keeping each cell with probability q and drawing survivors from Poisson(lambda /
q) is a zero-inflated Poisson. The mean is right, the variance is not:

```
var(count in a cell)   = lambda + lambda^2 (1/q - 1)
dispersion of a region = 1 + (1/q - 1) x sum(lambda^2) / sum(lambda)
P(2 stars in a cell)   = q x (lambda/q)^2 / 2 = lambda^2 / (2q)    # 1/q times the true lambda^2 / 2
```

For the test field below sum(lambda^2) / sum(lambda) is 0.01605, so q = 0.1
predicts 1.144 [C]; the harness measured 1.136 over 3,000 replicates. Stars
clump (survivors hold ten times their share). A cell's chance of exactly one
star falls (lambda x exp(-lambda/q) against lambda x exp(-lambda)) and of two or
more rises, so the pass keeps the expected count but not the distribution GEN.43
asks to preserve. The eye misses this at lambda of 10^-3; the tests below do not.

### 3.2 The tests

Field: 64^3 = 262,144 cells, lambda from 1.9 x 10^-4 to 0.086 (mean 1.14 x 10^-3,
total 300), 4,096-cell blocks, 3,000 replicates per run [C, `harness.py`]. Per
run: (a) replicate totals against Poisson(total): mean z, dispersion, chi-square
on 20 equiprobable bins; (b) per-chunk counts against expectation (128 chunks);
(c) multiplicity chi-square (0, 1, 2, 3 or more stars per cell) per half-decade
stratum of lambda against exact Poisson probabilities; (d) KS test of pooled
positions in cumulative-intensity coordinates.

| Method | Runs | p-values below 0.05 | Verdict |
|---|---|---|---|
| Per-cell Poisson (reference) | 4 | 1 of 20 | matches |
| Block skip-ahead (exponential gaps at the bound, thin) | 4 | 0 of 20 | matches |
| Block count-then-distribute (Poisson candidate count, uniform positions, thin) | 6 | 3 of 30, none repeated on other seeds | matches |
| Drop 90% of cells, boost survivors x10 (literal GEN.42) | 1 | 3 of 5 | **rejected**: mean right (z = +0.14, spatial p = 0.77), dispersion 1.136 (p < 0.001), multiplicity chi-square 743,172 on 14 |

The three exact methods gave 70 p-values, 4 under 0.05 (3.5 expected by chance).
Also checked [C, `trunc.py`]: truncated-Poisson inversion matches whole-count
rejection (chi-square 23.0 on 16, p = 0.12); marking gives independent type
counts (correlation -0.001); "at least K of type T" by truncated inversion plus
an independent Poisson for the rest matches joint rejection (means 8.085 and
8.094).

### 3.3 Done-test for GEN.42 and GEN.43

On the real `backfill_cells`, when block-first lands:

1. Stratified multiplicity and per-stratum chi-square against Poisson over 3,000
   or more replicates of a small region, including the lowest-lambda strata
   (GEN.43's "bright stars in odd places"); it has the power to catch 743,000.
2. Assert `rate <= bound` for every cell of the region, so a wrong bound fails
   loudly instead of biasing quietly.
3. Same star list in index order, random order, four workers, and as a sub-region
   of a larger one (section 2.2).
4. Compare with the existing per-cell path as well as with theory.

## 4. What the backfill costs, and the verdicts

### 4.1 Benchmark: per cell against block-first (C)

C, `-O2`, one thread, 16^3 blocks, mean lambda 1.14 x 10^-3 [C, `bench.py`]:

| Cells | Stars | Per-cell with density | Block skip-ahead | Block count + distribute |
|---|---|---|---|---|
| 8.8 x 10^5 | ~1,000 | 35 ms | 0.7 ms | 0.6 ms |
| 9.0 x 10^7 | 104,000 | 3.77 s | 17 ms | 14 ms |
| 1.02 x 10^9 | 1.17 x 10^6 | 42.6 s | 0.141 s | 0.123 s |

Block methods are 300 to 350 times faster at 10^9 cells because they evaluate
density only at candidates (1.48 million for 1.17 million stars, 79% accepted).
Scaled to 23.9 billion cells: about 1,000 s per-cell in C, 3.4 s block-first. In
Python the real `backfill_cells` costs 58 us per cell (three bands), about 16
days, dominated by `population_densities` (5.7 us), `relative_density` (2.8 us)
and stream creation (9.5 us per cell and band) [C, `micro3.py`].

### 4.2 The real backfill and a Python prototype

`backfill_cells` against `proto_bf.py` (block-first on the real geometry, with
`bound_relative_density_at`, `_densities` and `canonical_bands` unchanged). 100
ly sphere, tiers 100/250/500/750 L_sun, ceiling 1000 L_sun [C]:

| Place | Cells | Expected stars | Real mean (30 runs) | Prototype mean (400 runs) | Prototype dispersion / chi-square p | CPU real | CPU prototype |
|---|---|---|---|---|---|---|---|
| Solar radius (8.2 kpc, in plane) | 1,804 | 5.70 | 6.00 | 5.64 | 1.03 / 0.64 | 105 ms | 31 ms |
| Bulge edge (2 kpc) | 1,847 | 162.6 | 163.6 | 161.8 | 0.98 / 0.70 | 71 ms | 748 ms |
| Outer disk (13 kpc, 600 pc up) | 1,844 | 0.58 | 0.43 | not run | | 60 ms | 11 ms |

The prototype reproduces the expected count, a Poisson dispersion of 1 and the
radial distribution (shell p = 0.97 and 0.72). Work by centre [C, `tiers_tbl.py`]:

| Centre | Density | Cells | Expected stars | Cells per star | P(cell gets a star) |
|---|---|---|---|---|---|
| Outer disk | 0.14 | 1,844 (98 inside the outline) | 0.58 | 3,174 | 0.0003 |
| Solar radius | 1.55 | 1,804 | 5.70 | 317 | 0.0032 |
| Bulge edge | 33 | 1,847 | 162.6 | 11 | 0.088 |
| Bulge, 500 pc | 141 | 1,895 | 716.8 | 2.6 | 0.38 |

At the solar centre the outermost tier (50 to 100 ly, floor 750) is 88% of the
cells and 72% of the stars (1,584 cells, 4.09 stars of 5.70).

### 4.3 The real run against MariaDB 10.11

Scratch database, default plan, no ceiling, 4-core box at load 15 to 25 (wall
noisy, CPU steadier) [C, `dbbf2.py`]:

| Run | Sectors | Stars | Wall | CPU |
|---|---|---|---|---|
| First backfill in the process (cold) | 1,800 | 10 | 44 s | 40 s |
| Warm | 1,800 | 13 | 4.0 s | 2.7 s |
| Warm | 1,804 | 7 | 3.4 s | 2.6 s |
| Warm, bulge edge | 1,847 | 186 | 4.1 s | 2.8 s |

Warm, under cProfile (3.9 s): `backfill_cells` 2.3 s (58,118 `random.Random`
constructions, 32 per cell, because with no ceiling the bands run open-ended;
with a 1000 L_sun ceiling about 3 per cell and about 0.1 s) and 3,642 SQL
statements about 1.5 s. So **about 2 ms per sector warm**, about 0.8 ms of it
SQL, with one `sector_stats` row written per visited sector
(`lock_sector_stats`). The cold start is the surprise: `_bright_table` (an
`lru_cache` in `generation/star_population.py`) costs about 0.15 s per distinct
(band edge, population): 5 s with the ceiling (9 edges x 4 populations), 40 s
without (34 edges). RQ forks one work horse per job (`queue/redisqueue.py`
`worker_class`; `SpawnWorker` where there is no fork loses the cache too), so
every queued backfill job pays it again [S].

### 4.4 Verdicts

**GEN.41.** At the solar centre the backfill visits 1,804 sectors for 5.7 stars
(317 per star), at the bulge edge 11, at the core 2.6. A perfect pre-pass cuts
visits by 99.7% at solar density but saves only about 0.1 s of a 3.4 s run (the
rest is SQL and the ledger) plus the avoidable cold start. **GEN.42 as written:
no-go.** It is inexact, saves little at sphere sizes, and the large savings (the
ledger and the visits) come from exact block-first.

**GEN.42.** Replace with an optional exact block-first draw for a run over about
10^5 cells or one dominated by its sparsest tier. A literal drop pass would have
to be exact thinning of the candidate process, which is skip-ahead.

**GEN.43.** Satisfied by construction: every cell keeps its own mean, and
`MIN_RELATIVE_DENSITY = 1e-3` (GEN.78) plus the outline minimum of 0.069 mean no
cell has zero density. The test is section 3.3.

**PERF.18.**

- Make the task a chunk of many thousands of cells, not one 1,800-cell sphere:
  queue and fork overhead (tens of ms) equals a sphere's draw time. Discovery
  currently enumerates a 60-fold redundant cell list around every generated
  centre ([fill-order-curves-and-core.md](fill-order-curves-and-core.md) section
  3.5); enumerate once around the run's boundary.
- Pre-warm `bright_star_fraction` for every canonical band edge in the worker
  parent before it forks, or compute the tables at plan time and store them (they
  must be identical on every machine). This is the largest saving found.
- Determinism is already right: stars are keyed per (seed, cell, band).
- Draws (0.06 ms per cell) parallelise nothing worth having; SQL (0.8 ms per
  sector) and neighbour locks cap the gain at 1.87 times for 2 to 4 workers
  (timing report section 4). Measure again with PERF.31 before promising speed.
- In-request backfills should queue and return a job id; a web request should
  never pay the 5 s cold start.

## 5. GEN.102: filling the near-zero-density space

### 5.1 Arithmetic for the default galaxy

Outline: 2,041 layers (+-1,020), outer ring 3,763, **23,920,622,421 sectors**
(`candidate_sector_count`). E = 6.306, mean density 1.995, so about **3.01 x
10^11 systems**. Stars at or above 1000 L_sun: 2.70 x 10^7, matching the
26,925,945 the timing run wrote [C, `outline.py`, `void.py`, `void2.py`,
`galaxy_totals.py`: 4 x 10^7 sampled cells inside the outline, a numpy replica of
`relative_density` matching it to 10^-15].

Density quantiles inside the outline (1.0 is the solar neighbourhood): minimum
0.069, 1% 0.087, 5% 0.119, median 0.42, 95% 7.0. **Nothing is near zero.** The
outline keeps a cell if its density bound reaches 1/E = 0.1586, so the actual
density is at worst about 0.069; the 10^-3 floor is reached only outside the
galaxy, where no sector exists.

| "Void" definition | Share of sectors | Sectors | Systems in them | Share of all systems |
|---|---|---|---|---|
| density below 0.1 | 2.44% | 5.8 x 10^8 | 3.2 x 10^8 (0.56 per sector) | 0.11% |
| expected systems below 1 (density below 0.1586) | 12.4% | 2.97 x 10^9 | 2.3 x 10^9 (0.78 per sector) | 0.77% |
| expected systems below 2 | 40.4% | 9.7 x 10^9 | not computed | 3.9% |

Density below 0.1 sits at a median radius of 13.4 kpc (10th to 90th percentile
10.1 to 14.7) and height 0.44 kpc (0.08 to 0.96): the between-arms outer disk and
thick-disk fringe, not the halo. Below 1/E: 7.1 to 14.4 kpc, 0.1 to 1.4 kpc up.

### 5.2 Time and storage

Cost model from the timing report: about 40 ms fixed per sector plus about 16.4
ms per ordinary system (3.4 s for 205 systems, less the fixed part); 1.87 times
faster at 2 to 4 workers on its 4-core box. Storage about 40 KB per system (the
pre-placement plan's 3.5 PB for 88 billion systems; not re-measured) [R].

| Fill | Sectors | Systems | One worker | At 1.87 x | Storage |
|---|---|---|---|---|---|
| density below 0.1 | 5.8 x 10^8 | 3.2 x 10^8 | about 330 days | about 175 days | about 13 TB |
| expected systems below 1 | 2.97 x 10^9 | 2.3 x 10^9 | 5.0 years | 2.7 years | about 92 TB |
| the whole galaxy | 2.39 x 10^10 | 3.0 x 10^11 | about 185 years | about 100 years | about 12 PB |

Boss's premise ("those sectors have the least in them, so we can afford to fill
them") fails on the cost model: time is a fixed cost per sector, not per star. A
sparse sector costs 70 to 90 ms per system against 16 for an ordinary one.

### 5.3 Recommendation

Do not pre-fill voids. Keep lazy fill on visit (`ensure_sector_generated`,
`generation/run_galaxy.py`), which costs nothing until someone looks. If the
sparse rim is wanted for completeness, the only lever is a smaller fixed cost:
the timing report's levers (17.6 ms of serialised neighbour-lock steps, the
second uid pass, one registry upsert per sector) could plausibly halve the 40
ms, making the density-below-0.1 fill about 5 months of single-worker time. About
57% of those sectors are empty (exp(-0.56)), so "empty needs no row" would save
those writes, but then "unfilled" and "filled and empty" must be told apart (a
`sector_stats` row or an interval ledger). Define the threshold as expected
systems below 1, matching the existing rule `predicted_star_count >= 1`, not
"near zero".

## 6. GEN.99: backfilling a nebula's volume

A giant H II complex is at most about 200 ly radius (about 15,000 sectors,
[nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md)); at
58 us per cell that is about 1 s of draws, so exactness and "similar luminosity"
matter, not speed.

- **Skew the marks, never the intensity.** A cell's star count is Poisson(lambda)
  whatever the type law; draw each star's type from a tilted law p'(t) proportional
  to p(t) x w_t. Counts stay natural and only labels change. A weight of about 40
  to 50 on the required types turns the natural 8% O/B share of stars at 750
  L_sun and above [R, from the pre-placement plan] into about 80% (odds 0.087 to
  4).
- **Redrawing existing stars in place** conditions each on its own luminosity:
  draw the new type from p'(t | L within a tolerance of L_old) proportional to
  p(t, L) x w_t; if the support is empty leave the star alone, as Boss's rule
  says. Key the draw by (star seed, nebula id) so a second nebula over the same
  stars gives a defined result.
- The 750 L_sun floor uses the same band machinery, so a nebula larger than the
  sphere tiers can share the block-first draw.

## 7. Resumable runs: PERF.29 and PERF.30

### 7.1 What exists

A sector save is **one transaction** (`save_sector` inside `conn.batched()`:
sector row, systems, names, `sector_stats` row last, one commit by the caller), so
a half-written sector rolls back [S]. Everything after it is a separate pass:

| Step | Where its state lives today | Survives a crash? |
|---|---|---|
| sector, systems, phenomena | the sector row (atomic) | yes |
| bright-star level | `sector_stats.bright_level_sol` (-1 untouched, L_sun, 0 generated), written with the stars in one chunk transaction (`BACKFILL_CHUNK_SECTORS` = 200) | yes; `_draw_sector_bands` also wipes stray stars of a sector with no level |
| paths | `sector_paths` rows; which sectors need them is `created_at >= started_at` of the **running process** (`sector_ids_since`, `db/sector_paths.py`) | **no**: a restart gets a new `started_at` |
| backfill around the run | `backfill_after_run` uses `store.sector_centers_since(started_at)` | **no**, same reason |
| population | `population_state.scanned_planet_id` (a global watermark) | yes |

The gap is exactly the post-steps whose "which sectors" lives only in memory. The
run history exists (DB.6): `generation_runs` holds the seeds, version key and
`outcome` (`ok`, `failed`, `interrupted`), `generation_run_arguments` the command
line; a run killed hard leaves `finished_at` and `outcome` NULL [S, `db/store.py`].

Patterns [R]: idempotent units keyed by address; a watermark for ordered work
(the population pass); a transactional outbox (write the intent in the same
transaction as its effect); at-least-once delivery with idempotent handlers.
Minecraft stores a per-chunk `Status` so a half-generated chunk continues from
its status on the next load [R, names vary by version]; the bright-star level is
such a ladder, but PlanetGen's other steps are not (paths need neighbours'
final contents, population is global), so a bitmask fits.

**RQ semantics** [S, `rq` 2.12.0 wheel, `registry.py`]: a job whose worker died
stops heartbeating; `StartedJobRegistry.cleanup()` finds the expired entry and, if
`retries_left > 0`, calls `job.retry()`, otherwise files it as `AbandonedJobError`.
The project's `RQExecutor` enqueues with `retry=rq.Retry(max=TASK_RETRIES)`,
`job_timeout=-1`, and treats a dead worker as `WorkerDied`. Retry reruns the
whole job, so every task must be idempotent: sector fills skip existing
addresses and a backfill chunk rechecks its levels under the row lock. The lock
pins rq 2.12.0 for Python 3.10 and later, 2.8.0 for 3.9; only 2.12.0 was read.

### 7.2 Recommended design

1. **PERF.29: a `sector_pending(ring, layer, slot, mask, run_id)` table.** The
   sector save inserts its row in the same transaction, with a bit for each step
   still owed (paths, backfill tier, population, later steps); each step clears
   its bit in the same transaction as its effect. Rows exist only while work is
   owed, so the table stays small at 24 billion sectors, and "every partly filled
   sector and what it needs" is a `SELECT`. Where a step has a natural artifact
   (bright level, population watermark) derive it rather than store a copy; use
   the mask for steps with none (paths, the pending backfill request with its
   tier parameters). GEN.44's level rule stays the summary.
2. **PERF.30: re-run the unfinished command.** A run is unfinished when its
   `generation_runs` row has `finished_at IS NULL` or `outcome` is `failed` or
   `interrupted`; `--resume` re-runs its stored arguments and run seed (no new
   outcome value is needed). The fill skips existing sectors, then drains
   `sector_pending`. The address sequence is a pure function of the arguments
   (order changes only what exists when a run stops,
   [fill-order-curves-and-core.md](fill-order-curves-and-core.md) section 2.5), so
   a rerun reaches the same units. An integer frontier can report progress, but
   correctness should rest on per-unit existence checks because with parallel
   workers the done set is a prefix plus stragglers.
3. **Failure-mode tests:** kill between sector commit and settle; inside a
   backfill chunk (levels and stars stay consistent, as now); a forked work
   horse (RQ retries); a duplicate task (the unique address index and row lock
   make the second a no-op); a clock change (never use time as the frontier).
4. **Not checked:** whether the uid claim and name-registry confirmation (#657)
   sit inside the sector transaction. If either is a separate commit, a crash
   leaves an orphan reservation that the DB check (DB.8,
   [db-check-and-parity-repair.md](db-check-and-parity-repair.md)) should find.

## 8. GEN.96: directives

### 8.1 What the numbers say

On 6,000 systems from `StarSystem(SystemConfig())` [C, `hab2.py`]: habitable
worlds per system 0: 75.5%, 1: 18.95%, 2: 3.8%, 3: 1.2%, 4: 0.4%, 5: 0.1%; **at
least one: 24.5%**, mean 0.325. A system takes about 20 ms to generate; forcing a
habitable world (`HABITABLE_WORLD = True`, already in `StarSystem`) costs 29 ms and
never failed in 600 systems. A sector's count is a compound Poisson (systems ~
Poisson(n), each with the counts above); P(at least N) by exact recursion matches
simulation (0.327 against 0.3265 at n = 6.3, N = 3) [C, `hab_math.py`]:

| Expected systems in the sector | N >= 1 | N >= 2 | N >= 3 | N >= 5 | N >= 8 | N >= 10 | Mean worlds |
|---|---|---|---|---|---|---|---|
| 0.35 (sparse rim) | 0.082 | 0.021 | 0.0071 | 0.0007 | 0.0001 | 0 | 0.11 |
| 2.0 | 0.387 | 0.155 | 0.064 | 0.010 | 0.0008 | 0.0002 | 0.65 |
| 6.3 (density 1) | 0.786 | 0.530 | 0.327 | 0.104 | 0.014 | 0.0036 | 2.05 |
| 12.6 | 0.954 | 0.845 | 0.692 | 0.381 | 0.107 | 0.039 | 4.10 |
| 32 | 0.9996 | 0.997 | 0.989 | 0.939 | 0.738 | 0.548 | 10.4 |

Attempts for 99% success at per-attempt success p: 0.5: 7; 0.2: 21; 0.05: 90;
0.01: 459; 0.001: 4,603.

### 8.2 Methods

| Directive | Method | Exact? | Cost |
|---|---|---|---|
| must have density d | set the density multiplier | by definition | none (`--density` exists) |
| at least N stars (systems) | truncated Poisson by inversion | exact conditional law | O(mean), no rejection |
| at least N stars of type T | split the Poisson by type (independent counts): type-T count by truncated inversion, rest ordinary | exact conditional law | O(mean); needs a per-system forced type (`SystemConfig.STAR_TYPE` exists) |
| at least N habitable worlds | **bounded whole-sector rejection** | exact (a normal sector, given the event) | 1/P builds on average, K at most |
| same, fallback | forced top-up: swap habitable-free systems for `HABITABLE_WORLD = True` systems until N | no: over-represents forced systems | one build plus (N - H) x 29 ms |

Truncated inversion: `F = poisson_cdf(N - 1, lam); k = poisson_ppf(F + (1 - F) * u, lam)`
with u from the unit's stream, giving k >= N.

Rejection generates in memory and saves only the winner (a save is 40 ms plus
about 16 ms per system and must not be paid per attempt). At density 1: N = 3 3.1
attempts, 0.4 s; N = 5 9.6 attempts, 1.2 s; N = 8 71 attempts, 9 s; N = 10 278
attempts (K for 99%: 1,280), 35 s. Sparse sector (0.35 systems), N = 1: 12
attempts, 0.1 s; N = 3: 141 attempts (K 650), 1 s.

### 8.3 Recommendation

- Attempt 0 draws from the sector's plain unit seed, so no override means exactly
  the natural sector and a natural sector that already meets the directive is
  untouched. Attempt a >= 1 uses `SHA-256(galaxy seed || "sector:addr" || ":d:" ||
  H(canonical directive JSON) || ":a:" a)`.
- Try up to K attempts (default 200, user-settable), then the forced top-up for
  habitable-world directives, then report `met_naturally`, `met_after_k_attempts`,
  `met_forced` or `unmet` with the best achieved values.
- Store the attempt number, directive and mode with the sector (or as a GEN.59
  net-difference entry, [generation-determinism.md](generation-determinism.md)
  section 6): the sector is no longer a pure function of seed and address, and
  without the record `reproduce` (OPS.12) cannot rebuild it.
- Refuse an impossible request up front from the table above instead of burning K
  attempts.

Not recommended: generating each system conditional on a pre-drawn habitable
count. It is exact and cheaper, but the pmf depends on star type, age and the
nebula, so the table would be rebuilt for every setting.

## 9. GEN.97: N random neighbourhoods

- **Unbiased density-weighted centres.** Propose a cell uniformly (cells are
  within about 6% of equal volume; for exactness weight by pi x (2i + 1) x edge^3 /
  slots(i)), accept with probability (rho / rho_max)^gamma, rho_max = 154.2.
  gamma = 0 is uniform by volume (what `GalaxyBounds.random_address` does,
  [fill-order-curves-and-core.md](fill-order-curves-and-core.md) section 3.5), 1 is
  proportional to star density. Acceptance at gamma = 1 is 1.29%, about 77
  proposals per centre at about 10 us each. Never clip rho at a cap: it biases
  the weighting.
- **Cost depends on gamma, not on the sampler.** Systems in a 100 ly neighbourhood
  (1,887 cells) [C, `centres.py`]; time is 40 ms per sector plus 16.4 ms per system:

| gamma | Mean density at centre | Systems per neighbourhood | One worker |
|---|---|---|---|
| 0 | 2.0 | 2.4 x 10^4 | about 8 min |
| 0.25 | 4.1 | 4.9 x 10^4 | about 15 min |
| 0.5 | 8.7 | 1.0 x 10^5 | about 29 min |
| 1 | 28.1 | 3.4 x 10^5 | about 1.6 h |

  At gamma = 1 about 40% of centres land in the 1.7% of sectors with density
  above 20 (the bulge). Default gamma = 0: the least surprising reading and the
  cheapest.
- **Qualifying criteria**, in this order: inside the outline; the whole
  neighbourhood inside (3% of uniform centres have a 100 ly sphere crossing the
  edge, 7% at 250 ly); density at least 1/E (the project's "predicted star count
  >= 1"); optionally not within r of a filled sector (one indexed box query on
  `idx_sectors_center` per accepted candidate).
- **No overlap, by dart throwing.** Keep a hash grid of cell size 2r over accepted
  centres and reject a candidate closer than 2r plus a margin (Poisson-disk
  sampling; Bridson 2007, SIGGRAPH sketches) [R]. A 100 ly neighbourhood excludes
  9.7 x 10^5 pc^3 and random sequential packing fills about 30% of the volume, so
  it succeeds for up to about 4.8 x 10^5 centres (3.0 x 10^4 at 250 ly) galaxy-wide
  and 3.4 x 10^4 in the bulge alone. Cap attempts at 200 x N and report a shortfall.
- **Determinism.** Candidate j draws from a stream keyed by (galaxy seed, run
  seed, "neighbourhood", j), so the accepted list depends only on the seed and the
  arguments. If the filled-space test is on, record the chosen centres in the
  `generation_runs` row (or its arguments) so `--resume` and `reproduce` replay
  the same list.

## 10. Reconciliation with earlier documents

The earlier documents are not edited here except where stated.

| Earlier statement | Where | Now |
|---|---|---|
| "26B (x966)" column: 16 days scatter, 38 years bright-only fill | `bright-star-timing/report.md` | It scales the 26.9 million bright stars by 966 and matches nothing in this galaxy (23.9 billion sectors, 3.0 x 10^11 systems, 2.70 x 10^7 stars at 1000 L_sun). Drop or relabel it; a whole-galaxy fill is about 185 years on one worker (section 5.2). The 26.9 million column stands. |
| 500 L_sun, about 61 million stars, 88.2 billion systems, 9.6 GB table, 15 s scatter, 3.5 PB | `bright-star-preplacement-plan.md` section 2; `star-fix-spec.md` ("~88 billion systems") | 1000 L_sun since GEN.30; 26.9 million stars; table 7.4 GB (4.4 data, 3.0 index, about 164 bytes per row); scatter 24 min on 4 workers (74.6 CPU-min). GEN.118 and GEN.119 raised the galaxy to about 3.4 times as many systems, 3.0 x 10^11, so 12 PB at 40 KB per system. The plan's mix and algorithm are otherwise as built. |
| "The backfill's cost per cell is unmeasured" | [sky-view.md](sky-view.md) section 3.2 | Measured: 0.06 ms per cell for draws with a ceiling, 2 ms per sector warm with SQL, 5 s cold start (section 4.3). Its option (b), ring-by-ring angle bins, is the block-first idea; use section 2 for the keying rules. |
| `bright_star_blocks` as the backfill ledger | `docs/plan/notes.md` (GEN.44 note) | Replaced by `sector_stats` (one row per sector, schema v53: `bright_level_sol`, `level_before_fill_sol`, expected and actual counts). [reproducible-galaxies.md](reproducible-galaxies.md) does not name `bright_star_blocks` (checked 2026-10-09). OPS.25 records the dropped table. |
| The backfill is `backfill_bright_stars_around` in `generate.py` and `brightStars.backfill_cells`, "a Poisson draw for every unfilled sector of every block in range" | `docs/TODO.md` GEN.40 | Now `generation/run_galaxy.py` and `generation/bright_stars.py`; since GEN.44 it draws cell by cell in chunks of 200 sectors (`run_plan._draw_sector_bands`) with a row lock per chunk. |
| Space with "0 or nearly 0" density | `docs/TODO.md` GEN.102 | None inside the outline (minimum 0.069); use expected systems below 1 (section 5). |
| Drop sectors "so the region probability matches" | `docs/TODO.md` GEN.42, GEN.43 | Matching the mean is not enough; each cell's distribution must be kept (section 3). Superseded by the block-first draw. |

Cross-checked and standing: E = 6.31 and the outline in
[galaxy-disk-density.md](galaxy-disk-density.md); 26.9 million at 1000 L_sun in
[galaxy-coordinate-system.md](galaxy-coordinate-system.md); the 1,886-cell 100 ly
backfill in [sky-view.md](sky-view.md) (1,800 to 1,895 depending on the centre).

## Evidence notes

Computed [C] with scripts in the research scratchpad `r13`: outline and sector
count (`outline.py`); void counts and locations (`void.py`, `void2.py`); exactness
tests, benchmark and determinism (`samp.c`, `harness.py`, `bench.py`, `determ.py`,
`trunc.py`); prototype and real backfill (`proto_bf.py`, `proto_test.py`,
`realbf2.py`, `tiers_tbl.py`); the MariaDB run (`dbbf2.py`, scratch 10.11.14,
default shape, shared host); habitable-world statistics (`hab2.py`, `hab_math.py`);
centres (`centres.py`); star totals (`galaxy_totals.py`). The dispersion formula
(section 3.1) and the time columns in sections 5.2 and 9 were recomputed in the
integration pass.

Recalled [R], to verify when paper access is allowed (the research environment
could read only search-result text): Kingman 1993 and Devroye 1986; Lewis and
Shedler 1979; Salmon et al. 2011 (Philox, Threefry); Walker 1977; Vose 1991;
Bridson 2007; the Minecraft chunk status list and neighbour-radius rule (idea high
confidence, names medium); the outbox and idempotent-handler patterns; the 40 KB
per system storage figure and the 8% O/B share at 750 L_sun and above (both from
the pre-placement plan, not re-measured).

Not verified: uid claim and name confirmation inside the sector transaction
(section 7.2); database cost at scale (the run is 1,800 sectors on a shared 4-core
box); block-first in the bulge (worse there, by design); `SpawnWorker` behaviour
on Windows and macOS.

## Sources

- `rq` 2.12.0 wheel from PyPI: `rq/registry.py` `StartedJobRegistry.cleanup`,
  `rq/job.py` `retry`, `requeue`.
- Repository files read: `src/planetgen/generation/{bright_stars,run_galaxy,run_plan,star_population,system}.py`,
  `galaxy/{density,geometry,skeleton,seed,drill}.py`, `util/draw.py`,
  `queue/{redisqueue,work}.py`, `db/{models,store,control_schema,sector_paths}.py`,
  `tuning.py`; the design documents linked above.
- `/mnt/project-files/bright-star-timing/report.md` (PR #657, GEN.72, 2026-10-08) and
  `/mnt/project-files/galaxy-studies/{bright-star-preplacement-plan,star-fix-spec}.md`.
- NLTK `words` corpus from raw.githubusercontent.com/nltk/nltk_data (only so the
  project code would import in the scratch environment). No web search results
  were used.
