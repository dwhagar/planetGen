# Galaxy Disk/Spiral Density — Design Proposal (Revision 2: Preplanned Sector Table)

**Status:** superseded, not built as designed here. The density model
itself (section 1) is implemented as described
(`stellarObjects/galaxyDensity.py`). The preplanned-table idea this
revision's title refers to -- persisting every one of ~10.5 billion
qualifying sectors' position/density up front in a `galaxy_sector_plan`
table (section 6, ~1 TB per section 5's own measurement) -- was not built:
a sector's position, density, and vertices are all pure deterministic
functions of its `(shell_index, shell_slot_index)` address and this
model's own shape parameters, so none of it needs persisting per sector at
all. What was actually built instead -- a singleton shape row plus one row
per shell recording a *candidate* slot-index band, ~4,000 rows total for a
real Milky-Way-scale galaxy -- is documented in
`docs/design/galaxy-coordinate-system.md` section 9's storage-analysis
addendum, `stellarObjects/galaxySkeleton.py`, and `galaxyPlan.py`. Kept
here for the density model itself (still accurate) and as a record of the
storage-cost reasoning that led to the smaller design actually shipped;
section 5's feasibility numbers and section 6's schema no longer describe
what's in the database.

**Scope:** unchanged from revision 1 (a density model over galaxy-frame
position; how many sectors that implies and where; how `--density` gets
scaled at actual generation time). What's new in this revision: the density
model is evaluated and stored for *every* address that clears a real
threshold, in a new table, before any actual sector content exists — not
computed ad hoc per address at generation time.

## 0. Decisions carried over from the user

- **Preplan every sector, stored in its own table.** Not computed on the fly
  during `galaxyGen.py --shell`/`--center-sector` runs — a new table holds
  one row per qualifying address, built by a dedicated batch pass before any
  actual content generation happens.
- **Gate at "predicted less than 1 star possible," evaluated per sector,
  stopping outward once a radial shell no longer clears it.** This replaces
  revision 1's probabilistic occupancy gate with a hard deterministic
  threshold: a sector's own `expected_system_count() * relative_density`
  must be `>= 1` to exist in the plan at all (§2). Once an entire shell
  cannot clear this bar anywhere in it, radial generation stops (§3) — this
  *is* the galaxy's edge, derived from the model rather than picked as an
  arbitrary radius.
- **Every planned sector gets a sequential index** (§4).
- **Assume a spiral galaxy similar to the Milky Way** — real-astronomy-scale
  parameters (§1), not arbitrary ones.

## 1. Density model (unchanged shape, tuned parameters)

Same exponential-disk-plus-bulge-plus-spiral-arm model as revision 1 (see
that revision's §1 for the full derivation; restated concretely here with
final parameter choices and the normalization made explicit, which revision
1 had left underspecified):

```
r_cyl = sqrt(x^2 + y^2)
r_3d  = sqrt(x^2 + y^2 + z^2)          # galaxyGeometry.galactic_radius_pc
theta = atan2(y, x)

rho_bulge(r_3d)          = bulge_amplitude * exp(-r_3d / bulge_scale_radius_pc)
rho_disk_radial(r_cyl)   = exp(-r_cyl / disk_scale_length_pc)
f_z(z)                   = 1 / cosh(z / disk_scale_height_pc)^2
theta_arm(r_cyl)         = spiral_reference_angle_rad
                           + ln(r_cyl / disk_scale_length_pc) / tan(pitch_angle_rad)
arm_factor(r_cyl, theta) = 1 + arm_amplitude * cos(arm_count * (theta - theta_arm(r_cyl)))

raw_density(x, y, z) = rho_bulge(r_3d)
                        + rho_disk_radial(r_cyl) * f_z(z) * arm_factor(r_cyl, theta)

relative_density(x, y, z) = K_NORM * raw_density(x, y, z)
```

### Normalization, made explicit

Revision 1 said parameters were "chosen so `relative_density = 1.0`" at
Sol's position without naming the actual multiplier that requires. There is
one: `K_NORM = 1 / raw_density(R_sun_pc, 0, theta_interarm)`, evaluated at
Sol's real galactocentric radius (`R_sun_pc = utils.ly_to_pc(
physical_constants.GALACTIC_CENTER_DISTANCE_LY)`, ~7,910 pc) and at the
**inter-arm minimum** azimuth at that radius (`theta` where `arm_factor`
hits its floor `1 - arm_amplitude`) — chosen deliberately, not arbitrarily:
it's the popular-astronomy fact that the Sun sits between two spiral arms
(the Local/Orion Spur, not on Perseus or Sagittarius), and it's the
conservative calibration choice (if `relative_density = 1.0` on-arm
instead, everything off-arm at the same radius would read as
sub-realistic, which is backwards for "how dense is a typical patch of the
solar neighborhood"). `K_NORM` is **computed at runtime from the other
constants**, never hardcoded — it must be re-derived automatically if any
of the shape parameters below ever change.

### Parameters (Milky-Way-scale, tuned during this design pass)

| Constant | Value | Basis |
|---|---|---|
| `disk_scale_length_pc` | 2,800 | `R_sun_pc / disk_scale_length_pc ~= 2.8`, matching the real Milky Way's ratio |
| `disk_scale_height_pc` | 350 | Between real thin-disk (~300 pc) and thick-disk (~900-1,000 pc) |
| `bulge_scale_radius_pc` | 200 | See note below — smaller than revision 1's first guess (500 pc) |
| `bulge_amplitude` | 1.0 | See note below — smaller than revision 1's first guess (5.0) |
| `arm_count` | 2 | Grand-design two-arm spiral |
| `pitch_angle_deg` | 15 | Typical grand-design spiral pitch (real spirals run ~10-25 deg) |
| `arm_amplitude` | 0.4 | Arm/inter-arm contrast `(1+0.4)/(1-0.4) = 2.33x` |
| `spiral_reference_radius_pc` | = `disk_scale_length_pc` | Arbitrary but fixed anchor for `theta_arm` |
| `spiral_reference_angle_rad` | 0.0 | Arbitrary fixed orientation |

**Why the bulge parameters shrank from revision 1's first guess**: this
design pass's own feasibility investigation (§5) found that
`bulge_amplitude = 5`, `bulge_scale_radius_pc = 500` makes the bulge term
*alone* (independent of disk position, since it only depends on `r_3d`)
exceed the 1-star-per-sector threshold out to shell ~836 (`r ~= 2,950 pc`)
— meaning *every* direction in *every* shell inside that radius trivially
qualifies, not just a thin disk plane. That's a real "the galaxy has a
genuinely 3D bright core, not just a flat disk" fact, not a bug — but it
also means no phi-band pruning is possible anywhere inside that radius
(every one of those shells' full `N_k` slots must be visited, since every
one of them is a real qualifying sector needing its own stored density
value). At `bulge_amplitude = 1.0`, `bulge_scale_radius_pc = 200`, that
"trivially-solid" region shrinks to shell ~241 (`r ~= 850 pc`, cumulative
~60 million slots) — still a real, physically sensible bright bulge core,
just sized so the rest of the galaxy's pruning (§3) actually pays off. Both
values are still tunable by eye once implemented; this is the pair that
made the batch build tractable, not a rigorously derived astrophysical fit.

## 2. The occupancy threshold — deterministic, not probabilistic

Revision 1 proposed a hash-based Bernoulli roll to decide occupancy,
because it was treating "gate" as a soft probability. This revision's
brief is a hard, physical criterion instead, which is simpler and needs no
randomness at all:

```
E = SpaceSector(edge_ly=DEFAULT_SECTOR_EDGE_LY).expected_system_count()   # at density=1
    # = DEFAULT_SECTOR_EDGE_LY^3 * LOCAL_STELLAR_DENSITY_LY3 ~= 4.32 systems/sector

predicted_star_count(position) = E * relative_density(position)

qualifies(position)  <=>  predicted_star_count(position) >= 1.0
```

A sector's fate is a fixed fact of its position, computed once, exactly —
not a coin flip. This is strictly simpler than revision 1's design and
directly matches the brief ("gate density at a predicted less than 1 star
possible per sector").

## 3. Finding the edge, shell by shell, without visiting every slot

Brute-force enumeration (evaluate every one of a shell's `N_k` slots to see
which qualify) is correct but far too slow at outer shells — see §5's
numbers. Two facts, both already true of this model, make targeted pruning
possible:

1. **The bulge term is constant across an entire shell** (it depends only
   on `r_3d`, which every slot in shell `k` shares: `r_3d = r_k`). So
   `bulge_alone(k) = K_NORM * bulge_amplitude * exp(-r_k / bulge_scale_radius_pc)`
   is one number per shell. If `bulge_alone(k) >= THRESHOLD_RHO`
   (`= 1/E`), **the entire shell qualifies, unconditionally** — no pruning
   needed or possible there; every slot's own exact density must still be
   computed and stored (for the weighting value, §6), so this is a real
   `O(N_k)` cost for shells this deep in the core, not a shortcut.
2. **Past that regime, the true per-`theta` maximum at a given `(k, phi)`
   is exact and cheap**: `arm_factor`'s own maximum over `theta` is exactly
   `1 + arm_amplitude`, independent of `theta` itself, so fixing it there
   gives an exact (not merely conservative) upper bound on density at that
   `phi`:

   ```
   bound(k, phi) = K_NORM * (
       bulge_alone(k)
       + (1 + arm_amplitude) * exp(-r_k*sin(phi)/disk_scale_length_pc)
                              * f_z(r_k*cos(phi))
   )
   ```

   For each shell past the solid-core regime, sample `bound(k, phi)` on a
   fixed grid across `phi in [0, pi]` (500 points; this design pass used
   200 for speed and found it adequate, but recommends 500 for the real
   implementation as cheap extra safety margin — the grid is evaluated once
   per shell, ~4,000 shells total, trivial either way) and keep the
   contiguous run of grid points (plus one grid-step of padding on each
   side, the same safety-margin pattern `galaxyGeometry._candidate_shell_range`
   already uses) where it clears `THRESHOLD_RHO`. Feed that `[phi_min,
   phi_max]` into the *already-existing*
   `galaxyGeometry._slot_index_bounds_for_phi_range` (built for the
   neighborhood-query primitive, reused as-is here) to get the actual slot
   index range to visit. Every visited slot still gets its *exact* position
   and density computed and checked — the grid band is a safe (possibly
   slightly wider than necessary) net, never a final answer on its own.

   **Known imprecision, accepted deliberately**: this is a grid scan, not a
   closed-form inversion (this design pass found the true qualifying region
   isn't always a single band centered on the disk plane — for some
   mid-range shells, a *second* small qualifying patch exists near the
   poles too, since a slot there trades a favorable `r_cyl -> 0` against an
   unfavorable large `|z|`, and depending on the shell radius either term
   can win). A fixed grid handles this correctly (it doesn't assume
   unimodality) but could in principle miss an extremely narrow qualifying
   sliver that falls entirely between two grid points, right at the
   galaxy's outer edge where the qualifying band is thinnest. Accepted as a
   bounded, cosmetic imprecision (a handful of borderline sectors at the
   very edge, not a systemic error) rather than solved with exact root
   isolation, which would need calculus this model's shape doesn't obviously
   guarantee is well-behaved everywhere.

**Stopping the outward scan**: once a shell's grid scan finds no point
anywhere that clears the bar, both density terms are past their peak and
monotonically falling with radius, so every shell beyond it is empty too —
but as insurance against the imprecision just noted, the real
implementation should confirm a **run of consecutive empty shells** (e.g.
50) before concluding the edge has been reached, rather than stopping at
the very first one.

## 4. Sequential sector index

Every row inserted into the plan table (§6) gets `sector_index`, a dense
`1..N` integer assigned in the same canonical order the build walks: shell
`0, 1, 2, ...` outward, and within each shell in ascending
`shell_slot_index` order. This makes `sector_index` a stable "the Nth
star-bearing sector outward from the core" address — meaningful on its own
(a lower index is always closer to the galactic center) and stable across
a resumed/re-run build (§7), since it's assigned in a fixed traversal
order, never renumbered.

## 5. Feasibility investigation — real measurements, not estimates

Run as part of this design pass (see the session's own benchmarking, not
reproduced as scripts here, but summarized because the numbers drove every
decision above):

- **Full content generation** (`sectorGen.py`, real systems/planets, 1,000
  sectors measured): **22.9 ms/sector, 87.8 KB/sector** in the database.
  Confirms generating actual content for the whole plan is never on the
  table (see below) — expected, and fine, since a galaxy is always mostly
  unvisited.
- **Total addressable slots**, shells 0 through the computed edge (~4,100
  with revision-1-era parameters): **~289 billion.**
- **Qualifying sectors** (properly weighted by true per-shell slot counts,
  not naive uniform shell sampling, which badly overestimates the
  fraction): **~3.6% -> ~10.5 billion** sectors predicted to hold at least
  one star.
- **Brute-force compute** (evaluate every one of the 289 billion
  addresses): ~293 hours single-core. Not viable.
- **With §3's pruning**: total addresses actually needing evaluation drops
  to an *exact* (not sampled) **12.62 billion** — a ~23x reduction — and
  the pruning decision itself (which shells/bands to visit) costs under a
  second for the whole galaxy.
- **With NumPy vectorization** of the density formula: **241 ns/candidate**
  vs. ~3,000 ns in a pure-Python loop (~12x). Combined with pruning:
  **12.62B x 241ns ~= 51 minutes, single core** — and this step is
  embarrassingly parallel (every candidate is independent), so it drops
  further with more cores.
- **SQLite bulk-insert tuning** (`journal_mode=OFF`, `synchronous=OFF`
  during the one-time build, index created *after* loading instead of
  during): **3.6 us/row vs. 5.3 us/row** (~1.5x). This step does *not*
  parallelize as easily — SQLite allows one writer at a time — so
  inserting ~10.5 billion rows this way is genuinely the long pole: **~6.5-
  10.5 hours, serialized**, on top of the sub-hour compute step.
- **Storage**: ~80 bytes/row measured -> **~830 GB-1 TB** for the full
  table (before accounting for the post-build index's own size). This
  number is **not reducible by any of the above** — it's a direct
  consequence of a 3.5-pc sector next to a ~13,500 pc galaxy radius, not an
  algorithm inefficiency.

**Bottom line, and the basis for building the full table anyway (§0)**:
what looked like a 12+ day, clearly-infeasible job at first estimate is a
several-hour-to-overnight batch job once pruned and vectorized, requiring
~1 TB of free disk. Compute is no longer the constraint; one long,
necessarily-serialized write pass and the disk budget are.

## 6. Schema: new `galaxy_sector_plan` table

```sql
CREATE TABLE IF NOT EXISTS galaxy_sector_plan (
    sector_index              INTEGER PRIMARY KEY,   -- see section 4
    shell_index               INTEGER NOT NULL,
    shell_slot_index          INTEGER NOT NULL,
    center_x_pc               REAL NOT NULL,
    center_y_pc               REAL NOT NULL,
    center_z_pc               REAL NOT NULL,
    galactic_radius_pc        REAL NOT NULL,
    predicted_relative_density REAL NOT NULL,
    predicted_star_count      REAL NOT NULL,          -- E * predicted_relative_density, cached
                                                       -- to avoid re-deriving E at query time
    generated_sector_id       INTEGER REFERENCES sectors(id),  -- NULL until galaxyGen.py
                                                                -- actually generates this one
    UNIQUE (shell_index, shell_slot_index)
);
-- Created only after the bulk build finishes (see section 7) -- incremental
-- index maintenance during a multi-billion-row insert is the ~1.5x this
-- design pass measured leaving on the table otherwise.
CREATE INDEX IF NOT EXISTS idx_galaxy_sector_plan_shell ON galaxy_sector_plan(shell_index);
CREATE INDEX IF NOT EXISTS idx_galaxy_sector_plan_ungenerated
    ON galaxy_sector_plan(shell_index) WHERE generated_sector_id IS NULL;
```

`PRAGMA user_version` moves to the next version past whatever Track C left
it at, per `docs/database-schema.md`'s existing history-entry convention.

**A sector's own `center_x/y/z_pc`/`galactic_radius_pc` columns (already on
`sectors` since Track C) are left as-is, populated the same way they are
today once a plan row is actually generated** — genuinely redundant with
the plan row's own copy of the same numbers, but changing `sectors`' shape
or how existing code (the web UI, `queryDb.py`) reads it is out of scope
here; duplication is the accepted cost of not touching that surface.

## 7. `galaxyPlan.py` — the batch build tool

New root-level script, alongside `sectorGen.py`/`systemGen.py`/
`galaxyGen.py`, whose only job is building (or resuming) the plan table —
it never generates sector *content*, only the census of where content
*could* go:

- Creates `galaxy_sector_plan` (and its schema-version bump) if missing.
- **Resumable by construction**: on start, reads
  `MAX(shell_index)` already present and continues the outward shell scan
  from there — an interrupted multi-hour run loses at most the
  in-progress shell, not prior work. `sector_index` continuation follows
  from `MAX(sector_index)` the same way, preserving section 4's ordering
  guarantee across a resume.
- Per shell: the solid-core fast path or the pruned-band scan (§3),
  computed via NumPy in one vectorized batch per shell rather than a
  Python-level loop per slot (§5's 12x).
- **Write mode**: `PRAGMA journal_mode=WAL` + `synchronous=NORMAL` as the
  default (crash-safe, still fast) with commits every shell (bounding
  how much a crash mid-shell can lose); an opt-in `--unsafe-fast` flag
  drops to `journal_mode=OFF`/`synchronous=OFF` for a known-safe
  environment willing to trade crash resilience for the last ~1.5x, per
  §5's measurement. Indexes are created once, after the scan finishes
  entirely (or is deliberately stopped early via `--max-shell`, for
  testing on a laptop-sized slice before committing to the full run).
- Progress output every shell (or every N seconds): current shell index,
  cumulative sector count, elapsed time — this is a multi-hour job that
  should never run silently.

## 8. `galaxyGen.py` changes: consume the plan instead of computing occupancy

`--shell K` and `--center-sector ID --radius-pc R` keep their existing
meaning, but stop calling `galaxyGeometry.sector_position_pc`/
`shell_sector_count` directly to decide *what exists* — that question is
already answered, once, by the plan table:

```sql
-- --shell K, not-yet-generated slots:
SELECT * FROM galaxy_sector_plan
WHERE shell_index = :k AND generated_sector_id IS NULL;
```

For `--center-sector`, `galaxyGeometry.enumerate_sectors_within_radius`
still finds the *candidate* `(shell_index, shell_slot_index)` addresses
near the center sector (unchanged, still needed since that primitive
searches by geometry, not by plan membership) — the result is then
filtered against `galaxy_sector_plan` by `(shell_index, shell_slot_index)`;
a candidate with no matching plan row simply isn't a real sector (its
predicted density never cleared the threshold) and is skipped, same as
today's "already occupied" skip, just checking plan membership instead of
`sectors` occupancy.

For each plan row actually generated: `sectorGen.generate_sector(args,
galactic_center_dist_ly=..., density_multiplier=plan_row.predicted_relative_density)`
(the same hook revision 1 proposed, unchanged — §1's density value is what
flows into it), then `UPDATE galaxy_sector_plan SET generated_sector_id =
:id WHERE sector_index = :sector_index`.

No hash-based occupancy roll remains anywhere in this design — revision
1's `_occupancy_roll` is dropped entirely, superseded by §2's deterministic
threshold decided once at plan-build time.

## 9. Open questions carried forward

1. **Sharded/parallel plan build.** §5 measured the insert step as the
   long pole specifically because SQLite allows one writer per file, not
   because parallel *computation* isn't possible. Building into `N`
   separate files (one per worker, sharded by shell range) in parallel and
   merging them afterward (via `ATTACH` + bulk `INSERT ... SELECT`, likely
   far faster than the row-by-row Python API this design pass benchmarked)
   could cut the ~6.5-10.5 hour insert step significantly, at real added
   implementation complexity this design pass didn't measure. Left for a
   follow-up if the single-writer build proves too slow in practice.
2. **Exact final qualifying count and edge shell** depend on the precise
   `K_NORM`/bulge/disk constants finally chosen (§1's table is this pass's
   recommendation, not yet locked in) — §5's ~10.5 billion and ~4,100-shell
   figures came from evaluating this design's own formulas during the
   investigation and will shift somewhat (not by orders of magnitude) once
   the real implementation calibrates `K_NORM` at the precise inter-arm
   point specified in §1 rather than the placeholder calibration this
   pass's benchmarks used.
3. Everything revision 1's own §6 open questions listed (thin/thick disk
   split, metallicity-vs-radius gradients) still applies unchanged.

## 10. Testing plan

- `src/tests/test_galaxy_density.py`: `relative_density`/`raw_density`
  hand-checked at exact points (bulge center, solar calibration point
  evaluating to exactly `1.0`, on-arm vs. inter-arm contrast, in-plane vs.
  off-plane falloff), `K_NORM`'s derivation, `expected_star_count`
  threshold arithmetic.
- `src/tests/test_galaxy_plan_geometry.py`: the solid-core fast path and
  grid-band scan against small, hand-verifiable shells (a shell forced
  fully solid; a shell with a known, checkable band; a shell confirmed
  empty) — plus a slow/regression test comparing pruned output against
  brute-force full-shell enumeration on a handful of small shells, to
  guard the pruning's correctness (never excluding a true qualifier).
- `src/tests/test_galaxy_plan.py`: `galaxyPlan.py` end-to-end on a small
  `--max-shell` slice — resumability (kill mid-build, restart, confirm no
  gaps/duplicates in `sector_index`), schema creation, index-after-build
  ordering.
- `galaxyGen.py`'s existing `test_galaxy_gen.py` gains coverage for
  plan-table consumption (a shell with some plan rows already generated,
  some not; a `--center-sector` search whose candidates only partially
  exist in the plan).
