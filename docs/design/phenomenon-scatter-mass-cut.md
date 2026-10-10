# Phenomenon scatter: a lowest-mass cut, filled below at sector time

Boss (2026-10-09 19:54Z): "we need to have the phenomena scatter use a lower rate, this is
about mass, so recommend a lowest-mass-limit to still provide good gravitational data across
the unfilled galaxy but does not overload the system with creation, then when we make a
sector we'll fill it with the proper number below that mark, similar to how we do bright
stars." This note picks the cut with numbers, says how the sector fill draws what lies below
it, and what that does to reproducibility and to DB.19. It extends
[generation-performance-study.md](generation-performance-study.md) (which found 1.06e9
phenomenon rows) and uses the influence rule of [orbital-updates.md](orbital-updates.md) and the
potential of [galactic-potential.md](galactic-potential.md).

Informs: DB.19, GEN.100, GEN.104, GEN.109, GEN.115, PERF.31

Status: built, 2026-10-09 (GEN.166 to GEN.168, schema v71); extended 2026-10-10 by the star passes of GEN.185 (below); Boss accepted the 20 solar mass cut on a decision card (2026-10-09).

Update (Boss, 2026-10-10 03:01Z): the cut becomes a user setting with preset values between 8 and 20 solar masses (GEN.183; 20 stays the default), the bright-star luminosity floor becomes a preset from 2500 to 4,000,000 L_sun with a default of 3000 (GEN.184), and the scatter runs in five passes (GEN.185). The fixed 20 solar mass cut and the 1000 L_sun default below are superseded as defaults.

As built: the cut is `galaxy_shape.phenomenon_min_mass_solar` and `--phenomenon-min-mass` (in the settings file's plan settings); the scatter's classes are `phenomenon_scatter.SCATTER_CLASSES`, and a sector's below-cut draw is `phenomenon_scatter.below_cut_draws` with `run_sector.add_below_cut_remnants`. No per-sector level or band top-up was needed: a filled sector already holds every mass (its scattered rows above the cut, its own draw below), and a rescatter at a new cut (`planetgen plan --phenomena-only`) replaces the rows of every unfilled cell and leaves filled ones out, so nothing is missing or doubled whichever way the cut moves.

Evidence tags: [C] computed or measured here (`/mnt/project-files/research/scripts/genperf/`:
`default_totals.py`, `mass_cut_table.py`, `walk_cost.py`); [S] read in the repository; [R] recalled and unconfirmed.

## The answer in short

1. **Recommended cut: 20 solar masses.** The scatter then keeps the intermediate-mass black
   holes (100 to 100,000 solar masses) and the other kinds that are not drawn by mass (planetary
   nebulae, supernova remnants, hypervelocity stars, the nucleus): about **2.7e5 rows instead of
   1.17e9**, 37 MB instead of 161 GB [C]. Every neutron star (1.1 to 2.2) and every stellar-mass
   black hole (5 to 20) is drawn when its sector is made.
2. **Why 20.** An object's influence reaches its Hill (tidal) sphere, `r_t = 1.44 pc (m / Msun)^(1/3)`
   [C]; that equals one 4 pc sector edge at 21 solar masses. Below that, an object's pull ends within
   about one sector of where it sits, and the sector fill (with its neighbours) covers that. Above
   it, an object matters to sectors that may never be filled, so it needs a row. 20 is also the
   model's stellar-black-hole maximum, so the cut falls in the gap between "stellar" and "intermediate".
3. **The error is small for every cut, because it is flat.** The whole range of cuts, from keeping
   everything to keeping only the intermediate-mass holes, moves the missing share of the local mass
   density from 0 to 1.8%, and the share of objects with a missing remnant inside their influence
   sphere from 0 to 2.2% [C]. Cutting at 20 costs 2.2 points for 4,000 times fewer rows; no cut in
   between buys a point of error for fewer than 2e8 rows.
4. **The galactic field is untouched.** The field is the analytic Hernquist, Miyamoto-Nagai and NFW
   potential, which already holds the local mass density of 0.10 Msun/pc3 [S]; the scattered remnants
   are not in it and were never needed for it.
5. **Where the 1.6e8 came from:** an unverified estimate, not a result. See "Correcting the counts".

## Correcting the counts

Three figures disagreed; this fixes them.

- **The 1.6e8 rows** appear in `docs/plan/notes.md` (a line of 2026-10-09 08:15Z, written when GEN.100
  shipped) and in `db-check-and-parity-repair.md`, whose own text calls it "unverified" and derives its
  21 GB from "a 5.0 million row copy" scaled by 32 [S]. Nothing in the code derives it. The default
  galaxy's nominal counts are 1.18e9 neutron stars and **1.48e8 black holes** [C], so the figure matches
  the black holes alone, not the total.
- **The default galaxy** has 2.953e11 systems, 2,041 layers and an outer ring of 3,763 [C]. At the
  tuned rates (`tuning.phenomenon_rate_per_star`: neutron stars 4e-3, black holes 5e-4, planetary
  nebulae 1e-7, supernova remnants 7.1e-8) the nominal rows are 1.18e9 neutron stars, 1.48e8 black
  holes, 2.95e4 planetary nebulae and 2.1e4 supernova remnants [C].
- **What the scatter actually draws** is 64 times the quarter-scale run (18,241,331 rows, so about
  **1.17e9 rows and 161 GB at 138 B a row**), not the 1.06e9 and 146 GB that
  `phenomenon_scatter.layer_expected` gave the earlier note; that estimate is the progress bar's
  weight and runs 10% low [C]. The split is 9.5e8 neutron stars and 2.2e8 black holes.

**One thing for Boss.** Against the nominal rates the neutron stars come out 0.80 times and the black
holes 1.49 times [C], because the regional factors of GEN.132 (`remnant_distribution.placement_factor`,
vertical scale-height ratios 3.0 and 2.5, the black holes' core excess) multiply the bin weights without
renormalising. His 2026-10-09 retune said "1e9 neutron stars and 1e8 black holes"; the scatter yields
9.5e8 and 2.2e8, so the black holes are 2.2 times his stated total. If "the rate is the total" is the
intent, divide each kind's weights by their mean factor. The cut below does not depend on it, except that
the intermediate-mass black holes (0.1% of the black holes) would drop from 2.2e5 to 1e5.

## What mass the scatter has, and what it is for

Masses are drawn when a row is built into an object, from the row's seed [S]:
neutron stars uniform 1.1 to 2.2 (mean 1.65); stellar-mass black holes uniform 5 to 20 (mean 12.5);
with a 0.1% chance a black hole is intermediate-mass, log-uniform 100 to 1e5 (mean 14,500) [S]. So 99.9%
of the black holes are below 20 and the 0.1% above it carry **53% of all black-hole mass** (2.2e5 holes
of 14,500 mean: 3.2e9 solar masses against 2.7e9 in the stellar ones) [C]. If that share is not what Boss
meant by 0.1%, it is a tuning question (BLACK_HOLE_INTERMEDIATE_MASS_CHANCE), not a cut question.

Two jobs for a pre-placed mass:

1. **The galactic field.** The model is analytic and already contains the local density of
   0.1008 Msun/pc3 [S]. A neutron star or black hole is a former star that the density already counts.
   Local mass in the remnants is 9.2e-4 (neutron stars), 8.7e-4 (stellar black holes) and 1.0e-3
   (intermediate-mass) Msun/pc3 [C], 0.9%, 0.9% and 1.0% of the field. Dropping rows changes the
   field by zero, because the field is not summed from rows.
2. **The influence set.** For each object visited, the influence radius is the Hill sphere of the largest
   nearby object, and the point masses inside it are summed with the galactic gradient
   ([orbital-updates.md](orbital-updates.md) section 1) [S]. A missing remnant matters when it lies inside
   that sphere. That is the only place a row's mass changes a result.

## The tidal sphere makes the choice

For a body of mass m in the galaxy's flat rotation curve, `r_t = (G m / (4 Omega^2 - kappa^2))^(1/3)`
with Omega = 220 km/s over 8.2 kpc and kappa^2 = 2 Omega^2 gives **1.44 pc for one solar mass**, which
matches the 1.4 pc the orbital-updates note quotes [C][S]. The sphere's volume is proportional to the
mass (`12.5 pc3` per solar mass). Two consequences:

- The chance that a visited object lies inside the tidal sphere of some remnant is
  `12.5 pc3/Msun times the local mass density of those remnants`. This is the error measure below.
- A remnant's sphere is wider than a sector (4 pc) from `(4 / 1.44)^3 = 21` solar masses. Under that
  mass its influence stays within the sector it is in plus the ring of sectors around it.

## The cuts

At the default galaxy (counts: quarter-scale scatter times 64; mass laws from `tuning`) [C]. "In
reach" is the share of visited objects whose influence sphere holds a missing remnant, valid where all
the surrounding sectors are unfilled; inside a filled block it is zero.

| Cut (Msun) | Neutron-star rows | Black-hole rows | All rows | GB at 138 B | GB at 45 B | NS mass kept | BH mass kept | Missing mass, share of local density | Objects with a missing remnant in reach | Heaviest missing, tidal radius |
|---|---|---|---|---|---|---|---|---|---|---|
| none (today) | 9.47e8 | 2.20e8 | 1.17e9 | 161 | 53 | 100% | 100% | 0% | 0% | none |
| 1.65 | 4.74e8 | 2.20e8 | 6.93e8 | 96 | 31 | 58% | 100% | 0.38% | 0.5% | 1.7 pc |
| 2.0 | 1.72e8 | 2.20e8 | 3.92e8 | 54 | 18 | 23% | 100% | 0.70% | 0.9% | 1.8 pc |
| 2.2 to 5 | 0 | 2.20e8 | 2.20e8 | 30 | 10 | 0% | 100% | 0.92% | 1.2% | 1.9 to 2.5 pc |
| 10 | 0 | 1.46e8 | 1.46e8 | 20 | 6.6 | 0% | 91% | 1.09% | 1.4% | 3.1 pc |
| 15 | 0 | 7.3e7 | 7.3e7 | 10 | 3.3 | 0% | 75% | 1.38% | 1.7% | 3.6 pc |
| **20** | 0 | **2.2e5** | **2.7e5** | **0.04** | **0.01** | 0% | 54% | **1.78%** | **2.2%** | 3.9 pc |
| 100 | 0 | 2.2e5 | 2.7e5 | 0.04 | 0.01 | 0% | 54% | 1.78% | 2.2% | 3.9 pc |
| 1,000 | 0 | 1.5e5 | 2.0e5 | 0.03 | 0.01 | 0% | 53% | 1.79% | 2.3% | 14 pc |

The "other" rows (planetary nebulae, supernova remnants, hypervelocity stars, the nucleus) are 5.0e4
and are kept at every cut: they are few, and they are not placed by mass.

Reading it:

- **There is no knee until 20.** Each solar mass of cut between 5 and 20 removes 1.5e7 rows and adds
  0.06 points of error; the neutron stars (9.5e8 rows) are worth 1.2 points as a whole. The
  intermediate-mass holes are worth 1.3 points for 2.2e5 rows, so they are the only kind where a row
  buys much, and they are the kind with tidal spheres far wider than a sector (14 pc at 1,000, 67 pc at
  1e5).
- **A budget view:** the bright-star scatter writes 2.6e7 rows. A cut that spent as many rows on
  remnants (about 18 Msun) buys 0.13 points over 20. 20 is the cut.
- **If Boss wants every stellar black hole placed** (cut 5, 2.2e8 rows, 10 GB compact), the error falls
  to 1.2% and the scatter costs about 2 CPU-hours more; I do not recommend it, but the knob supports it.

## What the sector fill does

Same pattern as the bright-star band ([run_plan.py](../../src/planetgen/generation/run_plan.py)
`_draw_sector_bands`, `bright_share`): the plan stores the cut, every sector remembers the floor it was
filled to, and the fill draws exactly what the scatter left out.

1. **Count.** A kind's objects in a sector are Poisson with mean `rate_k × expected stars × regional
   factor`, split by Poisson thinning into independent counts above and below the cut. The scatter draws
   `lambda × S_k(c)` per ring bin, where `S_k(c)` is the share of the kind's mass law above the cut
   (neutron stars `(2.2 - c)/1.1`, stellar holes `(20 - c)/15`, intermediate `ln(1e5/c)/ln(1e3)`, clamped
   to 0 to 1). The sector draws `lambda × (1 - S_k(c))`, with `lambda` taken from the same expected-star
   count the scatter used (the ring-bin density), not from the sector's rolled star count, so the totals
   match the plan's statistics rather than adding the sector's own Poisson variance on top. Treat the two
   black-hole classes (stellar, intermediate) as separate kinds with rates `0.999 r` and `0.001 r`.
2. **Mass.** Inverse-CDF from the same law, truncated: above the cut `m = max(c, lo) + u (hi - max(c, lo))`
   (log-uniform for the intermediate class); below it `m = lo + u (min(c, hi) - lo)`. The constructors of
   `NeutronStar` and `BlackHole` need a `mass_range` argument in place of their fixed
   `draw.uniform(*tuning...)` ([compact_remnant.py](../../src/planetgen/generation/phenomena/compact_remnant.py)
   lines 255 to 265, 464) so that a stored row's seed rebuilds the same mass it was drawn with.
3. **Place.** Uniformly in the sector's cell, as a rogue planet is, with a sector-level regional factor
   (sectors are 4 pc; the scatter's bin weights are for the galaxy-wide walk).
4. **Level.** The sector records the cut it was filled under (a column beside the bright-star level in
   `sector_stats`), so a later lower cut is a band: the scatter draws `[c_new, c_old)` into unfilled cells
   and the top-up draws the same band into filled sectors, as `_needs_band` does for stars.

Cost: the same number of neutron stars and black holes is built across the galaxy as today; the sector
draws them instead of reading rows, so the fill gains no work and loses one SELECT. A 205-system sector
holds 0.8 neutron stars and 0.1 black holes on average [C].

## Reproducibility

- Above-cut rows keep per-layer streams (`{seed}:phenomena:{layer}`), so layers can still be drawn in any
  order or in workers. Truncating the means changes every downstream draw in those streams once: the same
  seed gives a different galaxy than before the change (no backward compatibility, Boss 2026-10-07).
- Below-cut objects come from a per-sector stream (`{seed}:phenomena-fill:{ring}:{layer}:{slot}`,
  SHA-256-derived like the other per-unit seeds, `reproducible-galaxies.md`), so a sector's remnants are the
  same whichever sectors are filled and in whatever order, and do not depend on whether a neighbour is filled.
- The cut is part of the plan: store it in the scatter settings and the settings file, and include it in the
  inputs of the reproducibility key, as the bright-star threshold is. Same seed and same cut give the same
  galaxy; a different cut is a different galaxy until the band is drawn.
- Rebuilding a row needs the cut that drew it (the mass range is part of the object's draw), so a changed
  cut means a rescatter of the layer, not a patch of rows.

## What this does to DB.19 and to the scatter

- **DB.19's two options become unnecessary at a cut of 20.** 2.7e5 rows of 138 bytes are 37 MB. Neither
  deriving neutron stars and black holes on demand nor a compact row is needed; keep the compact row as the
  fallback if Boss lowers the cut to 10 or below (1.5e8 rows, 6.6 GB compact).
- **The scatter's remaining cost is the ring walk.** With the rates scaled to nothing the quarter-scale
  layer 0 still takes 0.24 s (941 rings, 0.26 ms a ring) [C]. The default galaxy has 3,136,126 rings, so
  about 815 CPU-s (14 CPU-minutes, 4 minutes at 4 workers) against 10 to 12 CPU-hours of rows before (1.17e9 at 32 to 37 us); that walk is the
  next thing to cache if it matters (`population_densities` per ring bin).
- **The phenomenon scatter stops dominating `planetgen plan`:** about 14 CPU-minutes and 2.7e5 inserts
  against the bright-star scatter's 39 CPU-minutes warm.
- **What is lost.** Neutron stars and stellar-mass black holes of unfilled space no longer exist as rows,
  so a search for "nearest neutron star" or a map layer of remnants shows only the heavy ones, plus whatever
  sectors have been filled, exactly as ordinary stars below the bright-star threshold behave. The density
  tiles of the Galaxy Map come from the analytic model and are unaffected.

## Proposed items (the TODO thread assigns IDs)

1. A plan parameter, `--phenomenon-min-mass` (default 20 solar masses), stored with the scatter settings and
   in the settings file and the reproducibility key; the scatter draws only objects above it (GEN.100
   follow-up).
2. A truncated mass draw and `mass_range` argument on `NeutronStar` and `BlackHole` (and the intermediate
   class as its own kind), so a row's seed rebuilds its mass.
3. The sector fill's below-cut draw: per-sector stream, Poisson thinning from the scatter's expected-star
   count, level recorded per sector, band top-up when the cut is lowered (the bright-star pattern).
4. DB.19: reword to "implement the cut; derive-on-demand and compact rows only if the cut goes to 10 or below."
5. Decide for Boss: renormalise the regional factors so each kind's total equals rate times stars (today
   neutron stars 0.80 and black holes 1.49 of nominal), and whether 0.1% intermediate-mass black holes
   (53% of black-hole mass) is intended.
6. Correct `docs/plan/notes.md` and `db-check-and-parity-repair.md`: replace 1.6e8 by the computed counts
   (1.17e9 before the cut, 2.7e5 after) so DB.15 and the parity repair size their work for the right table.

## Evidence notes

`default_totals.py` builds the plan's skeleton in memory (defaults: 2600 pc scale length, 300 pc height,
1580 pc bulge radius) at scale 0.25 and 1.0 and sums the expected systems, bright stars and nominal
phenomenon rows. `mass_cut_table.py` prints the table from `tuning`'s mass laws and the quarter-scale
scatter counts. `walk_cost.py` times the ring walk with the rates at 1e-6. The scatter counts are from the
quarter-scale run of the performance study (14.8M neutron stars, 3.43M black holes, 445 planetary nebulae,
306 supernova remnants, 1,765 hypervelocity stars). The local density of 0.1008 Msun/pc3 and the flat
220 km/s curve are from `galactic-potential.md`; the 1.44 pc tidal radius is computed here and agrees with
the 1.4 pc of `orbital-updates.md`. The local remnant mass densities use the rates at the solar
neighbourhood (`PHENOMENON_DENSITY_PC3`), not the scatter's regional factors.


## The star passes (GEN.185, 2026-10-10)

Boss (2026-10-10 03:01Z) set the order of the whole scatter: (1) quasars, black holes and the other objects above the mass
limit, (2) the stars above the mass limit, (3) mark each sector with the luminosity of the star scattered in it, (4) the stars
above the luminosity limit, skipping every sector marked at or above that limit, (5) the other phenomena (no comets, no rogue
planets). How it is read and built:

- **One mass limit.** The cut of this note (`galaxy_shape.phenomenon_min_mass_solar`, `--phenomenon-min-mass`, 8 to 20 solar
  masses, default 20) now also divides the stars. Passes 1 and 5 are the phenomenon scatter, which draws nothing from the stars,
  so the order between them and the star passes cannot show; `planetgen plan` runs it first, then the star passes.
- **Pass 2, the mass pass** (`bright_stars.scatter` with `mass_range=(limit, None)`, luminosity floor at the brightest white
  dwarf): every star born with at least the limit, bright or not. At 8 solar masses a main-sequence star shines at about 2,000
  solar luminosities, so below a 3,000 floor the mass pass still places the stars the luminosity pass alone would have missed.
- **Pass 3, the marks** are not stored: `store.bright_star_marked_addresses` reads them from the mass pass's own rows (a star
  born at or above the limit and at least as bright as the floor) between the passes. Each sector keeps the brightest star it
  was given, so "equal to or above the setting" is simply "holds one".
- **Pass 4, the luminosity pass** draws the stars born lighter than the limit and at least the floor bright
  (`mass_range=(None, limit)`), in the unmarked sectors only. A marked sector therefore gets no second bright star; this thins the
  bright-star count a little where the heaviest stars sit (their sectors), which is the intent of the rule.
- **The sector's own draw** is lighter than the limit and dimmer than the floor (`MAX_STAR_MASS_SOL`,
  `MAX_STAR_LUMINOSITY_SOL`), its expected count is cut by `placed_star_fraction` (the mass pass's share plus the lighter bright
  share), and the backfill and the staged bands below the floor draw only lighter stars. `galaxy_shape.bright_star_mass_limit_sol`
  (schema v73) records the limit, so a galaxy scattered before it behaves as it did.
- **Cost:** the mass pass places about 2 million stars at 8 solar masses and 69,000 at 20 (the Research Lane 1 numbers),
  against 26 million for the old 1,000 L_sun floor; at the new 3,000 floor the whole star scatter is a fraction of the old one.
- The partition is exact: the luminosity fraction of the stars lighter than a mass plus that of the heavier ones is the whole
  fraction (`tests/test_star_scatter_passes.py`).
