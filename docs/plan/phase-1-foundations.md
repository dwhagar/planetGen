# Phase 1: Foundations and fixes

This is one phase of the plan drawn up on 2026-10-01 from Boss's list
of that evening and his research notes. [docs/TODO.md](../TODO.md) is
the master file: it holds each item's full text (what's wrong, where to
look, what "done" means), and its "Plan: phases" section indexes every
phase. This file holds what the phase needs beyond that: its goal, the
order and dependencies, how to split it into build threads, and the
research notes that apply. Where the two disagree, TODO.md wins; when an
item ships, it leaves TODO.md and its row here is deleted in the same PR.

## Goal

Fix what generation and storage get wrong, and lay the shared groundwork
the map and navigation work of phase 2 builds on, so that phase can run
in parallel threads without rewriting the same files twice. Nothing in
this phase changes how the site looks, except the System Map NaN fix.

## Items

### The parallel path (top priority)

Boss (2026-10-02) made this the top priority of the whole plan: it
comes first in this phase, before anything else starts.

| ID | Item | Parent |
|---|---|---|
| PERF.21 (new) | Generation works with any worker count: the parallel path is built, used and tested (bug) |  |
| PERF.22 (new) | On Python 3.12 a run hangs forever when a worker process dies (bug) | PERF.21 |
| TEST.74 (new) | Generation tests at more than one worker | PERF.21 |
| PERF.23 (new) | The bright-star progress bar can end at 101% (bug) |  |

TEST.74 (the generation tests at 1, 2 and 4 workers) lands first so the
14 tests that fail at more than one worker show what PERF.21 has to fix;
PERF.22 (the 3.12 hang in `workQueue._dispatch`) and PERF.23 (a late
layer report counted twice in `_LayerTracker`) are in the same parallel
code and go in the same thread. GEN.39 (reproducible seeds) touches the
same tests; whichever lands second keeps both working.

### Forcing options and prevalence

| ID | Item | Parent |
|---|---|---|
| GEN.48 (new) | Forcing options are impractical for whole sectors; replace them with prevalence controls (bug) |  |
| GEN.49 (new) | `+habitable_world` silently fails on hot stars (bug) | GEN.48 |
| GEN.50 (new) | `-planets +asteroid_belt` still makes an asteroid belt (bug) | GEN.48 |
| GEN.51 (new) | Forcing options only for single-system generation | GEN.48 |
| GEN.52 (new) | Prevalence controls for sector and galaxy runs | GEN.48 |
| ADM.16 (new) | Prevalence controls on the Generate page | GEN.48 |
| TEST.75 (new) | Tests for forcing and prevalence | GEN.48 |
| GEN.53 (new) | The two stars of a binary don't share one age (bug) |  |
| GEN.54 (new) | A `--star-type` secondary gets a mass that doesn't fit its type (bug) |  |

GEN.49 and GEN.50 (single-system forcing that fails or contradicts
itself) first, then GEN.51 (forcing only for single systems) and
GEN.52 (prevalence as a percentage deviation), then ADM.16 (the
Generate page fields); TEST.75 grows with each. GEN.53 and GEN.54 are
both in how a `--star-type` secondary is made (`systemData.py`) and go
in one PR.

### Planet classes and physics

| ID | Item | Parent |
|---|---|---|
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) |  |
| GEN.33 | One class per PR, each with its tests | GEN.28 |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there |  |
| GEN.38 | Rocky rogue planets over 10,000 km are still classed C (bug) |  |
| GEN.29 | Sweep every planet class for sense once the new ones are in (bug) |  |
| GEN.25 | A moon reclassified after its planet moves can be too large for its planet (bug) |  |
| GEN.34 | Gas and ice giants come out too light, so there are no super-Jupiters (bug) |  |
| GEN.35 | Rocky planets only ever get Class D moons (bug) |  |
| GEN.36 | Moon regeneration can produce gas-giant or blacklisted moon classes (bug) |  |
| GEN.37 | 97% of planets land in the cold zone (bug) |  |
| GEN.45 (new) | Check the rogue planet mix of terrestrial and gas giants (bug) |  |

GEN.28 lands first, one class per PR (GEN.33: R and S, then U, W, X, Y
and Z). GEN.38 (big rocky rogues) is likely solved by making S
rogue-eligible. GEN.27 and GEN.29 follow GEN.28. GEN.34 to GEN.37 and
GEN.25 are independent bugs in `planetPhysics.py` and the moon code and
can run alongside. GEN.45 (the rogue planet mix) touches
`roguePlanetData.py` and `ROGUE_PLANET_MASS_BINS` only; if GEN.28 adds
rogue classes, settle the mix after them.

### Galaxy generation

| ID | Item | Parent |
|---|---|---|
| GEN.24 | Generate the galactic core on layer 0 |  |
| GEN.31 | A point just under layer 0's top face lands in layer 1 (bug) |  |
| GEN.32 | Re-running an interrupted bright-star band draws it twice (bug) |  |
| GEN.39 | The same seed can't reproduce the same galaxy (bug) |  |
| GEN.44 (new) | Store each sector's backfill level so finished sectors drop out of any backfill | GEN.40 |
| GEN.46 (new) | Star system names of at most two words (bug) |  |
| GEN.47 (new) | Nebulae almost never appear (bug) |  |

GEN.31 (the layer 0 boundary) before GEN.24 (generate the core on layer
0). GEN.44 (each sector's backfill level) is this phase's galaxy schema
change; it replaces or sits beside `bright_star_blocks` and is what
GEN.40 to GEN.43 and PERF.18 in phase 3 use to skip work, so it should
land with or before GEN.32 (an interrupted band drawn twice), which
touches the same backfill code. GEN.39 (reproducible seeds) needs
Boss's answer first. GEN.46 (two-word names) changes `nameUniqueness.py`
and the name reservation in `_db.py`. GEN.47 (nebulae) is a new
galaxy-scale placement, not a one-line fix.

### Database

| ID | Item | Parent |
|---|---|---|
| DB.2 | Asteroid field and comet composition rows are written but never read (bug) |  |
| DB.3 | resetDb while another process holds id blocks can duplicate primary keys (bug) |  |
| DB.4 | A database with an emptied schema_migrations table is treated as current (bug) |  |
| DB.5 | Several first connections to an empty database race to create the schema (bug) |  |

Independent bugs. DB.4 and DB.5 both touch schema creation in
`_db.get_connection`; one PR. Any galaxy schema change in this phase
(GEN.44, possibly DB.2) takes the next version after main's (50 at the
time of writing), one writer at a time.

### Groundwork for phase 2

| ID | Item | Parent |
|---|---|---|
| MAP.63 | Shared map helpers in one module | MAP.61 |
| MAP.64 | One camera and input controller | MAP.61 |
| NAV.7 | One reference for every object, with its parents |  |
| NAV.8 | Pages and anchors for stars, planets, moons and belts | NAV.7 |
| NAV.9 | Search and locate return references for every kind | NAV.7 |
| TEST.70 | Tests for the map JavaScript |  |
| MAP.57 | The System Map writes NaN or infinite positions into its SVG (bug) |  |
| MAP.88 (new) | Parts of a star system run off the edge of the System Map (bug) |  |
| MAP.89 (new) | System Map: space orbits with a fitted scale and a minimum ring gap instead of plain log |  |

MAP.63 (shared helpers in `static/mapcore.js`) and MAP.64 (one camera
and input controller) change no behavior and must land before the
Galaxy Map redesign, which rewrites the same files. NAV.7 (one
reference for every object, with its parents) is what the picker and
courses of phase 2 build on, with NAV.8 and NAV.9. TEST.70 (tests for
the map JavaScript) should exist before the map rewrite starts, and
MAP.66 and MAP.68 need it. MAP.57 and MAP.88 are small independent
fixes to the System Map (`lib/systemmap.py`), one PR; MAP.89 (fitted
orbit spacing, from the orbit spacing study) changes the same file's
scale, so it goes in the same thread, after them.

### Operations and test fixes

| ID | Item | Parent |
|---|---|---|
| OPS.6 | Admin scripts accept impossible `--mysql-port` values (bug) |  |
| OPS.7 | Update asks to fill a wiped database with population data (bug) |  |
| OPS.8 | Update reloads Apache itself when run as root |  |
| TEST.71 | Intermittent failure in the admin planet-regenerate test (bug) |  |
| TEST.72 | Intermittent failure in the two-step (2FA) sign-in test (bug) |  |
| TEST.73 | Intermittent failure in the parallel galaxy-run interrupt test (bug) |  |
| TEST.76 (new) | A bright-star test breaks on Python 3.9 and 3.10 (bug) |  |
| MAP.90 (new) | The tile-level helper crashes on a subnormal view radius (bug) |  |
| UX.34 (new) | The sector summary calls white dwarfs "B-type" and "A-type" systems (bug) |  |
| OPS.9 (new) | Multi-line messages lose their prefix in the debug log (bug) |  |

Small and independent (TEST.76, MAP.90, UX.34 and OPS.9 are low-priority
finds of the debug-mode bug hunt). OPS.7 and OPS.8 both change `update.sh`,
`update.ps1` and `scripts/deploy-common.*`; one PR.

## Research notes

Boss's research notes (kept in the project's shared files under `todo-tasks/research/`) proposed fixes and numbers. They were checked against the code on 2026-10-01; where they were wrong about the code, the correction is given. Their numbers are starting points to tune, not requirements.

- **Backfill level (GEN.44).** The notes give the semantics as a table:
  -1 never backfilled, 0 fully generated, a positive value the dimmest
  luminosity reached. A backfill at floor L skips a sector at 0 or at a
  level at or below L, and afterwards sets each processed sector to L
  (or 0 once it is fully generated). Correction: the notes put the
  column on `sectors`, but an unfilled sector has no `sectors` row, so
  the level needs its own table keyed by address (or a reworked
  `bright_star_blocks`, which already keeps the level per 3x3x3 block,
  schema v49). The notes' file list (`generationLimits.py`,
  `spaceSector.py`) is not where backfill lives; it is
  `backfill_bright_stars_around` and `_backfill_block` in `generate.py`
  and the block table in `_db.py`.
- **Rogue planet mix (GEN.45).** The notes propose a power-law mass
  function, dN/dM proportional to M^-0.65 from 0.01 Earth masses to 13
  Jupiter masses, giving mostly terrestrial and few giant rogues, and
  suggest an 80/20 terrestrial/gas split. Today's arithmetic is about
  91/9. Correction: there is no `phenomena_config` table; the weights
  are `ROGUE_PLANET_MASS_BINS` and the gas giant threshold in
  `program_constants.py`.
- **Two-word names (GEN.46).** The notes suggest resolving collisions
  inside two words, with catalog numbers joined by a hyphen
  (`Kepler-452`) rather than a third word, and binary suffixes folded
  into the two words. The words come from the uniqueness decorations
  (Greek prefix, "Alpha ... IV", stacked diminutives), so that is where
  the fix goes.
- **Nebulae (GEN.47).** The notes say nebula generation was cut out of
  the sector loop in `generate.py`. Correction: it isn't; nebulae are
  rolled per sector, but at rates that make them almost never appear in
  a 4 pc sector. The notes' proposal still fits the real fix: a
  multi-octave 3D noise density field with Poisson-seeded centers at
  galaxy scale, with these classes as a starting point:

  | Class | Made of | Field threshold | Radius | Near |
  |---|---|---|---|---|
  | Emission (H II) | ionized hydrogen | > 0.72 | 15 to 45 pc | O and B supergiants |
  | Reflection | silicate and carbon dust | > 0.65 | 10 to 30 pc | B and A main sequence |
  | Dark | cold molecular gas | > 0.80 | 5 to 25 pc | low-density voids |
  | Planetary | ejected envelope | event | 0.1 to 2 pc | post-AGB stars, white dwarfs |
  | Supernova remnant | blast shell | event | 5 to 20 pc | pulsars, neutron stars |

  The existing classes and rates are in
  [nebula-and-asteroid-field-classes.md](../design/nebula-and-asteroid-field-classes.md).

## Build threads

Each thread is briefed with its exact item IDs and takes no others.

1. The parallel path (top priority, starts first): TEST.74, PERF.21,
   PERF.22, PERF.23.
2. Planet classes: GEN.28 (as GEN.33's PRs), then GEN.38, GEN.27,
   GEN.29.
3. Physics bugs: GEN.34, GEN.35, GEN.36, GEN.37, GEN.25, GEN.45.
4. Galaxy generation: GEN.31, GEN.24, GEN.44, GEN.32, then GEN.47;
   GEN.46 and GEN.39 (after Boss answers) alongside.
5. Database: DB.2 to DB.5.
6. Map groundwork: MAP.63, MAP.64, TEST.70, MAP.57, MAP.88, MAP.89.
7. Object references: NAV.7, NAV.8, NAV.9.
8. Ops and flakes: OPS.6 to OPS.9, TEST.71 to TEST.73, TEST.76,
   MAP.90, UX.34.
9. Forcing and prevalence: GEN.49, GEN.50, GEN.51, GEN.52, ADM.16,
   TEST.75 (parent GEN.48).
10. Binary pairs: GEN.53, GEN.54.

## Open questions for Boss

- GEN.39: should a seed reproduce a galaxy?
- GEN.46: are existing longer names renamed, or only new ones?
- GEN.44: a new table per sector address, or rework `bright_star_blocks`?
