# Comets: elliptical, parabolic and interstellar

**Status:** built in version 5.31.0 (schema v19), PR #66. Names follow the
IAU-style designations added in 7.31.0 (schema v40). Since 7.46.0 each
comet period class has a reference page (`/classes/comet/<class>`), and
since 7.57.0 an interstellar comet's page shows a rendered nucleus, with a
coma and tail when it is active (`planetgen/web/maps/phenomenonrender.py`). This
document started as the plan for that work; it now describes what the
code does and why.

## What exists

planetGen has three kinds of comet, split by orbital eccentricity `e`:

| Kind | Eccentricity | Bound to a star? | Class | Table |
|---|---|---|---|---|
| Elliptical (periodic) | 0 <= e < 1 | Yes, returns every orbit | `cometData.Comet`, `orbit_type = "elliptical"` | `comets` |
| Parabolic (one pass) | 0.995 to 1.0, propagated as an exact parabola | Marginally: from the system's own outer cloud, leaves after one perihelion | `cometData.Comet`, `orbit_type = "parabolic"` | `comets` |
| Interstellar (hyperbolic) | e > 1 | No, passes through a sector once | `roguePlanetData.InterstellarComet` | `interstellar_comets` |

Star-bound comets belong to a star system, like asteroid belts. Interstellar
comets are sector phenomena, drawn per star at the design rate in
`program_constants.PHENOMENON_DENSITY_PC3["comet"]` (see
`interstellar-object-rates.md`).

### Elliptical comets

- `period_class` is drawn first from `program_constants.COMET_PERIOD_CLASSES`
  (Jupiter-family 50%, Halley-type 30%, long-period 20%). The class sets the
  eccentricity range and the largest inclination.
- Perihelion `q` is drawn from `COMET_PERIHELION_DISTANCE_RANGE_AU`
  (0.05 to 5 AU). The semi-major axis is `q / (1 - e)`, and the period comes
  from Kepler's third law with the host star's mass.
- The class's `period_range_years` is not used by generation. Because the
  period is computed from `q` and `e`, a comet tagged Jupiter-family can come
  out with a period outside 3.3 to 20 years. The tag is descriptive only.
- Motion uses the mean anomaly `mean_anomaly_deg`, which advances linearly
  and wraps at 360. Position comes from solving Kepler's equation
  (`kepler.solve_eccentric_anomaly`, scipy Brent), so a comet moves
  fast near perihelion and slowly near aphelion.

### Parabolic comets

- About 30% of star-bound comets (`COMET_PARABOLIC_CHANCE`).
- `eccentricity` is stored (0.995 to 1.0) but the orbit is propagated as an
  exact parabola with Barker's equation (`keplerMotion.solve_barker_equation`).
- There is no period. `parabolic_mean_anomaly` advances linearly and does not
  wrap; it starts near zero (around the one perihelion passage) and can be
  negative (still approaching).
- Nothing removes or expires a parabolic comet after its pass. It keeps
  receding on every orbit update. The plan's "flag for removal" step was not
  built.

### Activity

Whether a comet shows a coma and tail (`is_active`) depends on perihelion
distance, not on orbit type: `cometData._activity_chance` falls linearly from
0.9 at `q = 0` to 0.05 at `q >= 3 AU` (`COMET_ACTIVITY_*`). Interstellar
comets keep a flat roll (`INTERSTELLAR_COMET_ACTIVE_CHANCE = 0.5`).

### Names

`cometData.comet_designation` (7.31.0, schema v40) names a star-bound comet
`P/<host>-<n>` when it is elliptical with a period under 200 years, and
`C/<host>-<n>` otherwise. `<host>` is the star it orbits and `<n>` counts
from 1. The name follows a rename of its star (`rename_comet_designation`).
Interstellar comets are `I/<sector>-<n>`.

### Storage and updates

- `comets` and `comet_composition` (schema v19, `_migrate_v18_to_v19`). The
  table stores `perihelion_distance_km` (km in the database, AU in Python),
  the full orientation (`inclination_deg`, `arg_periapsis_deg`,
  `ascending_node_deg`), `mean_anomaly_deg` or `parabolic_mean_anomaly`, and
  the current position and speed.
- `StarSystem` holds `comets` (and `secondary_comets` for a wide pair's
  second star), separate from `planets`, because a comet's distance changes
  continuously instead of sitting in an orbital slot.
- `SystemConfig.COMETS` is a tri-state flag; `planetgen system` exposes it.
- `planetgen.cli.orbits` calls `_db.advance_comet_orbits` after the planet and
  moon pass. It reads every comet row, advances the anomaly, solves Kepler or
  Barker in Python and writes all rows back with one `executemany`.
- Since 7.37.0 the whole system, comets included, also moves along its
  galactic orbit and can change sector.

### Tests

`src/tests/test_kepler_motion.py` covers the solvers (Kepler's equation,
Barker's equation, faster motion near perihelion, vis-viva speed).
`src/tests/test_comet_data.py` covers generation ranges, perihelion bounds,
round-tripping and the page text. `test_phenomena.py` and
`test_phenomena_plausibility.py` cover interstellar comets. The plan's
proposed comet cases in `test_orbital_motion.py` were placed in
`test_kepler_motion.py` instead.

## Why it works this way

- **Three kinds, split by eccentricity.** Real comets fall into these three
  dynamical groups, and they need different math: elliptical orbits use
  Kepler's equation, parabolic ones Barker's equation, and interstellar ones
  a fixed hyperbolic pass. Before 5.31.0 only the interstellar kind existed,
  so no comet could orbit a star (CHANGELOG 5.31.0).
- **A separate `Comet` class, not an option on `InterstellarComet`.** The
  plan followed the asteroid precedent: `AsteroidField` (standalone) and
  `AsteroidBelt` (star-bound) are separate classes. A bound comet needs the
  host's mass and a full orbit, which an interstellar comet does not.
- **Kepler propagation instead of a linear phase.** Planets and moons
  advance `orbital_phase_deg` linearly, which is close enough for nearly
  circular orbits. At comet eccentricities (up to 0.999) that would be badly
  wrong, because most of the orbit is spent far out. The cost is that the
  update cannot be one set-based SQL `UPDATE` as it is for planets; it runs
  in Python. The table is small (1 to 3 comets per system,
  `SYSTEM_COMET_COUNT_RANGE`), so this was accepted, and the writes are
  batched (CHANGELOG 5.31.0, "advance_comet_orbits now batches").
- **Activity from perihelion distance.** Ices sublimate when a comet comes
  within about 3 AU of a Sun-like star, whatever its orbit type. A flat roll
  would let a comet that never comes near its star show a tail.
- **IAU-style names.** Boss asked that comets and asteroid fields get names
  that say something about them in a standard way (GEN.13); the real IAU
  `P/`, `C/` and `I/` prefixes do that (CHANGELOG 7.31.0).

### Alternatives not taken

- **Gravitational capture and planetary slingshots.** Both need a real
  three-body close encounter. They were left out of 5.31.0 on purpose, to be
  revisited once the project has close-encounter primitives. None exist yet.
- **Expiring parabolic comets.** Planned, not built. The reason it was
  dropped is not recorded.
- **Plausibility checks for bound comets.** The plan asked
  `phenomenaPlausibility.py` to check eccentricity per orbit type and that
  only elliptical comets have a period. Those checks were not added; the
  ranges are enforced by generation and by `test_comet_data.py` instead.
  The reason is not recorded.

## Corrections made to the original plan text

The earlier version of this file was the pre-build plan. It named
`phenomenonGen.py` and `systemGen.py`, which were merged into `planetgen`
in 5.35.0; it cited line numbers in `roguePlanetData.py`, `schema.sql` and
`store.py` that no longer match; and it described a "time of perihelion
passage" field, where the code stores a mean anomaly instead.

## Ejected comets are by design (GEN.182)

Boss saw comets on the System Map that never came back. Checked: about 30%
of star-bound comets are parabolic on purpose (single apparition, see
above), and the 3D map draws them as an open path and moves them correctly at
any speed (`orbitpositions.js`, Barker's equation). Every elliptical comet is
bound: the widest is `q / (1 - e)` with `q <= 5 AU` and `e <= 0.999`, an
aphelion under 10,000 AU, well inside the star's Hill sphere; a test draws
many and checks it. A long-period comet's orbit can also simply be longer than
the run so far (a period of thousands of years at 100 years per second).
The info panel now says "Bound, returns every N years" or "Unbound: passes
the star once and does not return", with the perihelion and eccentricity.
