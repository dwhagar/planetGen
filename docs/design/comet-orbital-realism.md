# Comet Realism: Add Elliptical & Parabolic Comet Types

## Context

Comets in planetGen today are represented by exactly one class:
`InterstellarComet` (`src/stellarObjects/roguePlanetData.py:175-318`) — a
standalone, unbound phenomenon with a hyperbolic excess velocity, generated
only via `phenomenonGen.py`, never tied to a star system. There is no
concept of a comet that orbits a star, gets captured, or slingshots off a
planet — that vocabulary (`eccentricity`, `perihelion`, `capture`,
`hyperbolic`, `slingshot`) does not exist anywhere else in the codebase.

The user wants comets to feel more astronomically real. After discussion,
the scope for this pass is narrower than "full N-body capture/slingshot
simulation": add the two missing, physically distinct comet orbit types —
**elliptical (bound, periodic)** and **parabolic (marginal, single-pass)**
— as proper star-associated objects, alongside the existing standalone
interstellar (hyperbolic) comet. A separate, parallel thread is already
building trajectory/wobble-update logic for planets, moons, and binary
stars (recomputing trajectories roughly once per in-game month); true
gravitational capture/slingshot (which requires modeling a close encounter
with a second massive body) is out of scope here and should be revisited
once that thread's N-body/close-encounter primitives exist.

## Comet orbital taxonomy (research summary)

Real comets split cleanly by orbital energy/eccentricity `e`, and this maps
directly onto three implementation categories:

| Class | Eccentricity | Bound to star? | Real-world example | Status in codebase |
|---|---|---|---|---|
| Hyperbolic / interstellar | e > 1 | No — passes through once, never originated in this system | ʻOumuamua, Borisov | **Exists**: `InterstellarComet` |
| Parabolic / near-parabolic | e ≈ 1 | Marginally — originates in the system's own Oort-cloud-analog, but escapes after one perihelion pass | Many "great comets," single-apparition long-period comets | **Missing** — this plan |
| Elliptical | 0 ≤ e < 1 | Yes — periodic, returns every orbit | Halley (e≈0.967, P≈76yr), Jupiter-family comets (P<20yr) | **Missing** — this plan |

Elliptical comets further split by period/origin (Jupiter-family, P<20yr,
low inclination, shaped by repeated Jupiter encounters; Halley-type,
20–200yr, can be high/retrograde inclination; long-period, >200yr, e very
close to 1, isotropic inclination from the Oort Cloud) — useful as
flavor/subtype metadata but not a physics distinction, since all of them
propagate the same way (bound two-body Kepler orbit).

Key physical points that should drive the design:
- A comet's coma/tail activity is driven by solar heating, not by orbit
  type — activity should scale with **perihelion distance**, not with
  which class the comet belongs to (ices sublimate noticeably inside
  ~2.5–3 AU regardless of whether the orbit is elliptical or parabolic).
- Motion along an eccentric orbit is **not** uniform in angle (Kepler's
  second law — a comet near perihelion moves much faster than near
  aphelion). The existing planet/moon phase-advance in `updateOrbits.py`
  linearly advances `orbital_phase_deg` by elapsed time, which is a fine
  approximation for near-circular planetary orbits but breaks down at
  cometary eccentricities (elliptical comets can have e up to ~0.99;
  parabolic is e≈1). Realistic comets need proper Kepler-equation
  propagation (mean anomaly → eccentric anomaly → true anomaly), not the
  simple linear phase used elsewhere.
- Parabolic orbits have no periapsis-to-periapsis period (technically
  infinite) — they must be modeled as a one-shot event (time since/until
  perihelion) rather than with `period_years`, similar in spirit to how
  `InterstellarComet` already tracks a single passage via `is_active`.

## Recommended design

### 1. New star-bound comet types, following the existing belt/field precedent

The codebase already has a precedent for "same object, standalone vs.
star-bound" in asteroids: `asteroidFieldData.AsteroidField` (standalone)
vs. `asteroidData.AsteroidBelt` (star-bound). Mirror that pattern for
comets instead of overloading `InterstellarComet`:

- New module `src/stellarObjects/cometData.py` with a `Comet` class,
  associated with a star (like `AsteroidBelt`), carrying an
  `orbit_type` field: `"elliptical"` or `"parabolic"`.
- Shared physical fields (reuse as-is): `nucleus_diameter_km`,
  `composition` (sampled from `program_constants.COMET_COMPOSITION`),
  `is_active` — same fields already on `InterstellarComet`.
- New orbital fields:
  - `perihelion_distance_au` (q) — always present, drives activity.
  - `eccentricity` — 0 ≤ e < 1 for elliptical, fixed at/near 1.0 for
    parabolic.
  - `inclination_deg`, `arg_periapsis_deg`, `ascending_node_deg` — full
    3D orbit orientation (comets are not confined to a system's ecliptic
    plane the way planets are).
  - `orbital_period_years` — only for elliptical; `None` for parabolic.
  - `time_of_perihelion_passage` (or equivalent phase reference) — replaces
    the simple `orbital_phase_deg` used for planets, since position must be
    derived via Kepler's equation, not linear interpolation.
  - Optional `period_class` flavor tag (`"jupiter_family"`,
    `"halley_type"`, `"long_period"`) derived from `orbital_period_years`
    for description/plausibility text only.

### 2. Position updates: proper Kepler propagation, not linear phase

Add a small orbital-mechanics helper (e.g. in `planetPhysics.py` or a new
`keplerMotion.py`) that solves Kepler's equation for eccentric orbits:
mean anomaly `M = 2π·(t - t_perihelion)/P` (elliptical) or Barker's
equation for the parabolic case → eccentric/true anomaly → distance from
star. Hook this into `updateOrbits.py` as a comet-specific branch rather
than reusing the linear `orbital_phase_deg` advance used for planets/moons.
This is self-contained (pure two-body physics) and does not depend on the
parallel trajectory-wobble work, but should adopt whatever periodic
"recompute trajectory on update" structure that thread introduces, so the
two stay stylistically consistent once it lands — worth a follow-up sync
rather than blocking this work.

### 3. Activity driven by perihelion distance

Replace the flat `INTERSTELLAR_COMET_ACTIVE_CHANCE` roll (used for
interstellar comets, where perihelion is arbitrary/often irrelevant) with
a perihelion-distance-scaled activity chance for star-bound comets —
higher activity chance the closer `perihelion_distance_au` is to the star,
tapering off past ~3 AU. Add this as a new constant/function in
`program_constants.py` alongside the existing `COMET_COMPOSITION` and
active-chance constants.

### 4. Parabolic comets as one-shot events

Parabolic comets should behave like `InterstellarComet` in lifecycle (one
perihelion passage, then gone) but originate from the star's own system
rather than interstellar space — track only `is_active` /
time-since-perihelion, no `orbital_period_years`, and flag for removal or
"expired" status once sufficiently far past perihelion (mirrors the
existing `is_active` pattern rather than inventing a new lifecycle
concept).

### 5. Generation, schema, and plausibility — reuse existing wiring

- **Generation**: unlike `InterstellarComet` (only spawned standalone via
  `phenomenonGen.py`), the new `Comet` type should be spawnable from
  `systemGen.py` as part of normal star-system generation, the same way
  `AsteroidBelt` is.
- **Schema/DB**: add `comets` (+ `comet_composition`) tables in
  `schema.sql`, modeled on `interstellar_comets` (`schema.sql:1218-1253`)
  plus the new orbital columns, with a `star_id` FK; wire insert/select in
  `_db.py` alongside `insert_interstellar_comet` (`_db.py:1185-1223`).
- **Plausibility**: extend `phenomenaPlausibility.py` with checks specific
  to bound comets — eccentricity range per `orbit_type`, positive
  perihelion distance, `orbital_period_years` required iff elliptical.

### 6. Testing

- Unit tests for the Kepler-equation helper (mean→eccentric→true anomaly)
  against known values (e.g. verify Halley-like e≈0.967 orbit position at
  known mean anomalies).
- Extend `test_phenomena.py` / `test_phenomena_plausibility.py` with cases
  for the new `Comet` class and its plausibility rules.
- Extend `test_orbital_motion.py` with a comet-specific case showing
  non-uniform angular speed (faster near perihelion) to confirm the new
  propagation is actually being used instead of linear phase advance.

## Verification

- Run `pytest src/tests/test_phenomena.py src/tests/test_phenomena_plausibility.py src/tests/test_orbital_motion.py` and confirm new comet tests pass.
- Generate a test star system via `systemGen.py` and confirm elliptical/parabolic comets appear with sane orbital elements (e in valid range, perihelion positive, period present only for elliptical).
- Spot-check Kepler propagation manually: for a high-eccentricity test comet, confirm the time spent near perihelion is much shorter than time near aphelion (non-uniform angular speed), unlike the existing linear-phase planet motion.
- Confirm `InterstellarComet` and `phenomenonGen.py` are untouched/still pass existing tests — this is additive, not a replacement.
