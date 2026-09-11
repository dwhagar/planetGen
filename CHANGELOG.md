# Changelog

## [5.18.0] - 2026-09-11

### Added
- **`examples/maintenance/`: systemd timer for `updateOrbits.py`.** An
  Ubuntu/Debian-native alternative to the raw crontab line
  `docs/database-schema.md` already documented for running the periodic
  orbital-motion update ("once a month or so"). `planetgen-orbits@.service`/
  `.timer` are a systemd *template* unit -- the instance name (e.g.
  `planetgen-orbits@planetgen.timer`) selects which database gets
  updated, so a deployment with more than one `PLANETGEN_MYSQL_DATABASE_PREFIX`
  schema enables one timer instance per database rather than needing a
  separate script per database. `install-maintenance-timer.sh` installs
  both units, writes `/etc/planetgen/maintenance.env` (mode 600) from
  `maintenance.env.example` for the shared read-write MySQL credentials
  (skipped if that file already exists, so it never clobbers credentials
  already set up), and enables the timer for each database name given on
  its command line (defaulting to `$PLANETGEN_MYSQL_DATABASE`, or
  "planetgen"). Output is captured by journald automatically, so there's
  no logfile/logrotate entry to maintain the way the crontab example
  needs.

## [5.17.0] - 2026-09-11

### Added
- **Binary mutual-orbit position.** `updateOrbits.py`/`_db.advance_orbital_phases`
  now recomputes `star_systems.binary_mutual_position_x_km`/`_y_km`/`_z_km`
  (the secondary star's Cartesian position relative to the primary) every
  time it advances `binary_mutual_orbital_phase_deg` -- the same "position
  has no independent update of its own, it just has to move whenever
  phase does" treatment `position_x/y/z_km` already gets for planets/
  moons. Derived via `utils.orbital_position_au` from
  `binary_separation_km` and the mutual orbit's own (non-near-ecliptic,
  full `[0, 180)`/`[0, 360)`-range) inclination/ascending-node/phase --
  the pair's true "direction of orbit" was already fully captured by
  those v13 orbital elements; this just keeps the derived position
  correctly in sync with them as time passes, rather than only at
  generation. Schema v14, with a migration that backfills real position
  values for existing binary rows from their own already-stored
  separation/orbital-element columns.

## [5.16.0] - 2026-09-11

### Added
- **Star motion: galactic orbit phase, plus binary mutual orbit.** Stars
  now get the same floating-point update guard planets/moons already
  have, and binary pairs get a real orbit around each other, both
  actively advanced over time.
  - Every star gains `galactic_orbital_phase_deg` (its current angular
    position around the galactic center, randomly rolled at generation --
    the same role `orbital_phase_deg` plays for a planet/moon) and
    `galactic_min_update_interval_years` (the same guard formula,
    `utils.minimum_update_interval_years` applied to
    `galactic_orbital_period_gy * 1e9` years). `_db.advance_orbital_phases`
    now advances this phase too, guarded by its own interval, reversing
    v12's "stars have no periodic update mechanism" scoping note -- they
    do now. Both stars of a binary pair, and `star_systems`' own mirrored
    `binary_galactic_orbital_phase_deg`/`binary_galactic_min_update_interval_years`,
    always carry the identical value: a binary's AU-scale separation is
    negligible next to its light-year-scale galactic orbit, so the pair
    moves around the galaxy together, not independently
    (`StarSystem.__init__` rolls the phase once and threads it to both
    stars and the proxy).
  - Binary pairs also get their own **mutual orbit** -- the two stars
    circling their common barycenter, entirely separate from (and vastly
    faster than) the galactic orbit above:
    `star_systems.binary_mutual_orbital_period_years`/`_speed_kms`
    (Kepler's third law / circular-orbit speed --
    `planetPhysics.calculate_orbital_period_years`/
    `utils.circular_orbital_speed_kms`, the same formulas a planet's own
    orbit already uses, applied to the pair's separation/combined mass),
    `_inclination_deg`/`_ascending_node_deg`/`_phase_deg` (the same
    `utils.orbital_position_au` orbital-element convention planets/moons
    use for "direction", but drawn from the full `[0, 180)`/`[0, 360)`
    range -- a binary's mutual orbital plane has no protoplanetary-disk
    reason to prefer any alignment), and its own
    `binary_mutual_min_update_interval_years` guard. `advance_orbital_phases`
    advances `binary_mutual_orbital_phase_deg` the same way.
  - Schema v13, with a migration that backfills every real-derivable
    value (both `*_min_update_interval_years` guards, plus the mutual
    orbit's period/speed) from existing rows' own already-stored data;
    the phase/orientation columns with no derivable "correct" value get
    an arbitrary `0` placeholder, the same treatment v9 already gives
    pre-existing planets'/moons' orbital orientation.

## [5.15.0] - 2026-09-10

### Added
- **Floating-point update guard.** Every planet and moon now stores
  `min_update_interval_years` -- the shortest `elapsed_years` worth
  calling `_db.advance_orbital_phases` for. Below this threshold, the
  phase delta `elapsed_years` would add is smaller than
  `orbital_phase_deg`'s own IEEE 754 double-precision resolution, so
  `MOD(orbital_phase_deg + delta, 360)` is guaranteed to round right back
  to the exact value already stored -- a wasted write that changes
  nothing. Derived purely from `period_years`
  (`utils.minimum_update_interval_years`: `period_years *
  math.ulp(360.0) / 360`, the coarsest representable step anywhere in
  `orbital_phase_deg`'s `[0, 360)` range). Not a narrative/display stat --
  purely a guard value: `advance_orbital_phases` now skips a row's
  `UPDATE` entirely (not just a no-op write, no attempt at all) when a
  call's `elapsed_years` is below it. In practice this floor sits many
  orders of magnitude below any realistic elapsed time
  (`updateOrbits.py` runs "once a month or so"), so it exists for
  correctness against a caller advancing time in much smaller steps, not
  because today's actual usage pattern comes close to tripping it.
  Scoped to `planets`/`moons` only -- `stars`' galactic-orbit values are
  fixed forever at generation time, with no periodic update mechanism to
  guard. Persisted as `planets`/`moons.min_update_interval_years` --
  schema v12, with a migration that backfills real derived values (not a
  placeholder) for existing rows.

## [5.14.0] - 2026-09-10

### Added
- **Planet/moon position.** Every generated planet and moon now gets an
  actual 3D Cartesian position (`Planet.position_x/y/z`, in AU), relative
  to its orbital anchor -- the star (or, for a binary system, the
  `BinaryStarProxy` standing in for the system's combined center) for a
  planet, the parent planet for a moon -- continuing the same "each body
  positioned relative to its immediate primary" hierarchy
  `docs/design/galaxy-coordinate-system.md` already uses one level up for
  sectors/systems relative to the galactic center. Derived from the
  existing orbital-motion elements (`orbital_inclination_deg`/
  `orbital_ascending_node_deg`/`orbital_phase_deg`, schema v9) via the
  standard circular-orbit-to-Cartesian transform
  (`utils.orbital_position_au`). Also adds `Planet.orbital_speed_kms`, a
  circular orbit's constant tangential speed (`utils.
  circular_orbital_speed_kms`, `v = 2*pi*r/T`) -- exact here, unlike the
  galactic orbit's rotation-curve model, since a planet's/moon's period is
  already known exactly from Kepler's third law. Surfaced as a new
  "Speed" row in every planet's/moon's rendered data table. Persisted as
  `planets`/`moons`.`position_x_km`/`_y_km`/`_z_km`/`orbital_speed_kms` --
  schema v11, with a migration that backfills real derived values (not a
  placeholder) for existing rows.
- `updateOrbits.py`/`_db.advance_orbital_phases` now recomputes position
  in lockstep with `orbital_phase_deg` as real time passes (a single
  set-based SQL `UPDATE` per table, matching phase advancement's own
  performance characteristics); `StarSystem.validate_system` recomputes
  position/speed too whenever it corrects a planet's `distance` post-hoc
  to resolve an orbital overlap, closing the same kind of staleness gap
  its own docstring already flagged for `period` before this.

## [5.13.0] - 2026-09-10

### Added
- **Galactic orbit.** Every generated star system now gets a circular
  orbital speed and period around the galactic center
  (`Star.galactic_orbital_speed_kms`/`galactic_orbital_period_gy`,
  `utils.calculate_galactic_orbit`), derived from the system's actual
  distance from the galactic center where known (a sector-placed system)
  or the same fixed `physical_constants.GALACTIC_CENTER_DISTANCE_LY`
  fallback `system_perimeter`/`heliosphere_radius` already use otherwise.
  Uses a simple rotation-curve model
  (`GALACTIC_ROTATION_FLAT_VELOCITY_KMS`/`GALACTIC_ROTATION_CORE_RADIUS_PC`)
  rather than a Keplerian point-mass orbit around the galaxy's total mass,
  which would overshoot Sol's real ~220-240 km/s orbital speed by roughly
  4x -- calibrated against Sol's own distance (~206 km/s, ~236 million
  years per orbit, both close to the real Sun's measured values). Surfaced
  as a new "Galactic Orbit" row in every star's/binary pair's rendered data
  table, and persisted as `stars.galactic_orbital_speed_kms`/
  `galactic_orbital_period_gy` (plus the `star_systems.binary_galactic_orbital_*`
  equivalents for binaries) -- schema v10, see `docs/database-schema.md`.

## [5.12.0] - 2026-09-10

_Originally developed and released as `5.9.0` on a separate branch; renumbered
on merge since `5.9.0` was independently already used below for the same-day
Class W removal. No functional difference from the original release._

### Changed
- **`../src/html/` is now a thin frontend for the Flask API instead of a
  direct MySQL client.** Every CGI page (`index.py`, `browse.py`,
  `sector.py`, `system.py`, `nav.py`, `search.py`, `galaxy.py`) fetches
  its data from `GET /api/...` via a new stdlib-only HTTP client
  (`html/lib/apiclient.py`, `PLANETGEN_API_BASE_URL`) rather than opening
  its own read-only MySQL connection -- the database-querying logic those
  pages used to duplicate now lives once in `queryDb.py`, shared with the
  API. `html/lib/dbutil.py` is gone; its non-DB formatting helpers
  (`esc`/`linkify_location`/`format_density`) moved to a new
  `html/lib/fmt.py`. The API itself gained what this needed: `?db=` on
  every route (multi-schema, matching the browser's own picker --
  `stellarObjects._db.list_databases`/`resolve_database`, moved out of
  `html/lib/dbutil.py` into the package itself), `GET /api/databases`,
  `GET /api/galaxy/sectors`, `GET /api/search` (the full faceted-search
  query layer, ported from `html/search.py`), `sector_id=none` on
  `GET /api/systems` (standalone systems), and richer `GET /api/sectors`/
  `GET /api/sectors/<id>`/`GET /api/systems/<id>` responses (display-
  ready fields -- ids, quadrant/location, star summaries, galaxy
  placement -- distinct from `stellarObjects._db.load_sector`/
  `load_star_system`'s *generation* object graph, still reachable the
  same way). Apache's example vhost
  (`examples/apache/planetgen.conf.example`) mounts the API at `/api/`
  (`WSGIScriptAlias`, in its own `WSGIDaemonProcess`) on the same vhost
  that serves `../src/html/` -- landed alongside `html/api/`'s own move
  into the same `DocumentRoot` (`[5.10.0]`), resolving the open question
  `TODO.md` left from the 5.3.2/5.3.3 file-system cleanup from two
  independent, orthogonal angles at once (where the API's source lives,
  and how `html/` gets its data) that happen to combine cleanly. See
  `docs/api.md` and `docs/html-interface.md`.
- System Map: a body with `life_chemical` set (habitable) now gets a
  small green badge on its marker, surfaced in the info panel's "Life
  Chemistry" field too (`html/lib/systemmap.py`'s `has_life`,
  `static/systemmap.js`) -- the one piece of state the existing
  circle+letter marker couldn't show at a glance. `docs/TODO.md`'s long-
  stale "sprite-based graphical system view" item is resolved: the System
  Map (an interactive scaled-orbit SVG diagram, not bitmap sprite art)
  already satisfied it, per that module's own docstring.

## [5.11.1] - 2026-09-10

### Changed
- **Moon tidal locking is now physics-based, not a flat probability.**
  [5.11.0] gave a moon a flat 75% (`MOON_TIDAL_LOCK_PROBABILITY`) chance
  of being tidally locked; replaced with an actual tidal-despinning
  timescale estimate (new `planetPhysics._tidal_locking_timescale_seconds`,
  the standard simplified Murray & Dermott formula: `t_lock = (2Q/15k2) *
  omega0 * a^6 * m_moon / (G * M_primary^2 * R_moon^3)`, with fixed
  representative `Q`/`k2` values for a rocky/icy body -- new
  `physical_constants.MOON_TIDAL_DISSIPATION_Q`/`_LOVE_NUMBER_K2`).
  Verified directly against three real systems spanning many orders of
  magnitude: ~47 million years for the real Earth-Moon system (real
  estimates: tens of millions of years), ~90 years for Mars/Deimos
  (consistent with Deimos being locked given its tiny size), and ~1
  billion years for Saturn/Iapetus (consistent with real estimates of a
  billion-year-plus despinning time for that unusually slow case). A
  candidate (pre-locking) rotation period is drawn first, same as any
  planet; a moon actually ends up locked only if that timescale is
  shorter than the system's age (`star.age`, the best available proxy --
  planets/moons don't carry an independent age of their own).
- **Moon orbital distance is now drawn log-uniformly, not
  linearly-uniformly.** Surfaced by the physics-based tidal locking
  above: `generate_moons`' distance range can span many orders of
  magnitude (its outer bound reaches 1/5 of the parent planet's own Hill
  radius -- tens to hundreds of millions of km for a large planet, far
  beyond where any real large moon actually orbits, e.g. our Moon at
  384,400 km), and a plain `random.uniform` over that range spends almost
  all its density in the single largest order of magnitude -- nearly
  every generated moon landed implausibly far out, so almost none had
  time to tidally lock (measured: ~1.3% of moons locked). Real moon
  systems are much closer to log-spaced (e.g. the Galilean moons run
  421,700 / 671,100 / 1,070,400 / 1,882,700 km, each roughly 1.5-1.6x the
  last), which log-uniform sampling matches far better while still
  allowing occasional genuinely distant/irregular moons. Measured effect:
  ~1.3% -> ~14% of moons locked over the same sample size -- still a
  minority overall (this generator's moons span a much wider population
  than just the handful of large, close, well-known real moons that
  dominate popular intuition about "most moons are locked"), but an order
  of magnitude more of them landing close enough to plausibly have
  locked.

## [5.11.0] - 2026-09-10

### Added
- **Orbital and rotational motion.** Every planet/moon now has a real 3D
  orbital orientation and a rotation period, plus a live position that a
  new, separately-run script advances over time -- resolves
  `docs/TODO.md`'s "Introducing realistic orbital paths and speeds..."
  entry.
  - `Planet` gains four new attributes (`planetPhysics.
    generate_orbital_motion_properties`): `orbital_inclination_deg`/
    `orbital_ascending_node_deg` (fixed at generation time -- together
    they orient this generator's circular-orbit model in 3D; planets draw
    from a tighter, real-solar-system-like range than moons, which can be
    tilted much further), `orbital_phase_deg` (this body's current
    position angle around its orbit -- the one field that changes after
    generation), and `rotation_period_hours` (axial "day length" -- a
    static descriptive stat; no rotational phase is tracked, by design,
    since nothing consumes "which side currently faces the primary"). A
    moon has a real chance (75%, `MOON_TIDAL_LOCK_PROBABILITY`) of being
    tidally locked (rotation period equal to orbital period) -- the norm
    for real large moons, not a rare special case; otherwise rotation
    period is drawn from a `body_type`-appropriate range (`physical_constants.
    ROTATION_PERIOD_RANGE_HOURS`).
  - **New `src/updateOrbits.py`.** Advances every planet's/moon's
    `orbital_phase_deg` in the configured database based on real elapsed
    time since the last run (`stellarObjects._db.advance_orbital_phases`,
    a single set-based `UPDATE` per table, not a per-row Python loop) --
    meant to be run periodically (e.g. cron, "once a month or so"), not
    on every generation run. A new `orbit_simulation_state` singleton row
    tracks when it last ran, measured server-side via `TIMESTAMPDIFF`
    rather than trusting the calling process' own clock to agree with the
    database server's (`get_orbit_update_elapsed_years`). The first run
    against a database just establishes that reference point (nothing to
    advance yet) rather than guessing a start time.
  - **Orbital period fix, needed for this feature to mean anything.**
    `Planet.period` used to be `sqrt(distance_au^3)` unconditionally --
    correct Kepler's-third-law shorthand only for a 1-solar-mass primary.
    For an ordinary planet orbiting a star of very different mass, and
    *especially* for a moon (whose real primary is its parent planet, not
    the grandparent star `self.star` still points at for other purposes
    like life chemistry) this was wrong by orders of magnitude -- a moon's
    period came out as if it orbited its host star directly at that same
    tiny distance. New `planetPhysics.calculate_orbital_period_years(distance_au,
    primary_mass_kg)` takes the actual primary's mass
    (`Planet.__init__`'s new `primary_mass_kg` parameter, threaded in
    from `generate_moons` as the parent planet's own `mass` for a moon,
    defaulting to `star.mass` for an ordinary planet -- transparently
    correct for a binary system too via `BinaryStarProxy.mass`'s existing
    effective-mass property). Also fixed a related latent staleness bug
    found while testing this: `StarSystem.validate_system` adjusts a
    planet's `distance` after generation to resolve orbital overlap and
    already recomputed atmospheric conditions to match, but never
    recomputed `period` -- now it does.
  - Schema v8 -> v9: `orbital_inclination_deg`/`orbital_ascending_node_deg`/
    `orbital_phase_deg`/`rotation_period_hours` on `planets`/`moons`, and
    the new `orbit_simulation_state` table. `stellarObjects._db.
    _migrate_v8_to_v9` is the first real per-version migration step of the
    MySQL era (every database before this one started fresh, already at
    the then-current schema) -- existing rows default to `0`/`24` (an
    arbitrary but harmless placeholder), since this generator never ran
    its actual random generation for them; every body generated from this
    point on gets real values. See `docs/database-schema.md`.

## [5.10.1] - 2026-09-10

### Fixed
- **Subdwarf (Yerkes VI) age modeling could produce a pre-Big-Bang star.**
  `Star._calculate_initial_star_age_and_lifespan` routed every non-main-
  sequence, non-white-dwarf Yerkes class (giants, subgiants, bright
  giants, supergiants, hypergiants, *and* subdwarfs) through the same
  model: derive main-sequence lifespan from the star's own already-
  generated mass, then draw age from the post-main-sequence window. That
  model's mass-sampling side already got a reject-and-resample guard
  against pre-Big-Bang progenitor masses ([5.3.1]) for every class it
  applies to -- except Yerkes VI, deliberately excluded because its
  *entire* allowed mass range (0.1-0.8 Msun) implies a main-sequence
  lifespan longer than the age of the universe, so every draw would have
  been rejected. That left subdwarfs still running through the age *model*
  itself, unguarded, which routinely produced ages of hundreds of billions
  of years -- because real subdwarfs (sdB/sdO) are thought to form via
  binary mass-stripping near the tip of a lower/intermediate-mass
  progenitor's red-giant branch, not single-star post-main-sequence
  evolution, so this star's own post-strip mass was never a valid stand-in
  for a progenitor's main-sequence lifespan in the first place -- the
  model was wrong for this class, not just missing a guard rail (flagged
  as an open design question in `docs/TODO.md`'s "Future ideas" ever
  since).
  - Yerkes VI now has its own age-generation branch, entirely independent
    of the star's own mass: age is drawn directly from a dedicated,
    old-population-biased range (new `program_constants.SUBDWARF_MIN_AGE_GY`
    /`SUBDWARF_MAX_AGE_GY`, 1.0-13.5 Gy, with the same young/old
    `SystemConfig.AGE` bias every other branch applies), reflecting that
    real subdwarf progenitors are typically old, low-mass population stars
    (a short-lived, higher-mass star wouldn't have had time to reach the
    RGB tip and get stripped). Lifespan is that age plus a short
    remaining-phase window (new `SUBDWARF_REMAINING_PHASE_MIN_GY`/
    `_MAX_GY`, 0.05-0.3 Gy) rather than a separately-derived value, since
    the real core-helium-burning subdwarf phase is short relative to the
    age itself.
  - New regression test `test_subdwarf_age_never_exceeds_universe_age`
    (`test_star_matrix.py`) locks this in across every spectral letter.

## [5.10.0] - 2026-09-10

### Changed
- **Moved the Flask API (`src/api/`) into `src/html/`.** Resolves the
  open question `docs/TODO.md` had carried since the 5.3.2/5.3.3
  file-system cleanup: the API now lives at `src/html/api/`, served from
  the same checkout/deployment tree as the interim CGI browser instead of
  a second, separately-tracked location. `src/wsgi.py` moved alongside it
  to `src/html/wsgi.py` (same "lives next to `api/`, so `import api` just
  works via `sys.path[0]`" property as before, plus an explicit
  `sys.path` entry for `src/` now that `queryDb`/`stellarObjects` are a
  directory farther away). `pytest.ini`'s `pythonpath` gained `src/html`
  so `test_api.py`'s `from api...` imports keep resolving unchanged.
  `examples/apache/planetgen.conf.example` gained a `WSGIScriptAlias /api`
  pointing at `html/wsgi.py`, plus deny-all `<Directory>`/`<Files>` blocks
  for `html/api/` and `html/wsgi.py` itself (mirroring the existing
  `html/lib/` block -- both hold source that's imported, never meant to
  be requested directly). Along the way, corrected a pre-existing
  inaccuracy in `docs/api.md`'s Apache deployment guidance: the vhost's
  `SetEnv` directives configure the CGI browser (`mod_cgi`/`mod_cgid`
  copies them into each script's real process environment) but never
  reach `os.environ` under `mod_wsgi` -- `PLANETGEN_MYSQL_*` for the API
  needs to come from the Apache service's own process environment instead
  (e.g. `/etc/apache2/envvars`).

## [5.9.1] - 2026-09-10

### Changed
- **Class P climate tuning.** Gave Class P the same
  `atm_molar_density_range`/`atm_density_range`/`greenhouse_multiplier_range`
  treatment the M/O/H/K/L/N/E/F/G/V pass ([5.3.7]) gave the other habitable
  classes; P had stopped at just its own `albedo_range` that pass. No
  single real-world analog for P, so this isn't chasing a target delta the
  way M/K/N are -- instead the new ranges make the class's own "cold,
  glaciated"/"thinning with age" flavor text physically real: molar
  density stays near Earth's real value (P's atmosphere text names
  oxygen/nitrogen/argon, not a heavier CO2-like mix), `atm_density` is set
  thin, and `greenhouse_multiplier` is set weak so the cold comes from
  genuine physics on top of the class's already-tuned high albedo, not
  albedo alone. Verified via `climate_tuning_cli.py --class P`: mean
  surface_temperature ~219K (well below freezing, clearly colder than
  Class M's ~286K), mean atmospheric_pressure ~10.4kPa (~0.1 atm) over a
  400-sample run across the full host-star grid.

## [5.9.0] - 2026-09-10

### Removed
- **Class W ("a tidally locked world with extreme temperature
  variations") removed entirely.** Its day/night-split identity can't be
  produced from a single global `surface_temperature` scalar under this
  generator's climate model -- no per-class range (albedo, molar density,
  greenhouse multiplier, or atmosphere density) reaches it without an
  actual dayside/nightside model this generator doesn't have, flagged as
  a known gap during the [5.3.7] climate-tuning pass and left open in
  `docs/TODO.md`'s "Investigate Further" section ever since. Rather than
  build a whole day/night thermal model for one class, cut it entirely --
  removed from `PLANET_CLASSES`, `PLANET_CLASS_PROBABILITIES` (its
  0.0001 weight just dropped; these are relative weights, not a
  normalized distribution, so nothing needed redistributing),
  `HABITABLE_PLANET_CLASSES`, and `MOON_BLACKLIST`
  (`program_constants.py`), plus its color entry in
  `html/lib/systemmap.py`'s `_CLASS_COLORS`.

## [5.8.3] - 2026-09-10

### Changed
- **Provisional sector designation is now a single packed hex number.**
  `provisional_sector_designation` (added in [5.8.2]) dropped its
  `"R<ring>-Q<quadrant>-<slot>"` letter/dash format for a plain
  bit-packed hex integer -- `ring` in the high bits, `quadrant - 1` in
  the next 2, the raw `shell_slot_index` in the low 32
  (`DESIGNATION_SLOT_BITS`/`DESIGNATION_QUADRANT_BITS`), e.g.
  `"1500002EE0"` (was `"R5-Q2-2EE0"`). Genuinely reversible back to
  `(ring, quadrant, shell_slot_index)` now too -- the fixed-width bit
  fields have no ambiguous boundary the way concatenating separately-
  sized hex numbers would.

## [5.8.2] - 2026-09-10

### Added
- **Provisional sector designations for un-generated/unvisited addresses.**
  `stellarObjects.galaxyGeometry.provisional_sector_designation(shell_index,
  shell_slot_index, edge_pc, edge_ly)` builds a short, human-readable label
  for a `(shell_index, shell_slot_index)` sector address --
  `R<ring>-Q<quadrant>-<slot>`, Ring and slot index in uppercase hex, Quadrant
  a plain 1-4 digit (e.g. `"R5-Q2-2EE0"`) -- the way a real astronomical
  catalog gives a not-yet-fully-characterized object a provisional name from
  its position rather than waiting for one. `sector_ring`/`sector_quadrant`
  (new, duplicating `html/lib/galaxymap.py`'s identically-named Ring/Quadrant
  concept in the dependency-free `galaxyGeometry` module so `galaxyGen.py`'s
  CLI doesn't need to import the web front-end layer) back it, and are both
  O(1) -- no scan over a shell's other slots, consistent with this codebase's
  galaxy-skeleton design principle of never doing per-sector work
  proportional to a shell's slot count. `galaxyGen.py`'s `--shell`/
  `--center-sector` generation modes now print this designation alongside
  each newly saved sector's procedural name and raw shell/slot numbers.

## [5.8.1] - 2026-09-10

### Added
- **NAV Map: a rendered plot of the NAV feature's origin/destination/route.**
  `queryDb.nav_between` now also returns `origin_position`/
  `destination_position` (the light-year positions its `direct` course was
  computed from) and, when a route was found, each hop's own position
  (`route["positions"]`) -- surfaced from `GET /api/nav`'s JSON response
  too (see `docs/api.md`'s updated "NAV" section). `src/html/lib/navmap.py`
  builds a new "NAV Map" panel on `html/nav.py`'s results page from that
  data: a flat, static, top-down SVG plot of the galactic X-Y plane --
  origin and destination as labeled, clickable points, a dashed line for
  the direct course, and a solid polyline through the optimal route's
  intermediate hops when one exists. Modeled on `galaxymap.py`'s flat 2D
  SVG rather than `starmap.py`'s rotatable 3D CSS scene: like a galaxy
  Quadrant, it's deliberately blind to altitude (the course panel's own
  Altitude figure already covers that axis), auto-scaled to whatever
  points it's given (no fixed sector size to normalize against) with one
  uniform light-years-per-pixel ratio on both axes so azimuth angles
  aren't visually distorted, plus a `+X` compass tick tying the plot's
  orientation to the course panel's own azimuth convention and a
  light-year scale-bar legend. Closes the rendered-image gap the original
  NAV feature (CHANGELOG [5.8.0]) left open in `docs/TODO.md` and
  `src/api/routes.py`'s `systems_near` TODO comment.

## [5.8.0] - 2026-09-10

### Added
- **NAV feature: course, distance, and optimal routing between two
  systems.** Three new pure/query modules plus an API endpoint and a web
  page:
  - `stellarObjects/navigation.py` -- `course_between` (Euclidean
    distance plus galactic-plane-relative azimuth/altitude: azimuth from
    +X in the galactic X-Y plane, altitude as elevation above/below that
    plane) and `warp_travel_times` (`velocity_multiple_of_c = warp_factor
    ** (10/3)`, reported at warp 1/3/6/9, formatted via the existing
    `utils.years_to_time_string`).
  - `stellarObjects/navGraph.py` -- `build_knn_adjacency` (a symmetrized
    k-nearest-neighbor adjacency graph over a `{id: (x, y, z)}` position
    set) and `shortest_path` (Dijkstra) for the "optimal route via
    adjacent systems" half of NAV.
  - `queryDb.nav_between` -- resolves NAV availability between two
    systems (unavailable if either has no sector; same-sector always
    available; cross-sector only when both sectors have a galaxy
    placement), combining a sector's galaxy-frame center (parsecs) with a
    system's sector-local offset (milliparsecs) into one absolute
    position (no such combinator existed before this), then returns the
    direct course plus the optimal route.
  - `GET /api/nav?from=<id>&to=<id>` -- JSON endpoint over `nav_between`,
    `404` for an unknown system id, `400` (`NavUnavailable`) when NAV
    doesn't apply to the pair. See `docs/api.md`'s new "NAV" section.
  - `src/html/nav.py` -- a destination picker (a `<select>` of the
    origin's own sector-mates, plus a typed destination-id field when
    cross-sector NAV is available) and a results page (direct course,
    warp travel times, and the hop-by-hop optimal route, each hop linking
    to `system.py`). `system.py` links here ("Navigate from here")
    whenever a system has a sector.

## [5.7.0] - 2026-09-09

### Changed
- **Merged the parallel Phase 4 "galaxy density skeleton" work (schema
  v6-v8: `sector_vertices`, `galaxy_shape`, `galaxy_shell_band`,
  `galaxyGen.ensure_sector_generated` lazy generation -- see [5.4.7] and
  [5.4.8] below) into the MySQL backend from [5.5.0].** That work landed
  on `main` entirely against the pre-port SQLite backend while this
  branch's MySQL port was in flight from the same base commit, so every
  new table, `_db.py` function (`get_sector_id_at`, `save_galaxy_shape`,
  `get_galaxy_shape`, `replace_galaxy_shell_bands`,
  `get_galaxy_shell_bands`, the vertex-writing/-reading in
  `insert_sector`/`get_sector_galaxy_position`), and CLI script
  (`galaxyGen.py`'s rewrite, the new `galaxyPlan.py`) needed porting to
  MySQL/`pymysql`/`MySQLConfig` as part of reconciling the two lines of
  work, same conventions as the rest of the MySQL port: `?` placeholders
  (translated by the `Connection` wrapper), `--mysql-*`/`config=` instead
  of `--db-path`/`db_path=`, and MySQL's `INSERT ... ON DUPLICATE KEY
  UPDATE` instead of SQLite's `INSERT ... ON CONFLICT DO UPDATE` for
  `save_galaxy_shape`'s singleton-row upsert.

### Fixed
- **`stellarObjects/galaxyGeometry.py` called an undefined
  `_theta_for_index` in both `sector_position_pc` and
  `sector_wedge_vertices_pc`** -- a pre-existing bug on `main` (confirmed
  present there independent of this merge), only `_phi_for_index` had
  ever been defined despite the module's own `GOLDEN_RATIO` docstring
  describing the golden-angle azimuthal step it was supposed to compute.
  Every galaxy-placement test failed with `NameError` until this merge's
  full test run surfaced it. Added the missing function
  (`theta_i = (2*pi*i / GOLDEN_RATIO) mod 2*pi`), matching the exact
  formula `test_galaxy_geometry.py`'s own worked-example test already
  documented and asserted against.
- **`stellarObjects._db.Connection`'s `execute`/`executemany` could
  return a `tuple` instead of a `list` from `.fetchall()` when zero rows
  matched** -- confirmed by testing, `pymysql`'s own cursor returns `()`
  for no rows but a `list` when rows exist, unlike `sqlite3`'s cursor
  (always a `list`). Surfaced by a real MariaDB test run as a spurious
  `assert result == []` failure (`test_ensure_sector_generated_reports_no_content_outside_every_stored_band`)
  that never appeared without a live server. Added a `_Cursor` proxy
  wrapping every returned cursor to normalize `fetchall()` to always
  return a `list`, so no other call site (present or future) can trip
  over the same inconsistency.

## [5.6.0] - 2026-09-09

### Added
- **Rate limiting (Flask-Limiter), applied to the whole API out of the
  box.** Default limits (200/day, 50/hour per client IP -- Flask-Limiter's
  own quickstart example, overridable via `PLANETGEN_RATELIMIT_DEFAULT`)
  apply to every route except `/api/health`; every write endpoint (below)
  additionally layers a stricter 10/minute limit on top
  (`routes.WRITE_RATE_LIMIT`). Exceeding a limit returns `429
  {"error": "rate limit exceeded", "detail": "..."}` (never Flask-Limiter's
  own default plain-text body) with `Retry-After`/`X-RateLimit-*`
  headers. Storage backend defaults to in-memory
  (`PLANETGEN_RATELIMIT_STORAGE_URI`, correct for a single-process
  deployment; a multi-worker `mod_wsgi`/`gunicorn` deployment needs a
  shared backend, e.g. Redis, or each worker under-enforces the
  configured limit by tracking its own separate counters).
- **Write-stub endpoints**: `POST`/`PATCH`/`DELETE` on `/api/sectors` and
  `/api/systems`. Every one validates its JSON request body (sectors get
  a fully mapped-out `{"name", "edge_ly"}` schema; systems only check
  "is a JSON object" for now -- see docs/api.md for why that one's
  field-level schema is still an open design question) and applies
  `WRITE_RATE_LIMIT`, but always responds `501 {"error": "... is not
  implemented yet"}` -- no row is ever inserted, updated, or deleted.
  Exists now so the request/response contract is settled and testable
  before the real database logic lands; see docs/api.md's "Write
  endpoints" section for what filling them in for real will also need
  (a write-capable database account, and authentication/authorization --
  neither exists yet, and both are load-bearing on this API staying
  read-only in practice today).

## [5.5.0] - 2026-09-09

### Added
- **Flask API: pagination, input validation, health check, and JSON-only
  error handling.** `/api/sectors`/`/api/systems` now return a paginated
  `{"items", "total", "limit", "offset"}` envelope instead of a bare list
  (`limit` defaults to 100, clamped to 500; `offset` defaults to 0) --
  this project's own roadmap plans galaxy-scale generation, so an
  unbounded listing endpoint would eventually return an unbounded
  response. Every query parameter (`limit`/`offset`/`sector_id`/`radius`)
  is now validated and rejected with a `400 {"error": "..."}` rather than
  silently ignored or crashing. Added `GET /api/health` for
  liveness/readiness monitoring. Every error response -- 400, 404
  (including an unmatched route), 405, and 500 -- is now JSON, never
  Flask's default HTML error page, and a 500 never leaks exception detail
  to the client. See `docs/api.md`.
- **MySQL backend (`TODO.md` Phase 5), replacing SQLite entirely.**
  `stellarObjects/schema.sql` is now MySQL/InnoDB DDL (`BIGINT UNSIGNED`
  ids, `DOUBLE`/`VARCHAR`/`TEXT`/`LONGTEXT` typing, a real
  `schema_migrations` tracking table replacing `PRAGMA user_version`,
  every index/foreign key declared inline per table for
  `CREATE TABLE IF NOT EXISTS` idempotency). `stellarObjects/_db.py` now
  talks to MySQL via `pymysql` (pure-Python driver) through a small
  `Connection` wrapper that keeps every existing call site's
  `conn.execute(sql, params)`/`row["column"]` shape unchanged, backed by
  a real connection pool (`DBUtils.PooledDB`) per TODO.md's "add real
  connection pooling" note. Every tool that touches the database
  (`sectorGen.py`, `systemGen.py`, `galaxyGen.py`, `queryDb.py`,
  `migrateDb.py`, `src/api/`, `src/html/`) now takes `--mysql-*` flags/
  `PLANETGEN_MYSQL_*` environment variables (`stellarObjects._db.MySQLConfig`)
  instead of `--db-path`/`PLANETGEN_DB_PATH`; the CGI browser's database
  picker now lists MySQL schemas on the configured server (filtered by
  `PLANETGEN_MYSQL_DATABASE_PREFIX`) instead of `.db` files in a
  directory. A new one-time `src/migrateSqliteToMysql.py` script imports
  an existing pre-port SQLite database (already at schema v5) into
  MySQL. The SQLite-specific `v1`-`v5` in-place migration machinery
  (`migrate_database`'s per-version steps, gzip file backups,
  `BACKUP_MARKER`) is removed entirely, since every MySQL database this
  project creates now starts at the current schema directly. See
  `docs/database-schema.md`.

### Changed
- `setup.py`'s `install_requires` gained `pymysql`/`DBUtils` (core
  dependencies now, not just the `api` extra) -- every database-touching
  entry point needs them, not only the Flask API.

### Fixed
- **`test_mdconvert.py`'s own `sys.path` setup pointed at a
  nonexistent top-level `html/lib/` instead of `src/html/lib/`,
  silently masked by `test_db_migration.py` (collected first,
  alphabetically) inserting the correct path first.** Deleting
  `test_db_migration.py` (see above) exposed it; fixed the path
  computation directly.

## [5.4.8] - 2026-09-09

### Added
- **Galaxy-wide density "skeleton"** (schema v8: `galaxy_shape`,
  `galaxy_shell_band`; `stellarObjects/galaxySkeleton.py`; `galaxyPlan.py`):
  precomputes and persists where the galaxy has any content at all,
  without storing a single sector's position/density/vertices -- those are
  pure deterministic functions of `(shell_index, shell_slot_index)` plus a
  handful of galaxy-wide shape parameters, so they're recomputed on demand
  instead. `galaxy_shape` holds the galaxy's shape parameters as one
  singleton row; `galaxy_shell_band` holds one row per contiguous
  *candidate* slot-index band per shell (almost always exactly one, found
  via an exact upper bound over spiral-arm azimuth, bisected to precision)
  -- a safe superset, not a per-sector list. Reduces the galaxy's
  structural storage from the ~1 TB a naive per-sector plan table would
  need down to ~350 KB at real Milky-Way scale (~4,076 rows), built in
  parallel across shells (`multiprocessing.Pool`) in under a second. See
  `docs/design/galaxy-coordinate-system.md` section 10 for the full
  analysis and measurements.
- **`galaxyGen.ensure_sector_generated(shell_index, shell_slot_index)`**:
  the lazy, visit-triggered generation entry point built on the skeleton
  above -- returns an already-generated sector if one exists; otherwise
  checks the stored band and exact density to decide whether the address
  holds anything, and if so generates and persists it on the spot, using
  that position's own `relative_density` as the actual system-count
  multiplier (so a bulge sector and a sparse outer-disk sector generate
  proportionally different counts, not a uniform default). A new
  `sectors.UNIQUE (shell_index, shell_slot_index)` constraint (schema v8)
  turns a concurrent visit race into a recoverable `IntegrityError`
  instead of a duplicate row.
- **Galaxy disk/spiral density model implemented** (`stellarObjects/galaxyDensity.py`):
  exponential disk radial falloff x sech^2 vertical scale-height x
  logarithmic spiral-arm modulation, plus a spherical bulge, normalized so
  `relative_density == 1.0` at a calibration point. `build_galaxy_shape`
  auto-calibrates the normalization constant; `predicted_star_count` scales
  a caller-supplied baseline expected system count by relative density.
  Tested in `tests/test_galaxy_density.py`, and exercised end-to-end
  alongside sector geometry in a full (unfilled) small spiral galaxy
  simulation.

### Fixed
- **`galaxyDensity.relative_density` could raise `OverflowError` far off
  the galactic plane** relative to `disk_scale_height_pc`: its vertical
  falloff term computed `1.0 / math.cosh(x) ** 2` directly, which raises
  once `|x|` exceeds ~710 -- found by a new skeleton test using a
  toy-scale shape, reachable in production for any shape with a small
  scale height relative to its own radius. Fixed with an algebraically
  equivalent, overflow-safe rewrite (`_sech_squared`) that underflows to
  the correct `0.0` limit instead of crashing.
- **Two same-shell Voronoi tessellation bugs found by that end-to-end
  simulation, both specific to small/sparse shells** (`sectorGeometry.py`):
  (1) the same-shell candidate projection used the raw chord vector
  instead of a proper gnomonic (central) projection, understating real
  separation worse the farther away a candidate was; (2) the half-plane
  test applied after that projection used the flat-plane bisector formula,
  which only approximates the true spherical bisector in the small-angle
  limit and was too permissive on sparse shells, silently under-clipping
  cells. Both fixes are covered by new regression tests parametrized to
  include a small shell (shell 1), plus an exhaustive
  every-vertex-is-shared check for a fully-populated small shell. See
  `docs/design/galaxy-coordinate-system.md` section 9 for the full
  derivation. Re-running the simulation after both fixes: 0 unexplained
  gaps across 31,255 generated outer vertices (previously 4,798).

## [5.4.7] - 2026-09-09

### Added
- **Every galaxy-placed sector now has explicit vertices with genuinely
  zero gaps against its same-shell neighbors, stored as plain relational
  rows.** New `sector_vertices` table (schema v7) -- one row per vertex,
  no JSON or other serialized blob anywhere in the schema (matching this
  project's existing convention: variable-length structured lists always
  get their own child table, the same treatment `asteroid_belt_composition`/
  `planet_reflection_spectrum` already have). Each sector's vertices are
  built from an exact local spherical Voronoi tessellation among its
  same-shell neighbors (`stellarObjects/sectorGeometry.local_lateral_cell`)
  extruded radially between the shell's inner and outer bounding spheres
  (`prism_vertices`). Lateral sharing is exact, not approximate: each
  shared corner is the 3D circumcenter of a sector and two of its
  neighbors -- a plain geometric fact independent of which of the three
  computes it, so two real neighbors land on identical floating-point
  values (~1e-16 agreement, verified directly) rather than two nudged
  approximations. Vertex/face count varies per sector (typically 5-7, mean
  6.0) because it has to: a cube tiling of a sphere can't be gap-free in
  general (the same reason a soccer ball needs pentagons mixed with
  hexagons), so a fixed 8-vertex shape cannot exactly reconcile a sector
  with more real neighbors than it has faces -- confirmed by an earlier,
  never-released attempt at exactly that (fixed-cube corner averaging),
  which only ever reduced gaps (~35% aggregate), not eliminated them.
  Radially, adjacent shells match in total area covered, not
  vertex-for-vertex (a "non-conforming mesh interface," the same technique
  used where independently-meshed regions meet in finite-element/CFD
  meshing) -- shell k's outer bound and shell k+1's inner bound are the
  same sphere, and each shell's own sectors fully and independently tile
  it. Getting this fast required exploiting that the underlying placement
  is a Fibonacci lattice: true same-shell neighbors concentrate at
  Fibonacci-number index offsets (confirmed against real data: exact
  offsets of F_20 through F_23), which collapses same-shell neighbor
  search from ~0.4-0.5s/sector (the naive radius-search cost near a large
  outer shell's equator, exactly where the disk/spiral density model
  concentrates real generation) down to ~0.3-0.4ms/sector -- validated
  against a guaranteed-correct brute-force search across 840 cases
  spanning the full polar range and shell sizes from 3 to 211 million
  slots with zero mismatches, with small shells
  (`shell_sector_count(k) <= 2000`) falling back to unconditionally-correct
  brute force since the underlying asymptotic theory isn't reliable that
  close to the galactic core. `galaxyGen.py` computes and stores this for
  every sector it generates; `sectorGen.py`'s standalone (non-galaxy) CLI
  is unaffected, same as every other galaxy-frame data. New
  `_migrate_v6_to_v7` schema migration drops any v6 database's briefly-lived
  `vertices_pc` JSON column entirely rather than converting it into rows
  (a sector's vertices are cheap to recompute from its address if ever
  actually needed) -- `_migrate_v4_to_v5`/`_migrate_v5_to_v6` need no
  special-casing at all, since a v4 or v5 source's `sectors` table already
  has the same shape as the current one.

## [5.4.6] - 2026-09-07

### Fixed
- **System page table of contents was simply unavailable below a 90rem
  window width.** It lived exclusively in a fixed right-hand margin rail
  (`.toc` in `src/html/static/style.css`) that only existed at
  `min-width: 90rem`; narrower windows got `display: none` and no
  alternative. It's now a collapsible pulldown rendered inline above the
  description at any width, expanding on click, and only switches to the
  fixed sidebar (always shown open) once the window is wide enough for
  one. Built with a hidden checkbox + `<label>` rather than native
  `<details>`/`<summary>`: a closed `<details>`'s children turn out to
  stay out of layout/paint via internal browser state that isn't fully
  reachable through CSS overrides of `display`/`content-visibility`
  (confirmed by testing -- only toggling the element's `.open` property
  itself, not any stylesheet rule, restored it), which made "always open
  at the wide breakpoint" unreliable. A plain checkbox has no such
  internal state to fight.

## [5.4.5] - 2026-09-07

### Fixed
- **Sector Map: the 3D scene drew over surrounding page content instead
  of staying confined to its panel.** Zooming in (or rotating to an angle
  where the cube's diagonal grew past its nominal footprint) had nothing
  bounding where the map was visible, so it grew and drew over the rest
  of the page. Added `.starmap-viewport`, a fixed-size (`aspect-ratio: 1
  / 1`, matching the square footprint the old static image occupied),
  `overflow: hidden` window that the 3D scene now rotates/zooms/pans
  inside of -- clipped at its edges instead of spilling out. The 3D
  content inside renders/rotates/occludes exactly the same either way;
  this only bounds where it's visible from the outside.

## [5.4.4] - 2026-09-07

### Added
- **Sector Map is now interactive 3D: drag to rotate, scroll (or +/-
  buttons) to zoom.** Replaced the static server-baked isometric SVG
  projection in `src/html/lib/starmap.py` with a real CSS 3D scene
  (`transform-style: preserve-3d`, orthographic -- no `perspective`), so
  the browser's own compositor handles rotation and occlusion instead of
  a hand-rolled JS matrix routine; `src/html/static/sectormap.js` tracks
  two rotation angles and a zoom factor and feeds them straight to the
  scene's CSS transform. Each star dot is billboarded (counter-rotated
  every frame to keep facing the camera) so it stays a circle instead of
  going edge-on as the view turns, and dot-click detection is resolved by
  geometry (`getBoundingClientRect`) rather than native hit-testing,
  since the latter turns out unreliable for elements nested this deep in
  a rotated `preserve-3d` hierarchy.

### Changed
- **Sector Map star colors now come from the star's actual named
  spectral color** (`SPECTRAL_CLASS_COLORS`: Blue/Blue-White/White/
  Yellow-White/Yellow/Orange/Red, keyed off `star_type`'s leading letter)
  instead of a raw Kelvin-to-RGB blackbody approximation, so a "White
  Giant" reads white and a "Blue Giant" reads blue regardless of its
  exact temperature. Luminosity shades each color's vividness/lightness
  (brighter = more vivid, dimmer = more muted) and temperature nudges
  lightness within the star's own spectral band.
- **Binary systems draw two dots** (a larger primary and a smaller
  secondary, capped at 65% of the primary's radius and offset to its
  lower-right, overlapping) built from each component's own `stars` row,
  instead of one dot from the system-level `binary_type`/
  `binary_temperature_k` summary -- so each half of a binary is sized and
  colored from its own actual data.

See `docs/TODO.md` ("Near-term: interim `../src/html/` browser
enhancements") for where this started as a plan.

## [5.4.3] - 2026-09-07

### Changed
- **Web interface: system-page table of contents moved to a fixed
  right-margin rail.** The collapsed `<details>` TOC added in [5.4.2]
  still lived inline next to the description as a flex sibling, narrowing
  the prose column whenever it was open. `system.py`'s `_toc_html` now
  renders a plain, always-expanded `<nav>`, and `style.css` positions
  `.toc` fixed in the right margin (mirroring the left `.sidenav`),
  appearing only once the window is wide enough (`min-width: 90rem`) to
  hold it without crowding the main content -- narrower windows simply
  don't get one, rather than it floating over the page.

## [5.4.2] - 2026-09-07

### Removed
- **Class R ("an ejected, geologically active world") removed entirely.**
  It had `h`/`e`/`c` all `False` -- zero probability weight, unreachable
  outside a manual `zone_override`. A genuinely free-floating/rogue planet
  (no host star at all) is a real exoplanet category, but doesn't fit this
  generator's star-centric zone model -- "which zone" is the wrong
  question for an object with no star to be zoned relative to. Rather than
  leave it as permanent dead weight, cut from `PLANET_CLASSES`,
  `PLANET_CLASS_PROBABILITIES` (already 0.0000), and `MOON_BLACKLIST`
  (`program_constants.py`); the now-vacuous
  `test_known_issue_class_with_no_valid_zone_is_unreachable` regression
  test removed from `src/tests/test_planets.py`.

### Fixed
- **Naming: triple-consonant validation bug.** `utils.is_name_valid`
  enforced "no more than two consecutive vowels" but had no matching check
  for consonants, so a fixed denylist of specific clusters
  (`BAD_CONSONANTS`) was the only thing standing between a generated name
  and an arbitrary run of 3+ consonants -- empirically ~40% of generated
  star names had one. Added a `consonant_count` run-length counter
  symmetric to the existing `vowel_count` one. Verified 0/5000 across
  star/planet/moon/sector name generation post-fix.
- **Web interface: navbar had no way back to the database picker on a
  single-database deployment.** `lib/page.py`'s sidenav only linked to
  `index.py` when more than one `.db` file existed, since `index.py`
  itself auto-redirects past the picker for a single database -- but that
  meant the common single-database production deployment
  (`starmap.moltenaether.com`) had no Databases link at all. Fixed by
  always linking to `index.py?all=1`, and teaching `index.py` to render
  the full picker table (bypassing its own auto-redirect) when `?all=1`
  is present.
- **Web interface: system-page table of contents wasn't collapsible.**
  `system.py`'s `_toc_html` now builds a `<details>`/`<summary>` pair
  instead of a plain `<aside>`/`<h3>`, collapsed by default -- a native,
  no-JS toggle.
- Minor cleanups found in passing: an unused `program_constants` import in
  `sectorGen.py`, and a stray `f`-string prefix on a non-interpolated
  string in `evolution.py`.

## [5.4.1] - 2026-09-07

### Added
- **`sectorGen.py --density`: a controllable density value for sector
  generation.** Sectors previously always generated a flat count of systems
  (`--num-systems`, default 10) with no connection to sector volume -- there
  was no way to make one generated sector meaningfully denser or sparser
  than another. `--density` is a float multiplier on the real local stellar
  density this codebase already models (`physical_constants.LOCAL_STELLAR_DENSITY_LY3`
  via `SpaceSector.expected_system_count()`) -- 1.0 means a realistic sector
  this size, 2.0 twice as dense, 0.5 half. The actual per-sector count is
  drawn with the existing `_sample_poisson_count` helper (previously only
  used by `SpaceSector.grow_from_seed`), resolved fresh for each sector
  rather than once at parse time, so counts vary naturally sector to sector
  under `--num-sectors` and across `galaxyGen.py` runs, matching a real
  spatial Poisson process. Mutually exclusive with `--num-systems`; neither
  flag given keeps the original flat-10-systems default unchanged.
  `galaxyGen.py` gets the flag for free since it shares
  `sectorGen.add_shared_generation_options`/`validate_shared_generation_args`
  with no changes needed on its side.

## [5.4.0] - 2026-09-07

### Removed
- **Two brown-dwarf-scale gas-giant classes cut entirely.** Correcting
  their radius ranges to real sub-stellar physics in [5.3.9] left them as
  near-duplicates of each other (overlapping radius, differentiated only by
  an invisible density number no description text reflects), each carrying
  a vanishingly small (0.01%) generation weight, representing objects that
  aren't really planets in the first place (brown dwarfs are sub-stellar
  objects; this codebase already has a dedicated mechanism for stellar/
  sub-stellar companions via `BinaryStarProxy`/`BINARY_SYSTEM`). Their
  combined generation weight folded into the ordinary Jupiter/Saturn-class
  gas giant. All remaining references to the removed classes (and to the
  earlier-removed small-rocky-class pair) scrubbed from comments/docs,
  generically describing what they explain rather than naming
  now-nonexistent class letters.

### Added
- **Per-class `size_mode`: bell-curve size distributions for planets and
  moons.** Radius generation (`planetPhysics.generate_planet_properties`'s
  three radius-draw sites, and `generate_moons`' moon-radius draw) switched
  from a flat uniform draw across each class's declared radius range to a
  bounded Gaussian ("bell curve") draw peaking at a new per-class
  `size_mode` value (0.0-1.0, "what fraction through the available range is
  the statistically most common size") via new `utils.sample_bounded_bell`.
  For a moon, "available range" is its actual Hill-sphere/mass-capped
  window, not necessarily the class's full declared range -- `size_mode` is
  read relative to whatever range is actually being drawn from. Every
  surviving class's `size_mode` is set from a real-world single-body analog
  where one exists (Earth for M, Mars for K, Venus for N, Jupiter/Saturn
  for J, Uranus/Neptune for I/T, the real rocky-to-gaseous transition
  radius for V) or a reasoned default (small-body populations, both real
  asteroids/KBOs (Class C/D) and this generator's own hot-zone rocky
  classes (A/B), skew toward their smaller end, matching real
  size-frequency distributions). Verified empirically: generated radius
  means land within ~1% of each class's real-world anchor (e.g. Class M
  mean 6,337km vs Earth's 6,371km; Class K mean 3,410km vs Mars's 3,389.5km).
  The self-adjusting spread (`sample_bounded_bell`'s `spread_divisor`)
  keeps a legible bell shape even for a mode pinned near one edge of a
  range, via rejection sampling rather than clamping (which would pile
  spillover probability mass at the boundary).
- `PLANET_CLASS_PROBABILITIES`'s `J` weight increased 0.0529 -> 0.0531 to
  absorb the removed brown-dwarf-scale classes' combined weight.
- Atmospheric-pressure sanity bound (`test_full_matrix.py`/`test_planets.py`)
  lowered back down `1e9` -> `5e7` Pa now that the brown-dwarf-scale
  gravity/pressure extremes those classes produced are gone (observed max
  across the full star-type matrix is now ~1.1e7 Pa, from Class N).
- `test_planet_physics_fixes.py`'s gravity/pressure correlation test
  threshold recalibrated (Spearman `>0.5` -> `>0.1`) for the narrower
  gravity range the remaining gas giants span without the removed classes'
  two-orders-of-magnitude density spread; its density-blend-skip test
  rewritten to inject a temporary `density_range` onto an existing class
  via `monkeypatch` rather than depend on a specific class declaring one
  (no class currently does -- it's generic, reusable override
  infrastructure, same as every other per-class override this codebase has
  built up).

## [5.3.9] - 2026-09-07

### Removed
- **Classes X and Y removed, merged into B and A/B respectively.** Both
  were small, redundant variants of existing hot-zone rocky classes:
  - Class X ("a stripped core from a gas giant", no atmosphere, radius
    500-5000km) folded into Class B ("a small, molten world") as an
    alternate origin story ("occasionally the stripped core of a former
    gas giant") rather than a separate atmosphere-less class -- a
    Mercury-analog this close to its star already has only a negligible
    exosphere in reality, so B's existing thin atmosphere covers X's "no
    atmosphere" identity closely enough.
  - Class Y ("a 'demon' class world", radius 5000-7500km) folded into
    *both* Class A and Class B (both radius ceilings raised 5000->7500km
    to absorb Y's size range), its toxic/irradiated flavor folded into
    each class's description as a variant rather than kept as a third
    class. B's composition gained "and sulfur" for the shared
    volcanic/irradiated theme.
  - `MOON_BLACKLIST` and `PLANET_CLASS_PROBABILITIES` updated to drop X/Y
    (their combined generation weight folded into A/B rather than dropped).

### Changed
- **Gas-giant zones reworked using real exoplanet science.** Previously
  every gas/ice-giant class (I/J/S/T/U) was valid in zone `c` only --
  meaning no gas giant could ever appear close to its star or in the
  habitable zone, despite "hot Jupiters" and "warm Jupiters" being a
  standard, well-documented real classification (orbital period < 10 days
  / 10-365 days respectively) and Neptune-mass planets in or near a
  temperate zone being common too.
  - Class J now valid in `h`/`e`/`c` (hot/warm/cold Jupiter) -- a warm/cold
    Jupiter placed in zone `e` can also generate ordinary moons via the
    existing moon-generation path, including habitable-class ones (the
    "habitable exomoon around a giant planet" trope), verified working.
  - Class I now valid in `e`/`c` (warm/cold Neptune), deliberately **not**
    `h`: real close-in Neptune-mass planets are rare (the observed "hot
    Neptune desert") because stellar irradiation photoevaporates a
    Neptune-mass H/He envelope down to a bare rocky/metal core well before
    it could stay Class I -- that outcome is exactly Class B's newly-merged
    "stripped core" identity (see Removed, above).
  - Classes S/T/U left `c`-only: real directly-imaged super-Jovian/brown-dwarf
    companions are predominantly found at wide separations (formation and
    detection-bias reasons), and close-in high-mass companions, while known,
    are much rarer (the "brown dwarf desert").
  - Fixing `generate_moons` (planetPhysics.py) to respect
    `HABITABLE_WORLD=False` the same way direct planet generation already
    does -- unreachable before this change (no gas giant could ever be in
    zone `e`), but a warm/cold Jupiter placed there could otherwise roll a
    habitable-class moon even in a system explicitly configured to disallow
    habitable worlds.
- **Classes S ("supergiant") and U ("ultragiant") radius ranges corrected
  to real sub-stellar physics** -- previously 250,000-50,000,000km and
  25,000,000-60,000,000km respectively, i.e. up to ~86 solar radii, larger
  than most actual stars. Real brown dwarfs stay within ~15% of Jupiter's
  own radius (69,911km) across their *entire* 13-80 Jupiter-mass range
  (electron degeneracy pressure means more mass compresses them, R ~
  M^(-1/8)); even the smallest true hydrogen-fusing red dwarf stars are
  only ~0.1 solar radii (~69,600km, essentially Jupiter-sized). Corrected to
  S: 60,000-120,000km, U: 65,000-130,000km (deliberately overlapping --
  that overlap *is* the real physics, not an oversight). What actually
  differentiates a higher-mass sub-stellar object at essentially the same
  radius is **density**: real measured brown-dwarf densities run roughly
  10-200 g/cm^3 (~10-150x an ordinary gas giant's), so both classes gained
  a new `density_range` (S: 10.0-60.0, U: 60.0-150.0 g/cm^3). Class T ("gas
  dwarf") also corrected, 250,000-25,000,000km -> 15,000-55,000km -- its
  old floor matched S's own floor, letting a "dwarf" be exactly as large as
  a "supergiant"; now spans ice-giant-to-Saturn scale, meaningfully below
  both J and S.
- **New per-class `density_range` override** (`PLANET_CLASSES`, read by
  `planetPhysics.get_planet_mass_ranges`/`generate_planet_properties` and
  mirrored in `plausibility.theoretical_gravity_bounds_g`) -- same
  `.get(..., default)` pattern as `atm_molar_density_range` etc. Used by
  Classes S/U above.
- **Fixed a gas-giant density-blend bug that silently collapsed every gas
  giant's density to a near-zero, physically meaningless value**, previously
  flagged as necessary-but-deferred follow-up work in
  `test_gas_giant_sampled_densities_are_finite_positive_and_within_theoretical_bounds`'s
  own docstring. `generate_planet_properties`' core/envelope harmonic-mean
  density blend reused `planet.atm_density` (drawn from
  `ATMOSPHERE_DENSITY["g"]`) as the envelope term -- but that value is a
  thin, surface/pressure-layer density (a different physical layer, correct
  for the separate atmospheric-pressure/scale-height calculation), ~1000x
  lighter than a real gas-giant envelope's actual bulk density. A harmonic
  mean is dominated by whichever term is smaller almost regardless of mass
  fraction, so this collapsed density to ~0.001-0.003 g/cm^3 for every gas
  giant, independent of core density -- silently defeating the new S/U
  `density_range` work above (their much denser cores had almost no effect
  on the final blended density). Fixed by introducing
  `physical_constants.GAS_ENVELOPE_BULK_DENSITY` (0.06-0.3 g/cm^3, grounded
  in real measured "puffy" gas giants -- WASP-193b's ~0.06 g/cm^3 is the
  lowest confirmed bulk density known), used only for this blend. A class
  declaring its own `density_range` (S/U) now skips the blend entirely and
  uses that density directly -- real brown dwarfs don't have a meaningfully
  separate light envelope over a denser core the way an ordinary gas giant
  does, so blending toward a light "puffy" value would just dilute the
  elevated density right back down.
- Atmospheric-pressure sanity bound (`test_full_matrix.py`/`test_planets.py`)
  raised `5e7` -> `1e9` Pa -- Classes S/U's now-correct brown-dwarf-like
  gravity legitimately pushes pressure up to ~3.3e8 Pa across the full
  star-type matrix.
- `test_planet_physics_fixes.py`'s two gas-giant-blend tests updated to
  match the new formula (queuing an `envelope_density_gcm3` draw instead of
  reusing `atm_density`), plus a new
  `test_density_range_override_skips_the_blend` covering the S/U direct-
  density path.

## [5.3.8] - 2026-09-07

### Fixed
- **Habitable/life-bearing classes restricted to the ecosphere zone.** Class
  Q (`h`/`e`/`c` all `True`) and Class W (`h`/`e` `True`) both carry a
  `life_chemical`, but weren't restricted to zone `e` like every other
  life-bearing class already was — meaning a life-bearing Q or W world could
  be generated directly in the hot or cold zone. Both are now `e`-only
  (`h`/`c` `False`); Q's "eccentric orbit" flavor still holds fully confined
  to the ecosphere zone, and W's "tidally locked" identity is arguably more
  scientifically apt restricted to the ecosphere zone (real tidally-locked
  *habitable* worlds are a real, actively studied trope specifically because
  a cool star's habitable zone sits close enough in for tidal locking to be
  near-guaranteed — e.g. TRAPPIST-1's planets, Proxima b).
  New regression test `test_life_bearing_classes_are_ecosphere_only`
  (`src/tests/test_planets.py`) locks this in for every class with a
  `life_chemical`, not just the ones on the separately-curated
  `HABITABLE_PLANET_CLASSES` list.
- **Description/atmosphere-text cleanup across `PLANET_CLASSES`.** Several
  classes' `description` field repeated a word the render template
  (`planetData.to_paragraph_list`) already supplies via its own "with an
  atmosphere of {atmosphere}" or "with a composition of {composition}"
  clause, producing genuinely broken rendered sentences — e.g. Class B
  previously rendered "...a small, molten world **with a thin atmosphere
  with an atmosphere of** a mix of helium, sodium, and oxygen...". Fixed for
  Classes B, E, K, N, Q, W, X, and Y — descriptive detail that belonged on
  the atmosphere field (e.g. "thin", "dense, reducing") was moved there
  instead of dropped, and Y's atmosphere field (previously the only one not
  ending in a named gas mixture) reworded to "a turbulent, toxic, and
  irradiated mix of gases".

## [5.3.7] - 2026-09-07

### Changed
- **Greenhouse-formula fix and per-class climate tuning.** The prior
  `greenhouse_factor` formula (`planetPhysics.calculate_atmospheric_conditions`)
  used `atm_molar_density` as its only lever, scaled by a single shared
  `CO2_MAX_GREENHOUSE_FACTOR` cap (5) — real Earth air's own molar mass
  already produced `greenhouse_factor ≈ 3.33` under that formula, driving
  Class M's mean surface temperature to 362K instead of ~288K, and Mars and
  Venus (nearly identical real molar mass, ~100x different real greenhouse
  forcing) could never be told apart by molar mass alone.
  - `CO2_MAX_GREENHOUSE_FACTOR` (now 500) is a generous safety ceiling, not
    the calibration knob.
  - New per-class `PLANET_CLASSES` keys — `albedo_range` (extended beyond
    Class P, which introduced it), `atm_molar_density_range`,
    `atm_density_range`, and `greenhouse_multiplier_range` — give every
    tuned class independent control of composition (molar density),
    quantity (mass density), and potency (greenhouse multiplier), following
    the exact override pattern `albedo_range` established for Class P.
    Class N's old hardcoded `atm_density = 65` / `atm_molar_density = max`
    special case in `planetPhysics.generate_planet_properties` is folded
    into this same general mechanism.
  - **Tuned this pass** (via the new `climate_tuning_cli.py`, see below):
    M (Earth analog, ~286K/~99kPa), O (warm/wet ocean world, ~293K/~97kPa),
    H (hot/dry desert, ~325K/~42kPa), K (Mars analog, ~231K/~0.57kPa), L (K
    + vegetation, warmer/thicker than K, ~256K/~2.2kPa), N (Venus analog,
    ~740K/~9.4MPa), E/F/G (a young, cooling progression, ~373K -> ~329K ->
    ~292K, converging near M), and V (thick, hot Super-Earth, ~365K/~286kPa).
    Class P and W are unchanged this pass (P already had a working
    `albedo_range`; W's "extreme temperature variations" identity needs a
    day/night model this generator doesn't have, not just range tuning).
  - Text refinements: Class H's "and metals" -> "and mineral dust" (a real
    desert's atmosphere lofts particulate, not metal vapor); Class O's
    atmosphere text (previously byte-identical to Class M's) now mentions
    water vapor; Class E's vague "hydrogen compounds" now names a real
    Hadean/Archean-analog mix (water vapor, ammonia, methane); Class K's
    "carbon dioxide" -> "a thin mix of carbon dioxide and nitrogen" (names
    the real Mars-analog composition, not just class); Class V's atmosphere
    text now reflects its tuned CO2-retention (not primordial H/He) identity.
  - Class M's disabled gravity clamp (`planetPhysics.calculate_surface_gravity`)
    deleted outright (it had been commented out, inert, since the
    atmospheric-pressure formula fix; no longer worth keeping around "in
    case it needs restoring").
- **New: `src/tests/climate_tuning_cli.py`.** A human-driven tuning tool
  (mirrors `physical_plausibility_cli.py`'s batch-generate-and-report
  pattern, reusing `stellarObjects.plausibility`'s engine): generates N
  bodies of one class, reports summary stats, and — for classes with a
  direct real-world analog (M/Earth, K/Mars, N/Venus) — a delta line against
  that reference. `--albedo`/`--molar-density`/`--density`/`--greenhouse`
  flags temporarily monkeypatch that class's `PLANET_CLASSES` entry for the
  run only, so candidate values can be iterated without editing source
  between runs.
- **New: `src/tests/test_climate_tuning.py`.** Regression suite locking in
  the tuning above via generously-toleranced bands (Class M/N/K within
  Earth/Venus/Mars-like ranges) and relative orderings (N hottest/
  highest-pressure of the tuned classes; K colder/thinner than L; L colder
  than M; H hotter/drier than O; E > F > G cooling progression converging
  near M; M < V < N) rather than brittle exact-value assertions, since
  generation is inherently stochastic.
- `test_full_matrix.py`/`test_planets.py`'s atmospheric-pressure sanity
  bound raised from `1e7` to `5e7` Pa — Class N now legitimately reaches
  ~9-15MPa depending on host star luminosity (previously ~2.98MPa mean, per
  the greenhouse-formula bug above), and the old bound was sized for the
  un-tuned, incorrectly-cold N. `test_planet_physics_fixes.py`'s
  `test_class_p_has_own_albedo_range_distinct_from_default` updated to
  check P's range differs from M's own new tuned range, rather than
  asserting M has no override at all.

## [5.3.6] - 2026-09-07

### Added
- **Galaxy-scale coordinate system (Track C), merged.** Sectors can now
  be placed on a galaxy-wide, radial shell/Fibonacci-sphere tiling
  (`docs/design/galaxy-coordinate-system.md` sections 0-8) instead of
  existing only in isolation:
  - **Schema v3 -> v4** (`stellarObjects/schema.sql`): `sectors` gains six
    nullable galaxy-frame columns (`center_x/y/z_pc`, `galactic_radius_pc`,
    `shell_index`, `shell_slot_index`), NULL together for a sector never
    placed in a galaxy. `_db.migrate_database`'s new `_migrate_v3_to_v4`
    handles the upgrade (existing `_migrate_v1_to_v2`/`_migrate_v2_to_v3`
    updated in lockstep, per that function's "each hop maps straight to
    the current schema" design); `docs/database-schema.md`'s schema
    history documents the change.
  - **`stellarObjects/galaxyGeometry.py`** (new): the shell/Fibonacci-sphere
    tiling primitives (`shell_sector_count`, `shell_radius_pc`,
    `sector_position_pc`) and `enumerate_sectors_within_radius` — an
    exact, two-prune generation-unit primitive that finds every sector
    address within a radius of an arbitrary galaxy-space point without
    ever scanning a whole shell's slot count (design doc section 8).
  - **`galaxyGen.py`** (new, repo root): a CLI generating many sectors as
    one galaxy, reusing `sectorGen.py`'s own per-sector generation/save
    path. `--shell K` batch-generates a whole radial shell (guarded by
    `LARGE_SHELL_WARNING_THRESHOLD`, needing `--limit`/`--yes` above it);
    `--center-sector ID --radius-pc R` generates a local neighborhood
    around an already galaxy-placed sector. Either mode skips slots a
    sector already occupies.
  - **`GALACTIC_CENTER_DISTANCE_LY` is now per-sector**, not a single
    fixed constant: `Star.calculate_system_perimeter`/
    `BinaryStarProxy._calculate_system_perimeter_static` accept a
    `galactic_center_dist_ly` override, threaded from `galaxyGen.py`
    through `StarSystem`/`Star`/`BinaryStarProxy`, falling back to the
    old fixed constant for unplaced/standalone sectors
    (`sectorGen.py`'s own CLI unaffected).
  - `stellarObjects/utils.py` gains `mpc_to_pc`/`pc_to_mpc` (exact) and
    `pc_to_ly`/`ly_to_pc` (display-string conversions), per the design
    doc's unit-choice section.
  - New end-to-end coverage: `src/tests/test_galaxy_gen.py` runs
    `galaxyGen.py`'s actual CLI entry point (both `--shell` batch mode
    and `--center-sector` local-neighborhood mode) against a real
    temporary database, asserting correct shell/slot addresses, correct
    stored positions, no duplicate slots, and that already-occupied
    slots are skipped on re-run — the gap `docs/TODO.md`'s Phase 4 entry
    flagged as not yet done when this was paused mid-session. Developed on
    a branch cut before Track A's completion and the TODO/FIXME-comment
    migration (5.3.5); merged into `main` after both, with no functional
    changes needed beyond a `test_db_migration.py` assertion that had to
    decompress the (now gzip-compressed, per Track B) v3->v4 migration
    backup the same way the existing v1->v2 backup test already did.

## [5.3.5] - 2026-09-07

### Fixed
- **Atmospheric pressure is no longer independent of gravity**
  (`stellarObjects/planetPhysics.py`): the barometric-formula pressure
  calculation algebraically canceled gravity out entirely (`atmospheric_pressure
  = atm_density * R * T / atm_molar_density`), so a Neptune-gravity gas giant
  and a Jupiter-gravity one produced the same pressure. A new
  `_atmosphere_retention_factor(gravity_g)` (linear, normalized to 1.0 at
  Earth gravity) now scales an *effective* atmospheric density used only in
  the pressure calculation (not `planet.atm_density` itself, which also
  feeds the gas-giant density blend), reintroducing a real, tunable
  gravity/pressure relationship (Spearman correlation on a mixed
  terrestrial/gas-giant sample now > 0.5, vs. ~-0.11 before).
- **Class P ("cold, glaciated") is colder than Class M again**: gave Class P
  its own `albedo_range` (0.5-0.7, matching real ice/snow Bond albedo)
  instead of sharing the default rocky/Earth-like range (0.12-0.35) with
  every other terrestrial class. Previously P and M were statistically
  indistinguishable in temperature once the disabled clamp was removed (see
  `docs/analysis/habitability-atmosphere-sanity-review.md`); P's cold
  identity now comes from the unclamped physics instead of a post-hoc
  override.
- Completes Track A (see [5.3.4]'s gas-giant density/greenhouse fixes) --
  8682 tests passing, including 10 new regression tests in
  `src/tests/test_planet_physics_fixes.py`.

## [5.3.4] - 2026-09-07

### Fixed
- **Gas-giant density blend** (`stellarObjects/planetPhysics.py`): the
  core/atmosphere blend used a mass fraction as an arithmetic-mean weight
  between two densities, which is dimensionally wrong and could produce
  gas giants as low as 0.026 g/cm^3. Replaced with the mass-weighted
  harmonic mean, the physically correct way to combine two component
  densities via a mass fraction. `plausibility.py`'s independently
  reimplemented copy of this formula (used to derive analytical
  hard-invariant gravity bounds) was updated in lockstep, including its
  docstring's justification for corner-evaluation (still valid: the new
  formula is monotonic in each argument, just no longer multilinear).
- **Inverted greenhouse factor** (`stellarObjects/planetPhysics.py`): the
  formula rewarded an atmosphere for being *far* from CO2's molar density
  rather than for actually containing more CO2 — backwards from physical
  reality. Now scales directly with `atm_molar_density`, the only
  atmosphere-composition signal the data model has today.

### Changed
- **Database schema-migration backups are now gzip-compressed** and
  excluded from the web database picker and from a subsequent migration
  run (previously a plain `.db` copy that a naive `*.db` glob would both
  surface in the picker and silently re-migrate on the next run). See
  `docs/TODO.md`'s "File Management" section for detail.

### In progress, not yet merged (see `docs/TODO.md` for exact state)
- Atmospheric pressure is still algebraically independent of gravity, and
  Class M/Class P remain statistically indistinguishable — both scoped
  and partially started, paused mid-session in worktree
  `agent-a8acb02b98bed5b8d`.
- The galaxy-scale coordinate system's 8 open design questions were
  decided this session, and a schema v3->v4 migration,
  `GALACTIC_CENTER_DISTANCE_LY` fix, and a first `galaxyGen.py` were
  written but paused uncommitted in worktree `agent-a36f801e275fb2b71`
  before a full test pass — not part of this release.

## [5.3.3] - 2026-09-07

### Changed
- **Markdown consolidated into `docs/`, renamed for clarity.** Only
  `README.md`/`LICENSE.md`/`CHANGELOG.md` remain at the repo root; every
  other README and loose doc moved into `docs/` with a descriptive name:
  `TODO.md` -> `docs/TODO.md`, `src/html/README.md` -> `docs/html-interface.md`,
  `db/README.md` -> `docs/database-schema.md`,
  `src/api/README.md` -> `docs/api.md`, `apache/README.md` ->
  `docs/apache-deployment.md`, `examples/EXAMPLES.md` ->
  `docs/example-systems.md`, `examples/JSON.md` ->
  `docs/system-file-format.md`. Every cross-reference between them (and
  from code/scripts) was updated to match; a couple of pre-existing
  broken/mismatched links were caught and fixed along the way
  (`docs/system-file-format.md`'s "full command-line reference" link text
  didn't match its own target; `src/html/lib/dbutil.py`'s docstring still
  pointed at the pre-5.3.2 `stellarObjects/` path instead of
  `src/stellarObjects/`).
- **`wsgi.py`, `queryDb.py`, `migrateDb.py` moved into `src/`**, alongside
  `stellarObjects`/`api`/`tests`, so only the `*Gen.py` scripts
  (`sectorGen.py`/`systemGen.py`) are visible at the repo root as CLI
  entry points. Since these three now sit as direct siblings of
  `stellarObjects`/`api` under `src/`, Python's own sys.path[0] (the
  running script's own directory) already makes those packages
  importable -- the sys.path shims 5.3.2 added to them are gone, no
  longer needed (unlike `sectorGen.py`/`systemGen.py`, which stay one
  directory further away at the repo root and keep theirs).
  `migrateDb.py`'s `DEFAULT_DB_DIR` needed an extra `os.path.dirname()`
  level to still resolve to the repo-root `db/`, one directory deeper
  than before.
- **`physicalPlausibility.py` moved to `src/tests/physical_plausibility_cli.py`**,
  matching that directory's naming scheme, without a `test_` prefix so
  pytest doesn't try to collect it as a test module (it's a human-facing
  report generator, not a pass/fail check -- `src/tests/test_physical_plausibility.py`
  remains the actual automated test).
- **`apache/` moved into `examples/apache/`** (its `README.md` moved to
  `docs/apache-deployment.md` per the markdown-consolidation rule above);
  **`examples/*.json` moved into `examples/systems/`**, so `examples/`
  now holds two clearly-separated subfolders (`apache/` deployment
  config, `systems/` recipe files) instead of a flat mix of both kinds of
  example content. `.gitattributes`' LF-pinning rule for
  `apache/*.sh` was updated to `examples/apache/*.sh` to keep matching
  the actual file.
- `setup.py` and `pytest.ini` were **evaluated for a move into `src/` and
  kept at the repo root** -- both are genuinely not possible without
  breaking things, verified empirically rather than assumed:
  `pip install -e .` with `setup.py` moved silently built a bogus,
  empty `UNKNOWN-0.0.0` package instead of erroring (pip's PEP 517
  build only looks for `setup.py`/a full `pyproject.toml` project table
  at the invocation root, and this repo's `pyproject.toml` only declares
  a build backend, no project metadata of its own to fall back on); with
  `pytest.ini` moved, `pytest`'s own config-file discovery (which only
  searches the invocation directory and its parents, never a
  subdirectory) silently fell back to `pyproject.toml` and ignored
  `testpaths`/`pythonpath` entirely -- tests still happened to pass
  either way (coincidentally, via `pytest`'s own unrelated `__init__.py`-walkup
  sys.path behavior and an editable install's global registration of the
  root-level `py_modules`), which is exactly the kind of silent,
  environment-dependent fragility not worth introducing on purpose.

## [5.3.2] - 2026-09-07

### Changed
- **Repo layout: `stellarObjects`/`api`/`tests` moved under a new `src/`
  directory** (`src/stellarObjects/`, `src/api/`, `src/tests/`), so only
  the top-level CLI entry points (`sectorGen.py`, `systemGen.py`,
  `queryDb.py`, `migrateDb.py`, `physicalPlausibility.py`, `wsgi.py`) are
  visible at the repo root. `setup.py` now declares an explicit
  `package_dir` per discovered package rather than a blanket
  `package_dir={'': 'src'}`, since that would have also redirected the
  root-level `py_modules` lookups (`systemGen`/`sectorGen`) into `src/`,
  where they don't live. Every root entry script gained a small
  `sys.path` shim (inserting `src/` before its `stellarObjects`/`api`
  imports) so they keep working without requiring `pip install .` first,
  matching the no-install fallback `src/html/`'s CGI scripts already relied
  on -- those fallbacks (`src/html/lib/dbutil.py`, `src/html/sector.py`,
  `src/html/search.py`) were updated the same way. Two internal
  repo-root-relative path computations (`stellarObjects/webconfig.py`'s
  `_PROJECT_ROOT`, `stellarObjects/_db.py`'s `DEFAULT_DB_PATH`) needed an
  extra `os.path.dirname()` level to still resolve correctly one
  directory deeper; `pytest.ini` gained an explicit `pythonpath = . src`
  so both the entry scripts and the moved packages resolve during tests
  regardless of pytest's own import-mode heuristics.
- **`webconfig.json.example` moved into `src/html/`**; the real, gitignored
  `webconfig.json` stays at the repo root, outside Apache's `src/html/`
  `DocumentRoot`, for the same security reason `db/` already lives there
  (see `docs/webconfig.md`).
- **`WEBCONFIG.md` moved into `docs/`**, alongside this session's
  `docs/design/`/`docs/analysis/` additions, consolidating loose
  documentation in one place (`README.md`/`LICENSE.md`/`TODO.md`/
  `CHANGELOG.md` stay at the repo root, and per-directory READMEs
  `src/html/README.md`/`apache/README.md`/`db/README.md`/`src/api/README.md`
  stay next to the code they document).

## [5.3.1] - 2026-09-06

### Fixed
- **Evolved-star mass sampling could imply a pre-Big-Bang star.** An
  evolved-class star (`Yerkes != V`, e.g. a giant or supergiant) derives
  its required main-sequence lifespan from its own already-generated mass
  (`Star._calculate_initial_star_age_and_lifespan`'s evolved-star branch);
  for a sub-solar-mass progenitor (roughly under ~0.88 Msun), that implied
  lifespan alone already exceeded `UNIVERSE_AGE_GY` (13.8 Gy) -- meaning
  such a star couldn't actually have finished its main-sequence phase yet
  in the real universe. The [5.3.0] universe-age fix deliberately didn't
  paper over this by capping age below its own required floor, since that
  would produce a self-contradictory star (e.g. a red giant younger than
  its own progenitor's main-sequence lifespan); this is the deeper fix it
  called for. `Star.generate_star`'s evolved-star mass sampling now uses a
  new `_sample_evolved_star_mass_sol` helper (`stellarObjects/starData.py`)
  that rejects and resamples (not clamps, which would just pile an
  artificial spike at the cutoff) any candidate mass whose implied
  main-sequence lifespan would exceed `UNIVERSE_AGE_GY`, capped at
  `program_constants.EVOLVED_STAR_MASS_MAX_RESAMPLE_ATTEMPTS` (100)
  attempts before raising `ValueError` -- in practice a no-op resample for
  every class but III (Giant), whose 0.8-8 Msun range straddles the
  cutoff. Yerkes class VI (subdwarf) is deliberately excluded, since its
  entire allowed mass range (0.1-0.8 Msun) sits below the cutoff and would
  reject every draw; that's tracked as a separate, still-open modeling
  question in `TODO.md`.

## [5.3.0] - 2026-09-06

### Fixed
- **Atmospheric pressure was ~5 orders of magnitude too low for every
  planet class except M.** `planetPhysics.calculate_atmospheric_conditions`
  summed "shell" volumes derived from `planet.radius`/`scale_height`
  (stored in km) directly against `planet.atm_density` (kg/m³) without
  converting km³→m³, undercounting `atmospheric_mass` by ~10⁹×; an
  unexplained `* 7500` fudge factor only clawed back about 4 of those ~9
  orders of magnitude. This rendered as "0.0 kPa" for every atmosphere-
  bearing class but M (Class M never showed it, because a hardcoded clamp
  force-overrode its pressure into an Earth-like range regardless of what
  was computed). Replaced the whole shell-integration loop with the
  closed-form barometric formula for an isothermal, hydrostatic atmosphere
  (`P_surface = ρ_surface · g · H`), which needs no arbitrary zone count or
  fudge factor and lands within the right order of magnitude for both
  Earth-like and Mars-like test cases. `tests/test_planets.py`/
  `tests/test_full_matrix.py`'s pressure assertions were tightened from a
  no-op `>= 0` check into real physical sanity bounds (1 Pa – 10 MPa) so a
  regression like this can't silently pass again.
- **Stars could be reported as hundreds of billions of years old** (e.g.
  "918.77 Billion Years old"), which is impossible since the universe
  itself is only ~13.8 billion years old — even though the star's age never
  actually exceeded its own (very long, and realistically so: real M dwarfs
  are predicted to live trillions of years) lifespan. Added a
  `UNIVERSE_AGE_GY = 13.8` ceiling, applied in both the main-sequence and
  evolved-star branches of `Star._calculate_initial_star_age_and_lifespan`
  and in `adjust_age_for_planets`'s age-adjustment logic. (White dwarfs
  already had their own bounded 0.1–12.0 Gy cooling-age range and needed no
  change.) A sub-solar-mass evolved-star progenitor whose own
  main-sequence lifespan already exceeds the universe's age is a separate,
  deeper mass-sampling issue, noted in `TODO.md` rather than papered over
  here.
- **Rare `IndexError: string index out of range` crash in name
  generation.** `generate_phoneme_salad_name` (used for star, planet, and
  moon names) could crash when a base name containing a literal apostrophe
  (e.g. `PLANET_NAMES`'s "Hi'iaka") landed adjacent to a `UNIVERSAL_PHONEMES`
  chunk also ending in an apostrophe (e.g. "ch'"), producing a literal `''`
  in the assembled name; the final capitalization step split on `'` and
  indexed `part[0]` on every piece, crashing on the empty piece between the
  two apostrophes. Not specific to binary systems — any name generation
  call could hit it given enough attempts, which is why it only surfaced
  rarely across a large batch of generated systems.

### Changed
- Class M's (and Class P's) hardcoded gravity/pressure/temperature
  "forcing" clamps in `planetPhysics.py` are commented out, not deleted —
  the fixed atmospheric-pressure formula already lands close to realistic
  ranges on its own, so the band-aid that had been masking the pressure bug
  for Class M is no longer needed. May be restored later if Class M/P need
  tighter narrative guarantees again.
- `star_systems` gained a `location` column (schema bumped to v3): for any
  system placed in a sector, its sector's name plus distance (in
  light-years) to its 3 nearest neighboring systems, e.g. `"Voranthis
  Kelmoor -- nearest: Aldenar (4.2 ly), Brekthos (7.8 ly), Corvane (9.1
  ly)"`. Computed once at write time (mirroring how `quadrant` is already
  derived and persisted) via the existing `SpaceSector.nearest_neighbors`/
  `distance_between` helpers; `NULL` for systems never placed in a sector.
  Existing v1/v2 databases migrate automatically (backed up first, as
  usual) and have their `location` backfilled from already-stored position
  data.
- The web interface's Search link is now reachable from every page (it
  previously required detouring back through the database-browse page) —
  `src/html/lib/page.py`'s shared page header now includes it whenever a
  database is selected.
- The database-picker landing page (`src/html/index.py`) now redirects straight
  to browsing the one database present, if the configured database
  directory contains exactly one `.db` file, instead of always showing the
  picker.

### Added
- `webconfig.json` (repo root, gitignored — `webconfig.json.example` is the
  committed template): site-level configuration, currently `site_name` and
  `base_url`, plus unused placeholder fields (`db_username`, `db_password`,
  `db_name`) reserved for a possible future non-SQLite backend. Kept
  outside `src/html/`'s served document root, the same way `db/` already is.
  See [`docs/webconfig.md`](docs/webconfig.md) for full documentation.

## [5.2.5] - 2026-09-06

### Changed
- **Database schema bumped to v2**: moons split out of the shared
  `planets` table into their own `moons` table (`planet_id` FK to the
  planet they orbit, plus their own `moon_evolutionary_paragraphs`/
  `moon_reflection_spectrum` child tables), instead of self-referencing
  via `planets.parent_planet_id`/`is_moon`. This is what `src/html/search.py`'s
  attribute tags needed to actually tell a planet from a moon: previously
  a "Class D Planet" tag queried `planets` with no way to exclude
  `is_moon=1` rows of the same class, so it silently listed moons too.
  `stellarObjects/_db.py` gained a dedicated `insert_moon` (mirroring
  `insert_planet`, but writing to the new table); `insert_planet` no
  longer recurses into itself for moons.
- `src/html/search.py`'s tag facets and name search now follow the same
  split: "Planet Class"/"Planet Body Type"/"Planet Supported Life
  Chemistry" only ever match top-level planets, with an identically-
  shaped "Moon Class"/"Moon Body Type"/"Moon Supported Life Chemistry"
  set of tags (and a separate autocompleting "Moon name" field) for
  moons -- each still only rendering the tag buttons for values actually
  present in the chosen database. `src/html/system.py`'s planet/moon table
  now reads moons from the new table instead of a recursive
  self-join, and no longer needs to recurse (moons never generate their
  own moons -- confirmed by the existing
  `tests/test_moons.py::test_moons_cannot_themselves_have_moons`).

### Added
- `migrateDb.py` (repo root): converts every `*.db` file in a directory
  from schema v1 to v2, backing up each original first
  (`<name>.db.v1-backup-<timestamp>.db`) -- a no-op for a database
  that's already current. Backed by a new `stellarObjects._db.
  migrate_database`/`_migrate_v1_to_v2`, which builds the converted
  database at a temporary path and only atomically swaps it into place
  at the very end, so a crash or error partway through leaves the
  original file (and the backup already made) untouched. `install.sh`
  now runs this automatically (a new step, right after installing the
  Python package) over every database in `db/`, so `sudo ./update.sh`
  keeps an existing deployment's database working across the schema
  change with no manual step.
- `tests/test_db_migration.py` and `tests/test_db_persistence.py`: the
  first automated tests of `stellarObjects/_db.py`'s save/migrate paths
  (previously "manually smoke-tested via `sectorGen.py`" per `TODO.md`).
  Between the two, they cover the v1->v2 migration (including the
  moon-owned child-row split, which the real generated data used for
  manual smoke-testing happened not to exercise) and a genuine
  generate-then-save-then-query round trip through `insert_star_system` --
  the latter caught a real bug during development (`insert_moon`'s
  `INSERT` had one more column than value placeholder).

## [5.2.4] - 2026-09-06

### Added
- `src/html/search.py`: a faceted/name search page for the web interface.
  Two complementary ways to find an object in the chosen database:
  - Click-to-filter tag buttons for object type (star/planet/moon/
    asteroid belt), star spectral class and Yerkes luminosity class,
    planet class and body type, and supported life chemistry -- each
    group renders a button only for values actually present in that
    specific database (e.g. no Yerkes-luminosity buttons beyond "Main
    Sequence" if nothing else was generated), per `SELECT DISTINCT ...
    GROUP BY` queries against `stars`/`planets`/`asteroid_belts`, not a
    fixed enumeration. Multiple tags within one group OR together (any
    matching value); tags across groups AND together, except that
    selecting an explicit Object Type tag acts as a master filter (e.g.
    selecting only "Stars" hides the Planets panel even if a
    planet-class tag also happens to be selected). Clicking a tag
    toggles it via a plain link that rewrites the query string, so this
    works with JavaScript entirely disabled, same as every other page in
    `src/html/`.
  - A name search with one field each for sector, star system, star, and
    planet/moon names, each with its own HTML5 `<datalist>` for
    autocomplete (native browser suggestions, no JavaScript, populated
    from that entity's own distinct names). Asteroid belts have no name
    of their own (per `schema.sql`'s own note on this), so they're only
    reachable via the "Asteroid Belt" object-type tag.

  Every filter is a parameterized SQL query (`IN (...)`/`LIKE ... ESCAPE
  '\'`), and each result panel is capped (300 rows, with a "showing the
  first N" note) to keep a broad tag click from dumping an entire large
  database into one page. Linked from `browse.py`'s breadcrumb.

## [5.2.3] - 2026-09-06

### Changed
- `install.sh` no longer installs the Python package via the deprecated
  `python3 setup.py install`; it now runs a proper, build-isolated `pip
  install --upgrade --force-reinstall "$SCRIPT_DIR"` instead.
  `setup.py install` broke on a deployed Ubuntu 20.04/Python 3.8 host
  *twice* in a row, both times because setuptools' own vendoring shim
  (`extern`) prefers a real, already-installed copy of a dependency it
  vendors over its own newer bundled copy whenever one is importable —
  first `importlib_metadata` (`AttributeError: ... no attribute
  'EntryPoints'`), then, after that was patched around, `packaging`
  (`TypeError: canonicalize_version() got an unexpected keyword argument
  'strip_trailing_zero'`) — because that distribution's old apt-provided
  copies of each in turn shadowed the working vendored one. Chasing each
  shadowed dependency individually only fixes the one that broke that
  day, not the next one; `pip install`'s build isolation builds the
  package in a throwaway environment that can't see the system's
  site-packages at all, closing the whole bug class instead of patching
  it dependency-by-dependency, and means `install.sh` no longer needs to
  upgrade the system's global `setuptools`/`importlib_metadata` at all
  for this step. `--force-reinstall` is deliberate: a plain `pip install
  .` skips reinstalling when pip thinks the same version is already
  installed — true on every `update.sh` run between version bumps in
  `stellarObjects/_version.py` — which would otherwise silently leave the
  previous run's install in place instead of the source `update.sh` just
  pulled.
- `install.sh` and `update.sh` now re-`chmod +x` every `*.sh` file in the
  repo (not just `src/html/*.py` and one hardcoded `apache/set-permissions.sh`
  path), including themselves. A `core.fileMode=false` git config on the
  authoring machine drops the executable bit on any file type on
  checkout, not just `src/html/`'s `.py` scripts, and `update.sh` invokes
  `install.sh` directly (`"$SCRIPT_DIR/install.sh"`, not `bash
  install.sh`) — if a pull had dropped *its* executable bit, the shell
  would refuse to exec it before `install.sh`'s own permission fix ever
  got a chance to run. `update.sh` now fixes this itself, immediately
  before invoking `install.sh`, so a dropped bit on `install.sh` (or on
  `update.sh` for its own next run) no longer requires a manual `chmod
  +x` to recover from. (One unavoidable exception: the very first
  `update.sh` run to pull this fix is still executing the old script
  text it already had open when `git pull` replaced the file on disk out
  from under it, so that one migration may still need a manual `chmod +x
  install.sh` — every run after that starts from the fixed script.)

## [5.2.2] - 2026-09-06

### Added
- `src/html/system.py` now renders a system's description as actual HTML by
  default (`?view=rendered`), via a new small, purpose-built Markdown-to-
  HTML converter (`src/html/lib/mdconvert.py`) targeting exactly the narrow
  Markdown subset `StarSystem.__str__` generates (headers, pipe tables,
  paragraphs, `<sup>` exponents) — not a general-purpose parser, and no
  new dependency. The original raw wikitext/Markdown source is still one
  click away (`?view=source&format=...`) for copy-pasting into a wiki.
  Escapes every block in full before emitting markup, then narrowly
  re-enables only the one legitimate raw-HTML pattern generated content
  contains, so a mischievous `--name`/`--star-type` value can't inject
  live HTML into a rendered page (covered by new tests in
  `tests/test_mdconvert.py`).
- `src/html/static/style.css` rewritten as a small design system: CSS custom
  properties, automatic light/dark via `prefers-color-scheme`, card-style
  panels, badges, breadcrumbs, a proper type scale, hover states,
  focus-visible outlines, and a responsive breakpoint. Applied
  consistently across every page (`index.py`, `browse.py`, `sector.py`,
  `system.py`, and the shared error page in `lib/page.py`), not just the
  system page.
- `update.sh` (repo root, alongside `install.sh`): pulls the latest
  changes from git and re-runs `install.sh` so permissions stay correct
  afterward. Refuses to run over uncommitted local changes, and pulls
  with `--ff-only` (fails loudly rather than creating a surprise merge
  commit if history has diverged) instead of a plain `git pull`.

### Fixed
- `apache/set-permissions.sh` and `install.sh` both failed to make every
  `*.py` file under `src/html/` executable by `www-data` — `install.sh`'s own
  `chmod +x` step used `find -maxdepth 1`, silently skipping
  `src/html/lib/*.py`, and `set-permissions.sh` only `chgrp`'d (group
  ownership) rather than `chown`'d (user *and* group) the deployed
  directories. Both fixed: the `-maxdepth 1` restriction is gone, and
  `set-permissions.sh` now `chown -R`s to the detected Apache user:group
  and reports how many `.py` files it made executable, so a wrong path is
  obvious rather than silently matching nothing.

## [5.2.1] - 2026-09-06

### Added
- `src/html/`: a small, dependency-free web interface (plain Python CGI
  scripts, standard library only) for browsing the SQLite databases from
  `db/README.md` — pick a database, drill into its sectors and star
  systems, and view (or copy, via a `<textarea>`) the rendered
  wikitext/Markdown page saved for each one. See `src/html/README.md`.
- `apache/`: deployment tooling for the web interface —
  `planetgen.conf.example` (an example Apache2 virtual host) and
  `set-permissions.sh` (detects the user/group Apache2 actually runs as
  and sets ownership/permissions on the deployed `src/html/`/`db/`
  directories accordingly).
- `install.sh`: a one-shot Linux installer tying the above together —
  runs `setup.py install`, pre-fetches the NLTK `words` corpus into a
  shared, world-readable location, makes the CGI scripts executable,
  enables Apache's `cgid` module, and runs `apache/set-permissions.sh`;
  prints the one remaining manual step (copying/enabling the example
  vhost) rather than touching Apache's site configuration itself.
- `stellarObjects/names.py` gained `UNIVERSAL_PHONEMES`: a pool of ~100
  short, ASCII-7-bit-printable phoneme chunks romanized from roughly 18
  language families (Romance, Germanic, Slavic, Arabic/Hebrew/Persian,
  Turkish, South Asian, Mandarin, Japanese, Korean, Vietnamese,
  Austronesian, Polynesian, Bantu, Mesoamerican, Andean, Celtic,
  Finno-Ugric, Caucasian). `generate_phoneme_salad_name` now splices one
  of these into every generated name — star, planet, moon, and sector
  alike — with `UNIVERSAL_PHONEME_CHANCE` (40%) odds, widening the
  cultural range generated names are drawn from beyond each type's own
  base name list.

### Changed
- Removed `SECTOR_DESIGNATORS` from `stellarObjects/names.py` (dead code
  — defined but never referenced anywhere).
- `generate_phoneme_salad_name` gained an `allow_split` parameter
  (default `True`, unchanged for stars/planets/moons); see the sector
  name bug fix below.

### Fixed
- Sector names could come out as 3-4 words instead of the intended 2 (in
  a 500-name stress test, this happened 95% of the time).
  `generate_phoneme_salad_name` can split a long result into two
  space-separated words on its own (`split_long_word`), and
  `sectorGen.generate_sector_name` already joins two independent calls
  into one name — an internal split on either half silently produced
  3-4 words in the final result. `generate_sector_name` now passes the
  new `allow_split=False` for both halves.
- `stellarObjects/names.py` unconditionally called
  `nltk.download('words', quiet=True)` at import time; `nltk`'s
  `download()` always targets the *current user's* default download
  directory and attempts to create it, regardless of whether the corpus
  already exists elsewhere on `nltk.data.path` — this broke the web
  interface entirely under Apache's locked-down `www-data` account
  (`PermissionError: ... '/var/www/nltk_data'`), since that account has
  no writable home directory. It had only ever "worked" for the CLI
  tools because those ran as `root`. Now checks
  `nltk.data.find('corpora/words')` first and only downloads on
  `LookupError`; `install.sh` pre-fetches the corpus into a shared,
  world-readable path so that check always succeeds once installed.
- `apache/set-permissions.sh` hard-failed when it couldn't detect
  Apache's user/group (e.g. because apache2 isn't started/enabled yet,
  which is exactly the case during a fresh `install.sh` run) — now warns
  and defaults to the standard Debian/Ubuntu `www-data:www-data` instead.
- CGI scripts could be deployed non-executable regardless of a local
  `chmod +x`, because a `core.fileMode=false` git config on the
  authoring machine silently drops the executable bit before it reaches
  a commit. `install.sh` now `chmod +x`s the CGI scripts and
  `apache/set-permissions.sh` directly on every install, independent of
  whatever mode git happened to store.

## [5.2.0] - 2026-09-06

### Added
- Full SQLite database persistence: `stellarObjects/schema.sql` defines the
  schema, and the new `stellarObjects/_db.py` writes an already-generated
  `SpaceSector` (every system, star, planet, moon, and asteroid belt,
  plus a rendered copy of the wiki page in both wikitext and Markdown)
  into it in a single transaction. Documented column-by-column in the new
  `db/README.md`. `sectorGen.py` now calls this automatically on every
  run, saving to `db/planetgen.db` by default (overridable via the new
  `--db-path` option); the database file itself is gitignored.
- Sector name generation: sectors now get a random two-word name (e.g.
  "Voranthis Kelmoor") drawn from a new sector-flavored name list
  (`SECTOR_NAMES`/`SECTOR_PREFIXES`/`SECTOR_SUFFIXES` in
  `stellarObjects/names.py`, pulling from real galactic structures and
  well-known science-fiction sector names) instead of reusing star names
  with "Sector" appended.
- `Star`, `Planet`, `BinaryStarProxy`, and `AsteroidBelt` each gained a
  `get_table_properties()`/`get_composition_summary()` method that returns
  the same already-formatted values their `to_paragraph_list()` renders
  into text, so the new database layer can store exactly what was
  published without duplicating any formatting logic.

### Changed
- **Breaking:** `sectorGen.py`'s `-n` short flag now means `--name`
  (hard-sets the sector's own name) instead of `--num-systems`, which no
  longer has a short flag; the old `--sector-name` option was renamed to
  `--name`.
- Distances stored in the database use a two-tier unit convention
  (milliparsecs for sector-scale position/geometry, kilometers everywhere
  else); `stellarObjects/utils.py` gained `ly_to_milliparsecs`/
  `milliparsecs_to_ly` and `physical_constants.py` gained
  `AU_PER_PARSEC`/`AU_PER_MILLIPARSEC` to support the conversion, used only
  by the persistence layer.

### Fixed
- Atmospheric scale height (`planetPhysics.calculate_atmospheric_conditions`)
  was computed in meters but combined directly with `planet.radius` (in
  km) without converting units first, throwing off atmosphere thickness
  and volume for every planet with an atmosphere.
- `utils.calculate_object_mass` and `Planet.__init__` each independently
  (and inconsistently) converted radius-in-km to a volume, one of them
  mixing up its km/m factor; `Planet` now takes its volume solely from
  `calculate_object_mass`'s corrected calculation instead of recomputing
  it a second time.

## [5.1.0] - 2026-09-06

### Added
- `--version` option on both `systemGen.py` and `sectorGen.py`, printing the
  program's version, this repository's URL
  (https://github.com/dwhagar/planetGen), and a license summary, then
  exiting immediately.
- `stellarObjects/_version.py`: a single, dependency-free source of truth
  for the project's version number, shared by both CLI scripts and
  `setup.py`.

### Changed
- `setup.py` now reads its `version` from `stellarObjects/_version.py`
  instead of a hardcoded, never-updated placeholder.
- README updated with current version information and a link to this
  changelog.

## [5.0.1] - 2026-09-06

### Fixed
- Flavor text (both system-level and planet/moon-level) was being rolled —
  and shared `system_config` counters mutated — every time a system or planet
  was *rendered* (`__str__`/`to_paragraph_list()`) instead of once at
  generation time, so rendering the same object twice would silently re-roll
  and double-count it. It is now decided exactly once per object during
  generation (`StarSystem.__init__`), and rendering is a pure, idempotent
  read of the already-decided text.

### Changed
- Internal property cleanup across `Planet`, `Star`, and `AsteroidBelt`:
  removed five redundant star-property snapshots from `Planet` that
  duplicated data already reachable via its `star` reference (one of
  which — `star_radius` — had a latent double-unit-conversion bug whenever a
  moon was generated); renamed ambiguous or colliding attributes
  (`evolution` → `evolutionary_speed`, `hab` → `habitable_zone`, and
  `type` → `body_type` on `Planet`/`AsteroidBelt` so it no longer shares a
  name with the unrelated `Star.type`); filled in missing/stale docstring
  documentation.

### Added
- Regression tests proving rendering is idempotent (calling
  `to_paragraph_list()`/`__str__()` twice produces identical output and
  performs no further mutation).

## [5.0.0] - 2026-09-05 (evening)

### Added
- `spaceSector.py`: a new `SpaceSector`/`SectorSystemEntry` layer that places
  generated `StarSystem`s at `(x, y, z)` positions within a cubic sector, with
  save/load to JSON (`SystemConfig` gains `to_dict`/`from_dict` to support
  this).
- Poisson-disk-based sector growth for realistic inter-system spacing, plus
  named-location ("quadrant") formatting for a system's position within a
  sector.
- `TODO.md`: a long-term roadmap toward full database storage and a web
  interface.

### Changed
- Reworked sector growth to sample candidate positions and then fine-tune
  each one against every existing neighbor, rather than a single-pass
  placement.
- Sector density and the minimum allowed separation between two systems are
  now grounded in real astronomical data — density from real local stellar
  surveys, and minimum separation from each system's own gravitationally
  derived Hill-sphere radius — rather than arbitrary constants.

## [4.1.0] - 2026-09-05 (afternoon)

### Added
- JSON system-file support: a system's exact contents (star type, forced
  features, per-orbit slots) can now be specified via a JSON recipe file,
  with several example systems included under `examples/`.
- The project's first automated test suite (pytest), covering example
  systems, moon generation, and planet generation.
- `sectorGen.py` CLI entry point, ahead of the sector-generation logic it
  would call the next day.

### Changed
- Major standardization pass across CLI options and internal call
  signatures.
- Renamed the main script from `planetGen.py` to `systemGen.py`, reflecting
  that it now generates full star systems rather than a single planet.

## [4.0.0] - 2026-09-03

### Changed
- Major internal refactor: generation logic was split out of the `Planet`
  class into dedicated `planetPhysics.py` (physical/orbital generation) and
  `planetLife.py` (life chemistry and evolution) modules, so `Planet` itself
  holds state and presentation while generation logic lives alongside it as
  free functions.
- The monolithic `constants.py` was split into `physical_constants.py`
  (real-world physical constants) and `program_constants.py`
  (generation/tuning constants), for readability and maintainability.
- Import structure cleaned up throughout the package; install script
  (`setup.py`) updated.

## [3.0.0] - 2026-07-07

### Added
- Binary star systems: a new `doubleStar.py` module (`BinaryStarProxy`)
  represents a double star system as a single effective star (combined
  mass/luminosity/habitable zone) for the purposes of planet placement.

### Fixed
- Primary/secondary star identification and overall system age calculation
  for binary systems.
- An age/lifespan bug affecting planets.

### Changed
- Binary star output formatting refined; the large-star-forcing option
  updated to work correctly alongside binary generation; README updated.

## [2.2.3] - 2026-06-27

### Added
- More flavor text variants across different planet classes.

### Changed
- Refined repeat-prevention and selection logic for flavor text so the same
  text is less likely to recur in quick succession.

## [2.2.2] - 2026-06-26

### Fixed
- Flavor text selection options.

### Changed
- Switched random-number generation to Python's `secrets` module for
  higher-quality entropy, avoiding similar successive sequences.

## [2.2.1] - 2026-06-25

### Changed
- Wording and grammar adjustments throughout the generated text output.

## [2.2.0] - 2026-06-24

### Added
- Options to specify a system's age and to force (or forbid) intelligent
  life.
- Flavor text: a random chance of extra descriptive "sensor" text being
  appended to a system or a planet's description.

### Changed
- Centralized magic numbers into `constants.py` as named constants, making
  the physics/generation formulas easier to read.
- Unified the set of properties referenced across all object types onto a
  shared convention.
- README updated for the new CLI options.

## [2.1.1] - 2026-06-23

### Changed
- Asteroid belt data and logic split out into its own module
  (`asteroidData.py`).
- Expanded asteroid belt composition descriptions and grammar.

### Added
- Option to specify a custom name for the star system.

## [2.1.0] - 2026-06-21

### Added
- Stellar age is now generated and reported for every system.
- Evolutionary timeline narratives: systems now estimate the plausibility of
  life at every stage, from simple single cells through technological
  civilizations.

### Fixed
- Assorted output-formatting issues.

## [2.0.0] - 2026-06-20

### Added
- Life chemistry system: planets and star types now carry information about
  the chemical processes that could plausibly give rise to life, as a
  foundation for the evolutionary modeling that follows.

### Changed
- Command-line options reworked so multiple flags combine correctly together.

## [1.5.0] - 2026-06-12

### Added
- Option to generate a star system with no planets.
- Option to specify a star's exact spectral type from the command line.

## [1.4.0] - 2026-06-06

### Added
- Markdown export support, and the ability to write generated output to a
  file (in addition to the console).
- Heliosphere radius and stellar gravitational-influence ("system perimeter")
  calculations, describing the outer boundaries of a system.

### Fixed
- Name-generation issues.
- Scientific-notation formatting issues.

## [1.3.1] - 2026-06-05

### Changed
- Further cleanup of generated text output for clarity.
- Refinements to procedural naming.

## [1.3.0] - 2026-06-04

### Added
- Procedural name generation for stars, planets, and moons.
- Project README.

### Changed
- General code cleanup pass.

*(~23-month gap in development between 2024-07-15 and 2026-06-04.)*

## [1.2.2] - 2024-07-15

### Fixed
- Star age/temperature calculations corrected to consistently use Kelvin
  throughout.

## [1.2.1] - 2024-07-14

### Fixed
- Class M worlds are now reliably clamped to Earth-like gravity and
  atmospheric pressure.
- Class P (and other) planet classes now get appropriate temperature
  treatment for their class.

### Changed
- Asteroid belt sizing now makes better use of the habitable zone.

## [1.2.0] - 2024-07-13

### Added
- Options to force a habitable world, a large star, and/or an asteroid belt
  (forcing both a habitable world and an asteroid belt together
  automatically forces a large star).
- Options to control the overall size (object count) of a generated system.
- Expanded command-line help text.

### Fixed
- A long-standing infinite-loop bug in star generation.
- Several edge cases in generation logic.

### Changed
- Capped the maximum number of system objects at 500 to prevent runaway
  generations.
- Refined the system description text and radius output formatting.

## [1.1.1] - 2024-07-11

### Fixed
- Completed and debugged the moon system: moon orbital distances are now
  tracked and calculated correctly, and remaining orbital overlap between
  planets and asteroid belts was removed.

### Added
- Weighted probabilities for moon class selection.

### Changed
- Further wording, grammar, and formatting refinements.

## [1.1.0] - 2024-06-29

### Added
- First implementation of the moon-generation function (not yet wired into
  planet creation or tested at commit time).

## [1.0.2] - 2024-06-28

### Added
- Mass validation for a specified planet, laying the groundwork for moon
  generation.

## [1.0.1] - 2024-06-26

### Fixed
- Orbital placement now accounts for planetary position and Hill radius, so
  systems no longer generate overlapping or too-closely-spaced orbits.
- Asteroid belt inner/outer bounds no longer overlap neighboring planets.

### Changed
- Further readability tweaks to the generated text output.

## [1.0.0] - 2024-06-25

### Added
- First fully working end-to-end planet generator, validated against real
  Earth reference values ("I think it finally works!").
- Special-cased handling for Class N worlds (boosted atmospheric density and
  molar mass).

### Changed
- Multiple formatting and descriptive-text passes on the generated output.

## [0.2.0] - 2024-06-23

### Added
- Core stellar physics: working star generation and planet-count estimation.
- Initial atmospheric modeling: mass, gravity, and atmosphere calculations.

### Fixed
- Numerous early bugs in surface pressure and temperature calculations.

## [0.1.1] - 2024-06-22

### Fixed
- Continued debugging to get the initial prototype running end-to-end (moved
  development from VS Code to PyCharm along the way).

## [0.1.0] - 2024-06-21

### Added
- Initial project scaffolding: the `stellarObjects` package with the first
  `Planet`, `Star`, and `StarSystem` classes, and a `main.py` entry point.
- First rough pass at generating a single planet's basic properties.