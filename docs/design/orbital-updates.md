# Orbital updates, positions and spin

How planetGen will move things: one position object for every body, a
positional update that only touches what has visibly moved, gravity from
nearby masses, collisions, spin, and the limits where the math stops
being trusted. This condenses Boss's documents in this folder:

- "Orbital Update Full Algorithm.md" (the nightly batch: data structures,
  passes, helpers)
- "Orbital Position and Vector Update Algorithms.md" (thresholds, grid
  resolution, per-sector point-mass tables)
- "Orbital Update Collision Detection.md" (swept-sphere collisions,
  mergers)
- "Observational Kinetics for Rotational Vectors.md" (spin rates, axial
  tilt)
- "Orbital Position and Vector Mathematical Foundations.md" (the numerical
  methods underneath: root finding, optimization, linear solvers, eigen
  and singular value decompositions, and where each breaks down)
- "Computational Astrodynamics.md" (Boss's research of 2026-10-07: the
  galaxy's own gravity, the three frames, the orbit solvers, nonsingular
  elements, the time step and epoch, and collisions)
- `spacial-position.py` at the repo root (Boss's prototype position class)

## 1. What Boss asked for

- (2026-10-03) Only update and count objects that moved noticeably: 0.01
  milliparsec for star-like objects on the galactic scale, 0.01 AU for
  planets and binary stars in a system, 100,000 km for moons; store when
  each object is next due so a run finds them fast.
- (2026-10-03) Report how many objects changed sector or entered or left
  a nebula.
- (2026-10-03) Paths bent by the nearest 10 bodies at least as massive,
  only at updates; warn when an object is inside another's Hill radius,
  and email the admin when email is set up.
- (2026-10-03) Rogue planet collisions: a terrestrial one makes an
  asteroid field; two gas giants merge, and if the result can fuse it
  becomes a star with a planetary nebula and no planets, reported to the
  admin.
- (2026-10-07) One point-in-space object that keeps every coordinate
  system in step, used by every object, with its point-mass data; a spin
  vector and a realistic axial tilt for every rotating body; editable
  trajectories; trajectories shown in their own frame.
- (2026-10-07) "we need to ensure that reasonable limitations for when
  the math breaks down at the edge cases."

## 2. The position object

`planetgen/physics/position.py`'s `SpatialPosition3D` (Boss's prototype `spacial-position.py`, now removed) keeps one body's position in
three frames (galactic, sector, system) and three forms (Cartesian,
cylindrical, spherical). Setting any coordinate in any frame and form
updates all of them, with the velocity, and it works out how long until
the body moves enough to be observed. It moves into the physics package
with tests, and every positioned object holds one, with its mass and
gravitational parameter (mu = G M) beside it.

Built (GEN.74, part 1): all lengths are metres, speeds m/s, masses kg
(the storage layer converts at its edge). The position last set is kept
as given in its own frame and the others are derived from it, so a moon
set by its system offset keeps its metres exactly. The prototype's
placeholder minimum-observable constants are replaced by the thresholds
below (`THRESHOLDS_M`: galactic, system, planetary), with the next-due time
capped at a billion years. Lengths are in a unit the object chooses (`length_unit_m`: light-years for sector entries, AU inside a system), so no round trip through metres touches stored values. Values the maths cannot hold (non-finite
numbers, a negative radius, a polar angle outside [0, pi], a speed at or
above the speed of light, a negative mass) raise `ValueError`.

## 3. Thresholds and the next-due time

| Scale | Objects | Update when moved at least |
|---|---|---|
| Galactic | stars, black holes, neutron stars, brown dwarfs, rogue bodies (anything star-like that isn't in a system) | 0.01 mpc (about 2 AU) |
| System | planets, companion stars | 0.01 AU |
| Planetary | moons and other satellites | 100,000 km |

When an object is moved for any reason, its next-due time is set from its
speed and its threshold (time = threshold / speed, capped). An indexed
`next_update_due` column lets the run select only what is due. Objects
that didn't pass their threshold are not moved or counted.

## 4. The update run

From "Orbital Update Full Algorithm.md", adapted to planetGen:

1. **Galactic pass**: for each due root object (stars, remnants, rogue
   bodies), gather point masses from its sector and the neighbouring
   sectors, keep the top 10 influencers of equal or greater mass (plus
   the galactic anchors), and advance it with Velocity Verlet. Check its
   swept path against those influencers for collisions; record close
   encounters.
2. **Encounters**: resolve each close pair with a finer step.
3. **System pass**: move each system's children by the star's
   displacement, then advance their orbits with a Kepler propagator; check
   Hill-sphere crossings (capture if bound, deflection if not, disruption
   inside the Roche limit, ejection).
4. **Bookkeeping**: rebuild the point-mass table of every sector that
   changed, refresh Hill-sphere occupants for objects that moved far
   enough, and write the summary: objects moved, sector changes, nebula
   entries and exits (by the nebula shape test), Hill-radius warnings,
   collisions.

**Point-mass tables**: each sector stores a small table of its masses
(id, position, mass, velocity; about 128 bytes a row), rebuilt when a
sector is generated, when something moves in or out, and when something
is created or destroyed. The run loads tables for the sectors it needs and
builds a temporary in-memory tree; vector fields are never stored.

**Sector geometry**: the documents assume 11.5 ly cylindrical sectors
(and the anomaly documents 20 ly cubes). planetGen's sectors are 4 pc
(about 13 ly) cells of the ring, layer and slot grid
(galaxy-coordinate-system.md), so neighbour lookup uses that grid's own
address arithmetic.

## 5. Collisions

Continuous collision detection on swept spheres: two bodies collide in a
step when the quadratic for their closest approach has a root inside the
step at a distance below the sum of their radii (with gravitational
focusing). A collision is an inelastic merger that keeps momentum. Rogue
planet rules on top of that are Boss's (section 1); a merged gas giant
above the deuterium or hydrogen burning limit becomes a brown dwarf or a
star, reported to the admin.

## 6. Spin and axial tilt

Every rotating body stores a spin axis and rate, drawn by the kinetics
document's cascade:

| Body | Rule |
|---|---|
| Tidally locked (system age past the locking time) | rotation period equals orbital period; tilt 0, or a 3:2 resonance when e > 0.1 |
| Cool main-sequence stars (up to 1.3 M_sun) | gyrochronology: period from age and colour, with 10% log-normal scatter |
| Hot stars | equatorial speed log-normal around 180 km/s, capped at 85% of breakup |
| Small bodies over 200 m | period log-normal, never below the 2.2-hour spin barrier |
| Black holes | spin a* from Beta(1.4, 3.6), or about 0.69 for merger remnants |
| Stellar tilt | Rayleigh, sigma 15 degrees |
| Planet tilt | impact-modified distribution |
| Asteroid tilt | near 10 or 170 degrees (YORP) |

The axis is the orbit normal tilted by the obliquity at a random
precession angle.

## 7. Where the math breaks down

Each method has a range where it is trusted; outside it the code applies a
guard, logs it, and tests cover it:

- **Near-parabolic orbits** (e close to 1): Kepler's equation converges
  slowly; switch to the universal-variable (or Barker) form.
- **Steps longer than an orbit**: propagate analytically, never by
  integration.
- **Very close passes**: an encounter inside the bodies' radii is a
  collision, not a slingshot; inside the Roche limit, a disruption.
- **Inside a Hill sphere**: hand the body to the host's system pass
  instead of the galactic one.
- **The galactic centre**: the central mass's potential is softened, so
  nothing gets an infinite kick.
- **Precision at galactic distances**: positions are stored relative to
  the sector (and system) frame, not as absolute metres, so small moves
  aren't lost in rounding; the position object handles the frames.
- **Unbound results**: a body flung beyond the galaxy's bounds is kept
  and flagged, not deleted.

## 8. Light-travel positions

What an observer sees is where an object was, not where it is: its
apparent position is its position at (now - distance / c), solved by
iteration for moving bodies. This is the groundwork for the view from a
planet.

## 9. Numerical methods and their limits

The foundations document covers the solvers the update leans on, how each
one fails, and what to fall back to. In planetGen terms:

| Problem | Method | Where it breaks | Fallback |
|---|---|---|---|
| Kepler's equation (anomaly from time), Hill and Roche radii, the light-travel time in section 8 | Newton's method, quadratic near the root | the derivative 1 - e cos E goes to zero as e nears 1 near periapsis; a poor first guess wanders | Brent's method on a bracket that always holds the root ([M - e, M + e] for ellipses), which halves the bracket every step and cannot diverge; past e of about 0.99, the universal-variable form (section 7) |
| Stopping any iteration | tolerance 2 eps \|x\| + an absolute floor | a bracket narrower than one floating-point step never shrinks, so the loop never ends | stop at that width, cap the iteration count, and log the case |
| Editable trajectories and course fitting (phase 2 and 3) | BFGS (scipy.optimize) with a Wolfe line search | curvature turns non-positive and the step stops going downhill | Powell's damped update, then reset to steepest descent |
| Small linear systems (frame changes, encounter fits) | LU with partial pivoting (numpy and LAPACK) | element growth, or a matrix nearly singular | QR, then a regularized least-squares solve, and the result is flagged |
| Rotation matrices for frames and spin axes | products of rotations | round-off drifts them away from orthonormal after many updates | re-orthonormalize with the SVD (or a quaternion renormalize) on every save |

Precision: positions are 64-bit floats (epsilon about 2.2e-16). Across
the 30 kpc galaxy that is about 0.2 km in absolute coordinates, far below
the 0.01 mpc galactic threshold. Systems and moons still use their local
frames (section 7), because a moon's 100,000 km threshold must survive
being added to a star's galactic position. A problem whose condition
number nears 1 / epsilon (nearly equal eigenvalues, near-singular
matrices) is logged and solved at higher precision (mpmath), never
silently trusted.

Not needed: the Lanczos and dense eigenvalue methods. planetGen has no
large sparse matrix problem, and the 3 by 3 decompositions it does need
come from numpy.

## 10. Galactic gravity, frames, solvers and time

From "Computational Astrodynamics.md". Every number below was rerun
here; section 10.6 lists the places the document is wrong and the fix
the code uses.

### 10.1 The galaxy's own gravity (GEN.115)

A star feels a smooth background potential plus its nearby point masses
as perturbations on top:

`Phi_gal(R, z) = Phi_bulge(r) + Phi_disk(R, z) + Phi_halo(r)`

| Part | Model | Potential | Default mass | Default scale |
|---|---|---|---|---|
| Bulge | Hernquist | `-G M_b / (r + c_b)` | 1.0e10 M_sun | c_b = 0.5 kpc |
| Disk | Miyamoto-Nagai | `-G M_d / sqrt(R^2 + (a_d + sqrt(z^2 + b_d^2))^2)` | 6.8e10 M_sun | a_d = 3.5 kpc, b_d = 0.3 kpc |
| Halo | NFW | `-(G M_h / r) ln(1 + r / r_h)` | 5.4e11 M_sun | r_h = 16 kpc |
| Halo (alternative) | Flattened log | `0.5 v0^2 ln(R_c^2 + R^2 + z^2/q^2)` | v0 = 175 km/s | R_c = 2.5 kpc, q = 0.9 |

The default is the NFW halo; the log halo is a setting. With
G = 4.300917e-6 kpc (km/s)^2 / M_sun the circular speed in the midplane
is `sqrt(v_b^2 + v_d^2 + v_h^2)`. Rerun here:

| R (kpc) | 2 | 4 | 6 | 8 | 8.128 | 10 | 12 | 15 | 20 | 30 |
|---|---|---|---|---|---|---|---|---|---|---|
| v_circ (km/s) | 190.5 | 223.2 | 230.6 | 229.4 | 229.3 | 226.3 | 223.1 | 218.9 | 213.5 | 205.4 |

At 8.128 kpc the parts are v_b = 68.5, v_d = 163.6 and v_h = 145.3 km/s,
matching the document and the observed 229 to 232 km/s. These are the
values GEN.115's tests check at the default galaxy shape.

Point masses use Plummer softening, `a = G M r / (r^2 + eps^2)^(3/2)` with
eps = 1 pc, so the pull near the central black hole and in dense
clusters stays finite and falls to zero at the centre.

### 10.2 Frames

| Frame | Origin | Axes | To its parent |
|---|---|---|---|
| Galactic | The central black hole | +Z galactic north (along the disk's spin), +X the zero meridian | Master frame |
| Sector | See below | Physics uses the galactic axes, unrotated | Translation only |
| System | The host star's barycentre | +z along the star's spin, +x along the line where the star's equator crosses the galactic midplane (`Z_gal x z_sys`, or `X_gal` when the spin points along `Z_gal`), +y = z x x | `p_gal = p_star + R^T r_sys` |

The star's pole is given as a right ascension and declination measured
in the galactic frame (not Earth's sky), so
`z_sys = (cos dec cos ra, cos dec sin ra, sin dec)`,
`x_sys = (-sin ra, cos ra, 0)`,
`y_sys = (-sin dec cos ra, -sin dec sin ra, cos dec)`, and
`R_gal_to_sys` has those three as its rows (checked: orthonormal, right
handed). Planet inclinations and tilts are measured from the star's
equator.

### 10.3 Solvers

| Routine | What it does | Notes |
|---|---|---|
| PropagateKeplerianOrbit | Moves a two-body orbit by any time step through the universal variable chi and the Stumpff functions c0 to c3, solved by Halley's method to 1e-12, then Lagrange f and g | One routine for circles, ellipses, parabolas and hyperbolas. Use the corrected equation in 10.6. |
| ComputeTwoBodyHyperbolicDeflection | A flyby that stays unbound: e = sqrt(1 + (b v_inf^2 / mu)^2), turning angle 2 asin(1/e), the relative velocity rotated about h by Rodrigues' formula, the change split by mass | Momentum is conserved. v_inf must be the speed at infinity (10.6). |
| ResolveCloseEncounterMicroPass | Inside a mutual Hill sphere the step pauses for that pair and a 4th-order Hermite or Yoshida integrator sub-steps at `eta sqrt(r^3 / G(M1+M2))`, eta 0.01 to 0.05 | Fluid Roche limit `2.44 R1 (rho1/rho2)^(1/3)` disrupts the smaller body into ring debris; contact merges them. |

Orbits inside a system are stored as modified equinoctial elements
(p, f, g, h, k, L): no division by zero for circular or equatorial
orbits, singular only for an exactly retrograde equatorial orbit. The
conversions both ways were rerun here and round-trip to 1e-16 (unit
mu).

### 10.4 The update step and epoch

A run advances the galaxy by the real time since its last update
(Boss, 2026-10-07 17:11Z: "once set up and configured we follow orbital
paths in real time ... Default is 1 day = 1 day"), and an option adds a
stated extra span in one go. With `dt` that time:

- Phase A: each star moves by Velocity Verlet in the galactic potential
  plus its point masses, and everything it holds is shifted by the same
  displacement.
- Phase B: each bound child's mean anomaly advances by `n dt mod 2 pi`,
  exact for any number of orbits per step. A child that is no longer
  bound (e >= 1) uses PropagateKeplerianOrbit instead.

The document puts the epoch at J2000.0 (JD 2451545.0), which also matches
the spin formulas in "Observational Kinetics for Rotational Vectors.md".

### 10.5 Collisions

Continuous detection over the step: the squared separation of two
straight paths, `A t^2 + B t + C`, against an effective radius widened by
gravitational focusing,
`R_eff = (R1 + R2) sqrt(1 + v_esc^2 / max(dv^2, floor))`. A hit inside the
step merges the pair inelastically: mass and momentum conserved,
`R_new = (R1^3 + R2^3)^(1/3)`, the smaller body deleted, and the sectors
marked for their point-mass tables to be rebuilt. GEN.110's rogue-planet
rules decide what the merged body becomes.

### 10.6 Errors in the document, and the fix the code uses

- **The universal Kepler equation counts one term twice.** The document
  writes `... + r0 chi c1(alpha chi^2) - sqrt(mu) dt` alongside
  `(1 - alpha r0) chi^3 c3`. Since `c1 = 1 - z c3` that subtracts
  `alpha r0 chi^3 c3` twice. Rerun against a high-accuracy integration,
  the document's form misses by 0.62 (ellipse), 0.10 (hyperbola) and
  0.002 (near parabola) in unit distances and breaks `f gdot - fdot g = 1`
  by up to 15%; the corrected form below lands within 1e-13. The code
  uses
  `F = s0 chi^2 c2 + (1 - alpha r0) chi^3 c3 + r0 chi - sqrt(mu) dt`,
  `F' = s0 chi (1 - z c3) + (1 - alpha r0) chi^2 c2 + r0` (= r),
  `F'' = s0 (1 - z c2) + (1 - alpha r0) chi (1 - z c3)`,
  with `s0 = r0.v0 / sqrt(mu)` and `z = alpha chi^2`. The f and g
  formulas in the document are right.
- **The hyperbolic first guess drops a sign.** Inside the logarithm the
  `sqrt(-mu/alpha)` term needs `sign(dt)` (Vallado), or backward steps
  start on the wrong branch. Halley's method usually recovers, but the
  code uses the signed form.
- **The flyby needs the speed at infinity.** The deflection uses the
  current relative speed as v_inf. When the pair is caught close in,
  the code uses `v_inf = sqrt(v^2 - 2 mu / r)` instead.
- **The merged body's position.** The document puts it at body 1's
  position at impact; the code uses the centre of mass at impact,
  `(M1 p1(t) + M2 p2(t)) / (M1 + M2)`, which momentum conservation needs.
- **The focusing floor has no unit.** `max(dv^2, 1.0)` means 1 m/s in SI
  and 1 km/s in galactic units. The code works in SI and uses 1 m/s.

## 11. Open questions for Boss

The defaults below hold until Boss decides otherwise:

- **Sector shape.** The document proposes cubic 4 pc cells keyed by a
  Morton code with a 27-cell stencil, to avoid wedge-shaped cells near
  the core. planetGen's sectors are already 4 pc but sit in rings,
  layers and slots (galaxy-coordinate-system.md), and
  galaxy-drilldown-navigation.md already turned down Morton keys.
  Default: keep the current sectors; the stencil is "this sector and
  every sector touching it", found by the existing address math. The
  physics does not change because forces are summed in galactic
  coordinates.
- **Sector frame axes.** The document's sectors are unrotated so that
  forces add without rotation. Ours rotate (+X radially outward).
  Default: physics works in galactic coordinates; the rotated sector
  frame stays for display and navigation.
- **Scaling the potential to the galaxy's shape.** The masses and scales
  above are the Milky Way's. planetGen's shape is a setting (disk scale
  length 2,800 pc, bulge scale radius 200 pc, radius 15 kpc by
  default). Default: the document's values at the default shape; when a
  galaxy's disk scale length differs, every length scales by
  `disk_scale_length / 2800 pc` and the masses stay, so the rotation
  curve keeps its shape. Galaxy-to-galaxy differences beyond that are
  phase 3+ (several galaxies).
- **Epoch and step.** The earlier default was year 0 at the galaxy's
  generation. The epoch is J2000.0, as the document and the spin
  document both use. The step is not a fixed year: Boss (2026-10-07
  17:11Z): "No 1 year per orbital update or turn, rather, once set up
  and configured we follow orbital paths in real time.  The update
  script should have an option to update for more time in 1 go if
  specified.  Default is 1 day = 1 day." Each run advances by the real
  time since the last one, plus any extra span asked for.
- **The galaxy's own gravity.** Boss (2026-10-07 12:25Z): "We'll have to
  add a galactic gravitational gradient but we need to make sure that
  it's consistent with actual science." Settled by section 10.1.

### Where positions are held (GEN.74 part 2)

- A sector's system and phenomenon entries (`SectorSystemEntry`,
  `SectorPhenomenonEntry`) each hold one `SpatialPosition3D` in light-years
  (`entry.spatial`) with the object's mass and mu; `entry.position` is its
  sector-frame Cartesian. `SpaceSector.place_in_galaxy(center_ly)` sets the
  sector's center (generation does it from the cell's position, loading
  from the stored center) and carries every entry with it.
- Planets, moons and comets (`HoldsOrbitPosition`) and stars hold one in
  AU, anchored by `galaxy/system_position.place_system`; the "system" frame is
  the offset from the body's primary (star, barycenter, wide-binary
  secondary, or parent planet for a moon). `planet.position_x/y/z` and
  `comet.position_x/y/z_au` read and write through it, so the stored km and
  mpc columns map to the object. Asteroid belts have no point and hold none.

### Stored columns and the position object (GEN.74 part 3)

| Column | Object | Unit |
|---|---|---|
| `star_systems.position_x/y/z_mpc`, sector `center_*_pc` | entry `spatial` (galactic / sector frame) | light-years, via `place_in_galaxy` |
| `planets.position_*_km`, `moons` likewise | `body.spatial` system frame | AU on the object, km in the column |
| `comets.position_*_km` | `comet.spatial` system frame | AU on the object, km in the column |
| `*.mass_kg` | `spatial.mass_kg`, with `mu = G * mass` | kg |
| stars (`stars.mass_kg`) | `star.spatial` anchored at the system | AU |

`tests/test_spatial_position_db.py` saves a placed sector, loads it, and
checks every entry, star, planet and moon sits at the same galactic place
with the same mass and mu.

### Velocity, epoch and orbit elements (Boss, 2026-10-08)

Boss (2026-10-08 23:11Z): store the orbit as part of the coordinates, with
the velocity (the immediate vector of movement) relative to the center of
the star system, then the projected course; each orbital update refreshes
the vector (wobble included) and the ellipse from it. Bodies with no closed
orbit get a path through their sector from the vector, curved by the nearby
masses (23:16Z), as a spline.

- `SpatialPosition3D` keeps a velocity in two frames, as it does a
  position: **galactic** (the sector frame has the same axes) and
  **system**, relative to the nearest star. The one last set is the truth
  and the other is derived: galactic = `star_velocity` + system, so a
  planet's 30 km/s round its star stays apart from the star's 220 km/s
  round the galaxy. `set_velocity_cartesian(vx, vy, vz, frame=)` sets it,
  `set_star_velocity` keeps the body's galactic velocity when the star's
  changes and `carry_star_velocity` keeps its system velocity. `epoch_unix`
  is when the position and velocity hold. `get_time_to_observable_movement`
  measures a "galactic" move by the galactic velocity and the others by the
  system velocity.
- **Where the velocity comes from (GEN.121).** A planet's or moon's is
  the tangent of its circular orbit at its phase (`orbit_speed_kms` long,
  `orbits.circular_orbital_velocity_au_per_year`), a comet's is
  `state_from_elements` at its true anomaly (`kepler.comet_orbital_state`
  returns it), both relative to the body's primary and put on the body when
  its position is (`update_orbital_position`, `Comet.update_orbital_state`).
  A star's, a system's and a phenomenon's is its
  `galactic_orbital_speed_kms` along the tangent of the galaxy's rotation at
  its place, counterclockwise (`system_position.galactic_velocity_ms`, which
  `place_system` and `place_in_galaxy` apply, so it follows the object into
  another sector), plus, for a runaway or hypervelocity system, its
  `runaway_speed_kms` along `runaway_direction`
  (`system_position.peculiar_velocity_ms`); a body's galactic velocity is its
  star's plus its own. `planets`, `moons` and `comets` store theirs (schema
  v60, `velocity_x/y/z_kms`) and `advance_orbital_phases` and
  `advance_comet_orbits` move it with the position; `star_systems` stores
  the system's galactic velocity (schema v61) and `advance_galactic_positions`
  turns it with the position, so the motion update still follows the
  rotation curve and the runaway velocity is carried with it, not integrated. A loaded sector's
  objects carry the epoch `store.get_orbit_epoch_unix` gives (when the
  orbits were last advanced).
- **A body with no closed orbit gets a path through its sector**
  (`physics/sector_path.py`, GEN.123). `integrate_path` follows a test
  particle from the sector entry, with the velocity it has there, against
  the sector's point masses until it leaves the cell (`sector_inside`,
  `galaxy.geometry`). Masses that cannot move it by a tenth of the
  tolerance are skipped (`relevant_masses`: a pass at impact parameter `b`
  bends velocity by about `2 mu / (b v)`), at most 32 are used, and each is
  softened like the galaxy's own (1 pc). The path is cubic Hermite spline
  knots (time, position, velocity): the ends, then the integrated sample the
  spline is furthest from, until it is within 0.2% of the sector edge or 48
  knots, so a straight crossing is two knots. The exit knot is the next
  sector's entry (`SectorPath.restart`). Tests check the bend against the
  hyperbolic deflection `2 asin(1/e)`, `e = sqrt(1 + (b v^2 / mu)^2)`, to 2%.
  Saving the knots with the sector and the job that fills them follow.
- **The orbit is derived from the vector, never stored.** A planet's, moon's
  or comet's `orbit_from_vector()` works the osculating orbit out of its
  position and velocity relative to its primary
  (`state_vectors.elements_from_state`; mu from its circular radius and
  period, or 4 pi^2 per solar mass of a comet's star), so a vector that
  carries wobble or a flyby's pull gives the orbit it is really on, and the
  generated elements (`distance`, `period`, the comet's perihelion and
  eccentricity) stay what they were generated as. `projected_orbit_au()` is
  that orbit's closed ellipse for drawing (an open orbit raises
  `ValueError`; the sector path of GEN.123 covers those).
- With the sector edge known, the position also knows its **sector
  address**, `(ring, layer, slot)`, the cell of the galaxy's sector grid
  (`galaxy/geometry.py`) it is in. It is worked out again from the galactic
  position whenever anything moves it, so it is never stale; sector entries
  and every body of their systems get the edge when the sector is placed in
  the galaxy (`place_in_galaxy`). `set_sector_address` carries the body to
  another cell's center, keeping its offset from the sector's center (in
  galactic axes). A system near a cell face can have a body in the next
  cell: the address is where the body is.
- `physics/state_vectors.py` converts between a state vector (position and
  velocity relative to the primary) and orbital elements for every conic,
  ellipse, parabola or hyperbola: `state_from_elements`,
  `elements_from_state`. The orbital update takes the new vector to the new
  osculating ellipse with the second; `closed_orbit_points` samples a
  closed one for drawing. `orbits.circular_orbital_velocity_au_per_year` is
  the velocity of the circular orbits planets and moons are on today.

### Positions at any time (MAP.70)

`physics/body_positions.positions_at(scene, years)` and its browser twin
`static/orbitpositions.js` `positionsAt` give every star, planet, moon and
comet of a `/api/systems/<id>/scene` its place `years` after the scene's
epoch: circular orbits advance their phase by 360 degrees a period, comets
follow Kepler's or Barker's equation, a close pair balances on the
barycenter by its mass fraction. `tests/test_js_unit.py` checks the two
copies against each other. `static/orbitclock.js` is the view's time
control (real time, faster, pause, back to now); the 3D view (MAP.74) wires
it in.
