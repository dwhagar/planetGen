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

`spacial-position.py`'s `SpatialPosition3D` keeps one body's position in
three frames (galactic, sector, system) and three forms (Cartesian,
cylindrical, spherical). Setting any coordinate in any frame and form
updates all of them, with the velocity, and it works out how long until
the body moves enough to be observed. It moves into the physics package
with tests, and every positioned object holds one, with its mass and
gravitational parameter (mu = G M) beside it.

Note: the prototype's minimum-observable constants (1e-4, 1e-5, 1e-6)
are placeholders; the thresholds below replace them.

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

## 10. Open questions for Boss

None of the documents settle these. The defaults below hold until Boss
decides otherwise:

- **The galaxy's own gravity.** The documents pull on a star only from
  its nearest point masses and the central black hole. The smooth mass
  of the disk, bulge and halo, which keeps the Sun at about 230 km/s, is
  missing, and without it stars drift outward. Default: a fixed analytic
  potential (a Miyamoto-Nagai disk, a Hernquist bulge and an NFW halo,
  scaled to the density model's disk and bulge in galaxy-disk-density.md) added to every
  galactic step, with the point masses as perturbations on top.
- **Frame axes.** The galactic and sector frames stay as
  galaxy-coordinate-system.md defines them: galactic +Z is galactic
  north, and slots count counterclockwise from +X; a sector's local +X
  points radially outward, +Y toward increasing angle and +Z north, from
  the sector's centre. A system frame (missing today) has its origin at
  the system's barycentre and its +Z along the primary star's spin axis;
  planets' tilts and spins are measured from it.
- **Time.** One update run advances one in-game year, from an epoch of
  year 0 at the galaxy's generation. Both are settings.

