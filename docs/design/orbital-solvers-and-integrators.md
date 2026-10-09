# Orbital solvers, integrators, thresholds and the influence search

The numerical side of the orbital update: which Kepler solver to use and where
each one stops being accurate, how to guard the universal-variable propagator,
which integrator and which libraries (and licences) are acceptable, what Boss's
movement thresholds mean in time and write load, how to find the objects inside
a Hill sphere through the ring, layer and slot sectors, a guard for every edge
case, and the errors found in Boss's uploaded documents with the correction for
each. The overall design is [orbital-updates.md](orbital-updates.md); collisions
and mergers are in [collisions-and-mergers.md](collisions-and-mergers.md); the
galaxy's smooth potential is in [galactic-potential.md](galactic-potential.md)
(see also, not edited here).

Informs: GEN.105, GEN.106, GEN.107, GEN.108, GEN.109, GEN.104, GEN.111, GEN.115, ADM.36, MAP.62, MAP.70

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result or package page, [C] computed in the
research sandbox, [R] recalled and unconfirmed. The research environment could
only read search-result text and PyPI metadata, not the papers themselves; every
[R] is on the Evidence notes list. Numbers marked [C] were reproduced while
writing this document where they could be done with numpy and the standard
library (corrected universal equation, Stumpff errors, Hill against the exact L1
distance, Wilkinson matrix conditioning, softening factors, float64 spacing,
tidal-lock and gyrochronology examples, threshold arithmetic); the rest are from
the research scripts and carry their original +-30% timing noise.

## Decisions already taken

- Movement thresholds (Boss, 2026-10-03 05:38Z): update an object only once it
  has moved at least 0.01 mpc (stars and anything star-like outside a system),
  0.01 AU (planets, companion stars) or 100,000 km (moons and satellites); store
  when each object is next due.
- Influence rule (Boss, 2026-10-09): replaces "the nearest 10 bodies at least as
  large". An object's influence radius is the Hill sphere of the largest nearby
  object; the point masses of all objects inside it plus the galaxy's smooth
  gradient (GEN.115) give the new vector of each object visited. Section 6
  applies it.
- Time step (Boss, 2026-10-07 17:11Z): follow orbits in real time, "Default is 1
  day = 1 day", with an option to advance more in one go.
- Galactic gradient (Boss, 2026-10-07 12:25Z): it must be "consistent with
  actual science".
- Sectors (Boss, 2026-10-07 17:11Z): keep the ring, layer and slot sectors with
  forces summed in galactic coordinates, "but verify the algorithm will work
  with our sector geometry." Section 6.4 is that verification.
- Edge cases (Boss, 2026-10-07 11:47Z): "we need to ensure that reasonable
  limitations for when the math breaks down at the edge cases." Section 7.
- The project is CC0 1.0 (`LICENSE.md`, `setup.py`); a GPL dependency is not
  decided (section 4.3).

## Findings in brief

- The influence set under the Hill rule is a few to ten objects in the field,
  coincidentally about the old "top 10", and only explodes near heavy objects or
  in clusters (section 6.2). A radius search over the sector grid matches brute
  force exactly (540,000 stars, 27 site and radius combinations, 0 mismatches) if
  the search radius is widened by the cell reach, about 4 pc; without that it
  misses 21.6% of true neighbours [C].
- Real-time orbits remove most daily work. A star at 220 km/s moves 0.127 AU a
  day and leaves its 4 pc sector about once per 15 kyr. The literal 0.01 mpc rule
  on the galaxy-frame velocity makes every star due every 16 days (6.2e7 writes a
  day per 1e9 stars) because rotation, not drift, moves them (section 5).
- Planets and moons are pure functions of phase. Keep the phase, compute the
  position on read (`positions_at` of MAP.70 already does), and use the
  thresholds to schedule perturbation re-fits, not row writes.
- Kepler: Mikkola or Markley (non-iterative, accurate to the float64 floor for e
  up to 0.9999) for e < 1 - 1e-5; the universal variable with Stumpff functions
  by series near e = 1 and for any state-vector step; Halley from a cubic start
  for hyperbolas (section 2).
- The universal Kepler equation in "Computational Astrodynamics.md" is wrong in
  two terms, and its closed-form Stumpff functions lose up to 4 orders of
  magnitude near z = 0 (section 1, rows 1 and 2).
- Do not step stars by days. IAS15 integrates a star plus 10 neighbours for 1 Myr
  in 101 steps; day steps would need 3.65e8 (section 4). Own code (numpy plus
  `scipy.integrate.DOP853`) needs no new dependency; REBOUND is GPL-3.0-only and
  stays out of shipped code.
- Softening of 1 pc on every point mass cancels 46% of the pull at 1.4 pc, the
  Hill radius of a solar-mass star. Soften only the central black hole.

## 1. Corrections to the source documents

Each row was checked against the source text; [C] marks values rerun here.
Corrections to the collision procedure itself are in
[collisions-and-mergers.md](collisions-and-mergers.md) section 2 and are only
listed here.

| # | Source and place | Claim | Finding |
|---|---|---|---|
| 1 | Astrodynamics, "Universal Variable Keplerian Orbit Propagation" | `F = s0 chi^2 c2 + (1 - a r0) chi^3 c3 + r0 chi c1 - sqrt(mu) dt` and `F' = ... + r0 c0` | Wrong [C]. Since `c1 = 1 - z c3`, the `r0 chi c1` term counts `a r0 chi^3 c3` twice; `F'` has `r0 c0` where `r0` is needed (`(1 - a r0) chi^2 c2 + r0 = chi^2 c2 + r0 c0`). Rerun for a = 1, e = 0.5 from periapsis: the document's form places the body wrong by 0.12, 3.2 and 6.1 (units of a) at dt = 0.5, 1.5, 3.0 and breaks `f gdot - fdot g = 1` by 1.3e-2, 2.0, 0.24; the corrected form is within 2e-15 with the identity to 4e-16. The first and second derivative of the corrected form are in section 3. The document's `F''` is right |
| 2 | Astrodynamics, Stumpff functions | closed forms for c2, c3 "in floating point implementations" | Relative error of c2, c3: 8.9e-5 and 7.8e-5 at z = 1e-12, 8e-8 and 4e-6 at 1e-10, 5e-9 and 3e-8 at 1e-8, 3e-13 and 3e-12 at 1e-4, 1e-14 at 1e-2 [C]. This is exactly the near-parabolic range the universal variable is for. Series (12 terms) below `|z| = 1` is exact to 1e-16 everywhere. Near parabolic the propagated position error is 1e-11 to 6e-8 with closed forms against 8e-14 to 2e-10 with series [C, research script] |
| 3 | Astrodynamics, first guess | thresholds `alpha > 0.15`, `abs(alpha) <= 0.15` | `alpha` is 1/length, so 0.15 depends on the unit (AU or metres). Use the dimensionless `alpha r0` or e. A better start for every case is the real root of `((1 - alpha r0)/6) chi^3 + r0 chi = sqrt(mu) dt` (Cardano): 4 Laguerre iterations for ellipses, 3 to 8 for hyperbolas [C, research script]. It divides by zero for an exactly circular start at r0 = a, so guard `1 - alpha r0 == 0` |
| 4 | Astrodynamics, hyperbolic first guess | missing `sign(dt)` | Already corrected in orbital-updates.md 10.6 |
| 5 | Astrodynamics and Full Algorithm, local Hill radius | `dist * (m / 3M)^(1/3)` ("mutual Hill radius") | Against the exact L1 distance of the restricted three-body problem [C]: Earth +0.3%, Jupiter +2.4%, Moon/Earth +6%, m/M = 0.1 +14%, equal masses +39%. With `(m / (3 (M + m)))^(1/3)` the equal-mass error is +10%. Use `a (1 - e) (m / (3 (M + m)))^(1/3)` (pericentre side), not the instantaneous distance. For a star in the galaxy use the tidal radius (section 6.1) |
| 6 | Full Algorithm, `CheckPlanetaryHillCrossings` | negative two-body energy inside the Hill sphere means a captured moon | Necessary, not sufficient. Satellites stay stable only out to roughly 0.5 R_H (prograde) and 0.9 R_H (retrograde) at e = 0, and a three-body capture is temporary without a dissipative event [R: Hamilton and Burns 1991; Domingos, Winter and Yokoyama 2006]. Use `mu = G (m1 + m2)` and a stability fraction as a setting. Roche test placement is in the collisions document |
| 7 | Full Algorithm, Astrodynamics | velocity Verlet is "symplectic"; "energy-conserving 4th-order Hermite or Yoshida symplectic" | Verlet is symplectic only for a fixed step and a smooth force. Measured (Kepler, e = 0.5) [C, research script]: step P/50 gives a bounded energy error of 4.4e-2 over 2,000 orbits (no drift); step `0.03 r^1.5` gives 1.7e-3 after 200 orbits and 2.2e-2 after 2,000 (secular drift). A top-N force list that changes membership is a discontinuous force and breaks it further. Hermite is not symplectic [R: Makino and Aarseth 1992]; Yoshida's composition is, for fixed steps only. Verlet on a circular orbit goes unstable above a step of about P/pi (dE/E = 1.1 at P/4, 11.7 at P) |
| 8 | Astrodynamics, macro step | 1-year step, "Phase A" Verlet for the star in the galaxy | Superseded by real time. Verlet is also a poor galactic integrator: after 1 Gyr the position error is 0.18 pc (0.1 Myr step), 4.6 pc (0.5 Myr), 18 pc (1 Myr), 74 pc (2 Myr). `DOP853` at rtol 1e-13 needs 5,270 evaluations for 1 Gyr (energy error 3e-14) [C, research script] |
| 9 | Astrodynamics and Full Algorithm, softening | Plummer eps = 1 pc for all point masses | The pull is reduced by `(1 + eps^2 / r^2)^(-3/2)`: 15% less at 3 pc, 46% less at 1.4 pc, 1.5% less at 10 pc [C]. That removes the region GEN.109 cares about. Hernquist, Miyamoto-Nagai and NFW are already finite at the centre. Use eps = 1 pc for the central black hole only; between stars use none (or 1e-3 to 1e-2 pc) |
| 10 | Update Algorithms, worst-case example | `Delta x_min = 0.01 pc = 3.086e14 m`, `a_min = 0.62 m/s^2` | Boss said 0.01 mpc = 3.086e11 m, 1,000 times smaller. The "0.18 AU" cell near the black hole (and 4.3e19 cells, 1 ZB) belongs to `a_min` of about 1e-6 m/s^2, not 0.62; with 0.62 the resolution is 1.7e16 m (0.55 pc) [C] |
| 11 | Update Algorithms, `a_min = 2 dx / dt^2` filter | keep sources with `mu / r^2 >= a_min` | With Boss's day step and 0.01 mpc, `a_min = 2 x 3.086e11 / 86400^2 = 83 m/s^2`; a solar-mass star then counts only inside 1.3e9 m (its radius is 7e8 m), so every neighbour is filtered out [C]. The filter only makes sense against the time to the next due date, which is years to Myr. Drop it under the Hill rule |
| 12 | The three sector definitions | Update Algorithms: 10 ly cubes; Full Algorithm and `ResolveNeighborSectors`: 11.5 ly cylindrical cells; Astrodynamics: 4 pc cubes with Morton keys | planetGen's sectors are 4 pc ring, layer and slot cells. The stencil of `resolve_neighbor_sectors` takes at most one ring and one layer each way whatever `R_search` is, so it is too small for a Hill search (spheres of 1.4 to 60 pc). Replace it with `geometry.enumerate_sectors_within_radius` plus the reach margin (section 6.4) |
| 13 | Full Algorithm, Phase B | `PropagateKeplerianOrbit(child_pos, child.velocity, dt)` about a moving star | `Entity.velocity` is global, so the child's velocity includes the star's 220 km/s. Propagating it as the star-relative velocity is wrong; subtract the star's velocity, propagate, add it back (planetGen stores the system-frame velocity apart, orbital-updates.md "Velocity, epoch and orbit elements") |
| 14 | Full Algorithm, `ExecuteInelasticMerger` | merge into the primary | The consumed body's children keep a deleted `parent_system_id`; hand them to the survivor or the host. Position, radius and the rest are in the collisions document (section 2.2) |
| 15 | Collision document, Full Algorithm, Astrodynamics | straight-line sweep, discriminant form, `R_eff` as contact, 2.44 for every body, candidate set of 10 | See [collisions-and-mergers.md](collisions-and-mergers.md) sections 2.1 to 2.4: the discriminant form is wrong for 19% to 50% of geometries at long steps, `R_eff` is a capture cross-section, 2.44 is for fluid bodies only, partners come from a mass-blind distance query |
| 16 | Kinetics, YORP | prose: obliquity driven "toward +-90 degrees"; procedure: near 10 or 170 degrees; conclusion: poles "perpendicular to orbital planes" | The document disagrees with itself. The usual statement is YORP toward 0 and 180 degrees, some toward 90 [R: Vokrouhlicky et al., Asteroids IV, 2015]; the 10 and 170 degree procedure matches it and orbital-updates.md section 6 follows it |
| 17 | Kinetics, small-body period | text `LogNormal(mu = 2.1, sigma = 0.65)` hours; table `LogNormal(13 h, 0.65)` | Median e^2.1 = 8.2 h against 13 h [C]. Pick one and make it a setting |
| 18 | Kinetics, gyrochronology | constants `n = 0.55, b = 0.31, c = 0.495, a = 0.40` credited to Barnes (2007) | They look like the Mamajek and Hillenbrand (2008) fit (0.407, 0.325, 0.495, 0.566) [R]. With those, a Sun-age star (4,600 Myr, B-V 0.65) has P = 26.3 d, not the document's 25.4 d [C] |
| 19 | Kinetics, tidal lock | worked example `tau_lock ~ 1.2 Myr`, `P_orb = 2.1 d` at 0.03 AU | With the document's formula, Q = 100, k2 = 0.3, a one-day initial spin and I = 0.33 M R^2, tau is about 80 yr, and P_orb at 0.03 AU round 1 Msun is 1.90 d [C]; "locked" holds, the numbers do not reproduce. Gladman et al. (1996) have `omega a^6 I Q / (3 G m_p^2 k2 R^5)`, so 4/9 differs by a convention-dependent factor [R] |
| 20 | Kinetics, citations | merger-remnant spin 0.69 cited to Gilbert et al. 2020 and Gladman et al. 1996 | A planet-tides paper and a tidal-synchronisation paper cannot support a black hole spin; Lousto and Zlochower (2012), also listed, is the fitting source. The list was not checked against the papers |
| 21 | Foundations, Wilkinson matrix | "well-conditioned ... `kappa_inf(W_n) ~ n 2^(n-1)` is modest" | Self-contradictory: `kappa_inf(W_n) = n` [C]. Growth 2^63 and the 1.31e5 backward error figure are right |
| 22 | Foundations, summary table | "Reliability Bound Threshold" cells for Lanczos, tridiagonal QR and bidiagonal SVD end mid-sentence | Truncated in the source; nothing to carry |

Also checked and correct: the Lagrange f and g formulas; `F''` and the corrected
`F` and `F'`; `n = sqrt(mu / a^3)`; `mu = 4 pi^2 M` in AU, year, solar-mass
units; the fluid Roche coefficient 2.44 (rigid about 1.26); the turning angle
`2 asin(1/e)`; the Plummer force law; the rotation curve numbers (229.27 km/s at
8.128 kpc); the ellipse bracket `[M - e, M + e]`; the Foundations formulas for
Newton, Brent, BFGS and LU growth. Roche, Skumanich and breakup formulas in the
Kinetics document check: `v_crit = sqrt(2 G M / (3 R_p))`, and
`P_crit = sqrt(3 pi / (G rho))` is 2.33 h at 2,000 kg/m^3 [C].

## 2. Kepler solvers

Test: 11 eccentricities from 0 to 0.9999, 300 mean anomalies each (log-spaced
down to 1e-12, uniform, near 2 pi, and the edge values), error against mpmath at
50 digits; hyperbolic e from 1.0001 to 10, |M| from 1e-10 to 1e12 [C, research
scripts `kep_solvers.py`, `ell_accuracy.py`, `hyp_accuracy.py`].

**Accuracy (ellipse).** Every method that converged (Newton, Halley, Danby,
Laguerre-Conway, Markley, the project's scipy code) has maximum |E error| of about
5e-16 up to e = 0.3, 1e-15 at 0.7, 2.6e-15 at 0.9, 2.5e-14 at 0.99, 2.5e-13 at
0.999 and 2.4e-12 at 0.9999. That growth is the conditioning of Kepler's equation
(`dE = dM / (1 - e cos E)`), not a solver defect. In position terms the E-based
formulas lose digits as e nears 1 (`cos E - e` cancels): relative error 3e-14 at
1 - e = 1e-3, 7e-12 at 1e-5, 5e-10 at 1e-7. The universal variable with series
Stumpff functions gives 8e-14 at 1 - e = 1e-5 and 2e-10 at 1e-12. Mikkola's
starter alone is not enough: one Halley polish leaves up to 9e-9, two reach 9e-16.

**Cost (numpy, 1e6 orbits a call, fixed iteration count reaching residual
5e-15; +-30% noise).**

| Method | e <= 0.3 | 0.3 to 0.9 | 0.9 to 0.9999 |
|---|---|---|---|
| Newton / Halley / Danby (iterations) | 4 / 3 / 2 | 8 / 4 / 3 | 13 / 7 / 5 |
| Markley, one pass (us per orbit) | 0.5 to 0.7 | 0.6 | 0.5 |
| Mikkola plus 2 Halley (us per orbit) | 0.3 to 0.5 | 0.5 | 0.3 to 0.4 |
| Current code (Halley plus `brentq` fallback, us per orbit) | 0.6 | 1.5 | 11 |
| Universal variable, Laguerre x4 (us per orbit) | 3 (closed Stumpff) or 7 (series) | | |

Mikkola plus 2 and Markley are fastest with a fixed operation count (no
data-dependent loop); the current code is 20 to 30 times slower for comets (19 per
million fall to `brentq`); the universal variable is 4 to 10 times slower but
returns the whole state.

**Hyperbolic.** Every converging method reaches 1.5e-15 absolute error in H;
relative error is 1.9e-12 at e = 1.0001, 1.6e-13 at 1.001 and under 2e-15 from
1.1. The textbook start `H0 = M / (e - 1)` needs 32 to 40 iterations for
e <= 1.5 and misses 5 of 240 cases at 1.0001. The cubic or log start (the real
root of `(e/6) H^3 + (e - 1) H = M`, or `sign(M) ln(2 |M| / e + 1.8)` when smaller)
with Halley needs 3 to 4. Newton on the universal equation failed for some very
large dt at e = 1.0001; use Laguerre. Gooding and Odell (1988) is the classic
reference [R]; not reimplemented.

**Switch points.**

| Regime | Method |
|---|---|
| 0 <= e < 1 - 1e-5 | Mikkola cubic start plus 2 Halley, or Markley plus one Newton polish; reduce M to [0, pi] by symmetry; vectorise, no loops |
| 1 - 1e-5 <= e <= 1 + 1e-5, and any step from a state vector | universal variable (corrected equation, series Stumpff, Laguerre); works from a position and velocity, so it covers perturbed orbits |
| abs(e - 1) < 1e-12 | Barker's closed form (`solve_barker_equation`, already in `kepler.py`) |
| e > 1 + 1e-5 | Halley (or Danby) from the cubic and log start; any abs(M); H stays under about 30, so sinh cannot overflow |

Keep `brentq` as the last-resort fallback, as `kepler.py` does. The
non-iterative solvers are easy to mirror in `static/orbitpositions.js`, which
`tests/test_js_unit.py` keeps in step with the Python twin (MAP.70, MAP.62).
Swapping `solve_eccentric_anomalies` is optional and low priority.

## 3. The universal-variable propagator and its guards

The corrected equation (section 1, rows 1 to 3), with `s0 = r0.v0 / sqrt(mu)`,
`alpha = 2 / r0 - v0^2 / mu` and `z = alpha chi^2`:

```
F   = s0 chi^2 c2 + (1 - alpha r0) chi^3 c3 + r0 chi - sqrt(mu) dt
F'  = s0 chi (1 - z c3) + (1 - alpha r0) chi^2 c2 + r0          (= r)
F'' = s0 (1 - z c2) + (1 - alpha r0) chi (1 - z c3)
```

Then `f = 1 - chi^2 c2 / r0`, `g = dt - chi^3 c3 / sqrt(mu)`,
`r = f r0 + g v0`, `fdot = sqrt(mu) chi (z c3 - 1) / (r r0)`,
`gdot = 1 - chi^2 c2 / r`, `v = fdot r0 + gdot v0`, as in the source.

Guards:

- Stumpff `c0` to `c3` by their series (12 terms) for `|z| < 1`, the closed
  forms for `|z| >= 1`.
- Start from the Cardano root of the cubic above, not the `alpha > 0.15` split;
  skip the cubic's division when `1 - alpha r0 == 0` (circular at r0 = a) and use
  `chi = sqrt(mu) dt / r0`.
- Iterate with Laguerre (Newton can stall for very large dt near e = 1) to a
  relative tolerance of 1e-13, cap at 50 iterations, and on failure bisect on a
  bracket and log it.
- `dt = 0` returns the input state; negative `dt` is allowed (tests cover both).
- Check `f gdot - fdot g = 1` to 1e-12 in tests, not in production.
- Hyperbolic deflection uses `v_inf = sqrt(v^2 - 2 mu / r)` when the pair is
  caught close in (orbital-updates.md 10.6).

## 4. Integrators and libraries

### 4.1 Measurements

REBOUND 5.2.2 and scipy `DOP853` in the research sandbox [C, research scripts
`rebound_test2.py`, `dop_galaxy.py`]:

| Case | Result |
|---|---|
| Sun plus 5 planets, 1,000 yr | WHFast at P_min/20: 83,074 steps, 0.55 s, dE/E 5.7e-11. WHFast at 1 day: 365,251 steps, 2.4 s, 2.9e-12. IAS15: 149,702 steps, 11 s, 1e-15 |
| WHFast throughput | N = 2: 1.1 us a step; N = 11: 10.6 us; N = 101: 294 us |
| Sun-like star plus 10 neighbours at 3 pc, 30 km/s, 1 Myr | IAS15: 101 steps, 0.02 s, dE/E 2e-16. WHFast or Mercurius at 1e-3 Myr: 1,000 steps, 2e-9 |
| Same 11 bodies in the project's galactic potential (own numpy force, `DOP853`) | 1 Myr: 197 evaluations, 0.07 s, dE/E 1e-14; 100 Myr: 785 evaluations, 0.29 s, 9e-14 |
| Same at a day step | 3.65e8 steps a Myr x 10.6 us, about 3,900 s a Myr per object |
| Two Jupiter masses 0.03 AU apart (inside their mutual Hill sphere), 200 yr | WHFast dE/E 0.25 (fails), IAS15 7e-3, Mercurius 2.8e-3 in 11 s. Point masses without radius: an encounter needs a physical radius and a stop rule |

### 4.2 What each is for

- **Velocity Verlet or leapfrog:** fixed step and a smooth force only; not for the
  galactic orbit (section 1, row 8), variable steps, or steps near the period.
- **Wisdom-Holman (WHFast):** hierarchical systems with a dominant mass; step
  still limited by the shortest period and encounters (about P_min/20); not for a
  star among similar neighbours. **MERCURIUS** switches to IAS15 inside a multiple
  of the Hill radius; slow in a hard encounter.
- **IAS15:** adaptive, 15th order, round-off energy error in most cases; best for
  the galactic influence set and for encounters.
- **Hermite (4th order):** needs the force derivative, not symplectic;
  `DOP853` does the same job with no extra code.
- **`scipy.integrate.DOP853`:** 1e-14 for 1 Myr in 197 evaluations, in a
  dependency the project already pins (scipy 1.13.1 on Python 3.9,
  `requirements.lock`); not symplectic, but at rtol 1e-13 the 1 Gyr drift is
  3e-14. This is the recommended integrator.

Perturbation by the influence set is done per due object, over the whole interval
since its last update, with one adaptive call, never in day steps: the clock
advances a day a day, the integration happens when an object is due. A pair
inside a few mutual Hill radii, or within 3 R_hill of a planet, leaves the cheap
path: sub-step at IAS15-class accuracy and stop at contact or Roche
(collisions-and-mergers.md section 5).

### 4.3 Libraries and the licence decision

PyPI JSON API, 2026-10-09 [S]:

| Package | Latest | Licence | Verdict |
|---|---|---|---|
| rebound | 5.2.2 | GPL-3.0-only | wheels for Python >= 3.9 on all platforms; the licence blocks bundling |
| reboundx | 5.1.0 | GPL-3.0-only | sdist only (needs a C compiler); no |
| poliastro | 0.17.0 (2022) | MIT | abandoned; fork `hapsira` 0.18.0 (2023) MIT |
| astropy | 8.0.1 | BSD-3 | needs Python >= 3.11 (6.0.1 is the last for 3.9); not needed |
| galpy | 1.12.0 | BSD | Milky Way potentials, but needs Python >= 3.10 and is heavy for one formula set |
| gala | 1.12.0 | MIT | Python >= 3.12 and no Windows wheel |
| heyoka 7.13.2, pykep 3.0.1 | | MPL-2.0 | manylinux wheels only; no |
| pyorb | 0.6.3 | MIT | tiny element converter; `state_vectors.py` already does it |
| numba | 0.68.0 | BSD | optional speed-up if numpy becomes the bottleneck |

Licence reasoning (not legal advice [R]). The project is CC0 1.0. CC0 code can
sit inside a GPL work, but a project that imports or bundles REBOUND and
distributes the result is, on the FSF reading, a combined work to be offered under
GPL-3.0, which ends CC0's "anyone can reuse this for any purpose" for that
distribution; whether a dynamic import creates a derivative work is disputed. Own
code avoids the question and is small. Recommendation: do not add REBOUND to
`requirements*.lock` or the installers; if Boss wants a cross-check, one optional
`pytest.importorskip("rebound")` test that is not shipped or installed by
default, noted in the README.

## 5. Movement thresholds and scheduling (GEN.106)

Status: GEN.106 is built (schema v66, PR #802): each moving row has an
`epoch_unix` and an indexed `next_update_due`, and a run moves only rows past
their threshold ([orbital-updates.md](orbital-updates.md) section 3). The
recommendations here that are not built: letting a star follow its analytic
orbit between rewrites, computing planets and moons from their phase on read,
and defining the moon rule on path length.

**Thresholds in time** (`THRESHOLDS_M` in `physics/position.py`; 1 km/s = 1.0227
pc/Myr; 0.01 mpc = 3.086e11 m = 2.06 AU [C]):

| Speed (km/s) | 1 | 5 | 20 | 30 | 50 | 100 | 220 (galactic rotation) |
|---|---|---|---|---|---|---|---|
| Time to move 0.01 mpc | 9.8 yr | 1.96 yr | 179 d | 119 d | 71 d | 36 d | 16.2 d |

A star at 220 km/s moves 0.127 AU in a day, 46 AU in a year, 0.0225 pc (0.56% of
a sector) in 100 years and 225 pc in 1 Myr.

**Time to leave a 4 pc sector** (straight line from a uniform start in a 4 pc
cube, Gaussian velocities; dispersions from Bensby et al. 2003 as recalled [R]:
thin disk (35, 20, 16) km/s, thick (67, 38, 35), halo (160, 90, 90)) [C, research
script `thresholds.py`]:

| Population and frame | Median first exit | 10th percentile | Crossings per Myr per star |
|---|---|---|---|
| thin disk, peculiar motion only | 42 kyr | 7.1 kyr | 14.5 |
| thin disk, inertial (+220 km/s) | 8.1 kyr | 1.5 kyr | 66.7 |
| thick disk, inertial | 7.3 kyr | 1.3 kyr | 66.8 |
| halo, no rotation | 8.9 kyr | 1.5 kyr | 69.3 |

The grid is fixed in the inertial frame, so rotation dominates. With 1e9 stars
that is 182 sector changes a day (1,825 for 1e10), against 6.2e7 writes a day from
the literal 0.01 mpc rule at 220 km/s (8.4e6 at 30 km/s).

**The code today.** See the head of [orbital-updates.md](orbital-updates.md).
`advance_galactic_positions` moves every placed system, phenomenon and
stand-alone facility along a circular galactic orbit by one angle and refiles it
when the cell changes (peculiar and runaway velocity carried, not integrated);
`advance_orbital_phases` rewrites a planet's or moon's phase and position in one
set-based UPDATE against one global reference time
(`orbit_simulation_state.last_updated_at`, `get_orbit_epoch_unix`). Each body has
a `min_update_interval_years` floor (where a phase delta falls below float
resolution), which is not a `next_update_due`. There is no `next_update_due`
column and no per-object epoch in the database (`SpatialPosition3D.epoch_unix`
exists in memory), so with objects updated at different times a global timestamp
loses time for the rows that were skipped.

**Planets and moons** (1 Msun host, circular):

| a (AU) | Period | v (km/s) | Time to move 0.01 AU | Updates a day |
|---|---|---|---|---|
| 0.05 | 4.1 d | 133 | 0.13 d | 7.7 |
| 1 | 365 d | 29.8 | 0.58 d | 1.7 |
| 5.2 | 11.9 yr | 13.1 | 1.3 d | 0.75 |
| 30 | 164 yr | 5.4 | 3.2 d | 0.31 |

Log-uniform a from 0.05 to 50 AU gives 2.2 updates a day per planet: at two
planets a star and 1e9 stars, about 4e9 writes a day, 50,000 a second sustained.
Moons at 100,000 km: the Moon every 27 h, Callisto 3.4 h, Titan 5 h, Io 1.6 h;
Phobos (2a = 18,800 km) and Deimos (47,000 km) never move 100,000 km from their
last position, so the moon rule must be on path length or a fraction of the orbit,
never on displacement.

**Cost.** The arithmetic is cheap (0.3 to 0.7 us per orbit; 1e9 objects is 6 to 12
minutes of numpy on a core); the cost is I/O. A single MariaDB writer manages
perhaps 1e4 to 1e5 updated rows a second with batching [R: an assumption, not
measured]. At 1e5 rows a second, 1e9 rows take 2.8 hours and 4e9 planet and moon
writes a day take 11 hours. The literal thresholds are not sustainable at 1e9 to
1e10 objects with a row write per update; they are if bound bodies are not written.

**What "visible" means.** Assuming 1,000 pixels across a view (the map code's zoom
limits were not read): the 30 kpc galaxy view shows 30 pc a pixel (0.01 mpc is
3e-7 pixel); a 4 pc sector view 4 mpc a pixel; a 100 AU system view 0.1 AU a pixel
(0.01 AU is 0.1 pixel); a 4 AU inner-system view 0.004 AU a pixel (2.5 pixels); a
Jovian system of 5e6 km 5,000 km a pixel (100,000 km is 20 pixels). The system and
planetary thresholds are about a pixel where those objects are drawn; the galactic
one is three orders below any pixel. Half a pixel at the 4 pc view is about 2 mpc:
a 220 km/s star due every 9 years, 3.7e5 writes a day per 1e9 stars.

**Recommended scheme.**

1. Three kinds of object. Analytic (planets, moons, comets, close binaries on
   Keplerian orbits; stars on the smooth galactic orbit): store elements, phase and
   epoch, compute positions on read. Sector-bound (stars): rewrite the DB position
   at sector exit or when the influence set changes. Perturbed (a non-empty
   influence set): re-fit at the due time.
2. `next_update_due` is the earliest of (a) exit from the current sector (a
   straight line over a horizon, `sectors_along_segment`, bisected in time; a
   face-distance bound works but needed a median of 16 re-checks, 39 at the 90th
   percentile, no violations in 348 cases [C]); (b) the time for the influence set to
   change, `R_infl / (10 v_rel)` (about 7 kyr for R = 3 pc, v = 40 km/s); (c) one
   orbital period for bound bodies with a perturber; (d) the literal threshold time
   if Boss keeps it.
3. An indexed `next_update_due` needs no bucketing. Parallel workers partition by
   due-day and sector id (ring range) so none shares a sector table.
4. Keplerian propagation has no step limit (the phase wraps exactly). Re-fit
   osculating elements often enough that the secular change between re-fits is
   under the threshold, roughly every `P / (m_pert / M)^(1/2)` for planets [R:
   untested heuristic]. N-body steps are at most P_min/20 (WHFast) or IAS15's own
   control; a pair inside 3 to 5 mutual Hill radii leaves the cheap path.
5. Star positions feed the nearest-systems search, sector tables and the map, and
   change by 0.0067 sectors a century: update them lazily at sector exit.

## 6. The influence search under the Hill-radius rule (GEN.109)

### 6.1 Which Hill radius

For a star in the galaxy the primary is the galaxy, so the sphere is the tidal
(Jacobi) radius `r_J = (G m / (4 Omega^2 - kappa^2))^(1/3)` [R: Binney and
Tremaine, *Galactic Dynamics*, 2nd ed.], using the rotation curve of GEN.115.
For a flat curve (`kappa^2 = 2 Omega^2`) it is `(G m / (2 Omega^2))^(1/3)`, the
form in collisions-and-mergers.md section 3.1. It is about 12% larger than the
point-mass Hill radius `R (m / 3 M_enc)^(1/3)` at 8 kpc. Radii in pc, point-mass Hill / Jacobi, midplane, GEN.115 potential [C, research
script `hill_galactic.py`]:

| m (Msun) | R = 0.5 kpc | 2 kpc | 8.128 kpc | 15 kpc |
|---|---|---|---|---|
| 1 | 0.24 / 0.29 | 0.54 / 0.68 | 1.22 / 1.37 | 1.89 / 2.10 |
| 10 | 0.52 / 0.62 | 1.16 / 1.47 | 2.62 / 2.95 | 4.07 / 4.53 |
| 100 | 1.13 / 1.33 | 2.51 / 3.16 | 5.65 / 6.36 | 8.76 / 9.76 |
| 1e4 | 5.24 / 6.19 | 11.7 / 14.7 | 26.2 / 29.5 | 40.7 / 45.3 |

`r_J(m, R) = k(R) m^(1/3)` with `k(R)` depending only on the galactocentric
radius, so `k` is one number per ring (about 4,000 rings in the default
galaxy). Inside a system the rule runs per level: the host star's sphere (1 to 3
pc) holds the whole planetary system, so members see neighbouring stars only
through the host's motion and the galactic tide; planets and moons use their
siblings' mutual Hill radii and the host (`a (m / 3 (M + m))^(1/3)`).

### 6.2 Cost of the sphere

Sectors touched (about volume / 64 pc^3 plus a one-cell shell) and objects
inside, for the largest nearby object of the given mass at 8.128 kpc [C]:

| Largest m (Msun) | Radius (pc) | Sectors | Objects at 0.1 pc^-3 | at 1 pc^-3 | at 100 pc^-3 (cluster) |
|---|---|---|---|---|---|
| 1 | 1.4 | 8 | 1 | 11 | 1.1e3 |
| 10 | 3.0 | 18 | 11 | 108 | 1.1e4 |
| 100 | 6.4 | 63 | 108 | 1.1e3 | 1.1e5 |
| 1e4 | 29.5 | 2,359 | 1.1e4 | 1.1e5 | 1.1e7 |

In the field (largest neighbour up to 10 Msun at solar density, 0.1 pc^-3) the
set is a few to ten objects. The cost explodes only near a heavy object (an
intermediate-mass black hole, a cluster treated as one mass) or at cluster
density.

### 6.3 Finding the largest object, cover lists and the cap

- Add `max_mass_kg` and `max_mass_id` to each sector's point-mass table (the
  Full Algorithm already keeps `max_mass`). Take the maximum over the cells
  within a first radius of 1.5 edges (6 pc covers everything with a Hill radius
  under 6 pc), compute `R = k(ring) M^(1/3)`, and grow the scan until the largest
  mass found no longer enlarges R (a monotone fixed point; a few passes).
- **Gap.** A heavy object beyond the scanned radius can still hold the target in
  its own sphere (1e3 to 1e4 Msun spheres are 14 to 30 pc; the target scans 6
  pc). In a test with 400,000 objects (Salpeter 0.1 to 100 Msun plus 40 of 300 to
  1e4 Msun, five sites including core, rim and a slot-0 seam) the fixed point
  alone missed 74 of 243 pairs where the target lay inside another object's
  Hill sphere (30%). Registering every object whose sphere is wider than the
  first scan (73 objects) in a per-sector cover list of the sectors it touches
  (23,467 rows, built in 0.4 s) and adding the covering objects to the target's
  set gave 0 misses [C]. This is the `hill_occupants` / `r_hill_max` idea of the
  Full Algorithm; only objects with a Hill radius above 1.5 edges need it.
  Maintain the lists when such an object moves (a sector change or a fraction of
  its Hill radius) or is created or destroyed. If Boss prefers the simpler rule of
  clipping the radius to the sector stencil (about 4 pc, as collisions-and-
  mergers.md section 3.1 suggests), the cover lists are the improvement to make
  later; the galactic gradient term covers the rest either way.
- **Cap.** In the same test the radius had median 3.1 pc, 95th percentile 27 pc,
  maximum 29 pc; members median 12, 95th percentile 4,936, maximum 8,845 (all from
  the injected cluster-scale masses). Pairwise sums over thousands of members for
  each of thousands of targets is N^2 work. Setting `max_exact_members`, default
  64 (the research draft proposed 256; 64 is the collisions document's number):
  the nearest are summed exactly, the rest by one monopole per sector (the
  sector's total mass at its centre of mass), and the truncation is logged.
  Plummer softening of 1 pc is kept for the central black hole only. Skipping a
  contribution whose acceleration moves the object by less than the threshold over
  the interval to the next due time (`a < 2 dx / T^2`, T that interval, not one
  day) is a valid optimisation under the new rule. Members bound in a cluster
  that is itself an object should count as members of that object: a design
  question for Boss.
- Collision candidates are not the influence set. They come from a mass-blind
  distance query of radius `R1 + R2 + |dv| dt` (collisions-and-mergers.md
  section 2.4).

### 6.4 Neighbour search through the sector bounds (verification of the sector geometry)

`geometry.enumerate_sectors_within_radius(center, radius, edge)` returns cells
whose centre lies within the radius, not cells that intersect the sphere. A query
must therefore search with `radius + reach`, where reach is the largest distance
from a cell's centre (`sector_position_pc`) to any point of the cell: 4.00 pc at
ring 0, 3.71 at ring 1, 3.64 at ring 2, 3.47 at ring 5, 3.49 at ring 10, 3.48 at
ring 100, 3.49 at ring 2,000 and 3.48 at ring 3,700 (edge 4 pc; 4.2 used to be
safe). The query reads the point-mass tables of those cells and keeps members
within R. With the margin the results were identical to numpy brute force and a
`cKDTree` at: the core (R under 60 pc), the ring-0 and ring-1 axis, mid radius
(8,128 pc), the rim (14,900 pc), either side of the slot-0 seam (y = +-0.01 pc),
top and bottom layers (z = +-600 pc) and a point on a ring face (4,000 pc), at
query radii 2, 6 and 20 pc, 15 queries each: 540,000 stars, 0 mismatches [C,
research script `geom_search_test.py`]. The pure-Python grid search took 0.4 to
22 ms a query (11 to 907 cells) against 0.14 to 8.8 ms for a prebuilt `cKDTree`
and 50 to 140 ms for brute force; the tree is for tests only. In production the
table reads come from the database per sector, so cost is cell count times the
read.

**Test plan for the TODO's "same nearest bodies as brute force" test.**

1. Seeded synthetic catalogue at the sites above (plus a ring-0 pie wedge and
   points exactly on slot and ring faces), Salpeter masses and a few heavy
   objects.
2. Grid radius search against numpy brute force and `cKDTree` for radii 2, 6, 20
   and 60 pc: identical sets.
3. The full rule (fixed point plus cover lists) against the exact pairwise
   criterion "target inside j's Hill sphere": zero misses.
4. Negative control: the same search without the reach margin must fail (21.6% of
   neighbours missed here), so the test cannot pass vacuously.
5. Cap test: truncating to `max_exact_members` keeps the total force within a
   stated tolerance of the full sum in the dense case.
6. Reproducibility: the same seed gives the same influence sets in any worker
   order.

## 7. Edge cases and guards (GEN.108)

| Edge case | What happens | Guard |
|---|---|---|
| Near-parabolic (abs(e - 1) < 1e-5) | E-based Kepler loses digits (7e-12 relative position error at 1 - e = 1e-5, 5e-10 at 1e-7) | universal variable with series Stumpff; Barker below 1e-12; the round trip at e = 1 +- 1e-12 already works in `state_vectors.py` |
| e crosses 1 after a perturbation | A 4e-12 relative kick at periapsis flips ellipse to hyperbola (verified) | classify from the state vector each time; log "orbit became unbound or bound"; no hysteresis |
| Step longer than an orbit | Verlet and any integration fail (unstable above P/pi) | analytic phase wrap with `fmod(t, P) / P`, exact in IEEE arithmetic: for t = 1e9 d and P = 0.319 d the error is 0 rad against 1.7e-6 rad for `(n t) mod 2 pi` [C] |
| Phase accumulation over 1e9 days | wrapped daily additions differ from extended precision by 3.5e-11 rad after 1e6 steps (Io-like); rounding of `n` alone gives 4e-7 rad over 1e9 days (Io), 2e-9 (Earth), 1e-11 (Neptune); a FLOAT (24-bit) column would give hundreds of radians | DOUBLE period and phase columns (check schema types); the `fmod` form; refit `n` from the semi-major axis when it changes |
| Float64 at galactic distance | coordinate spacing 32.8 km at 8.128 kpc, 65.5 km at 15 kpc, 131 km at 30 kpc, 4 m at 1 pc, 3e-5 m at 1 AU [C]; 0.01 mpc is 9.4e6 spacings at 8 kpc, 100,000 km only 3,052; adding one 1.9e10 m step repeatedly to an absolute coordinate biases it 4 km a step | sector-local and system-local offsets (as `position.py` does); never accumulate in absolute coordinates |
| Retrograde equatorial; radial (h = 0) | equinoctial elements singular at i = pi; `elements_from_state` raises for h = 0 | propagate from state vectors (universal f, g); elements for display; catch h = 0 and treat as rectilinear or add a tiny tangential component and log |
| Backward or zero step | first guesses need `sign(dt)` | test dt < 0 and dt = 0 |
| Hill sphere against sphere of influence | SOI / Hill = 0.62 (Earth), 0.91 (Jupiter), 0.75 (Neptune), 0.46 (Moon) | Hill (or Jacobi) for membership and warnings, SOI only for patched-conic estimates; never mix |
| Capture inside a Hill sphere | two-body negative energy is not a stable capture | stability fractions as settings (0.5 prograde, 0.9 retrograde of R_H [R]) and a dissipation event |
| Roche limit; close pass in a long step | break-up inside the limit; a straight sweep misses curved paths | factor 2.44 fluid, 1.26 rigid, tested on every update for pairs already inside; substep any pair within a few mutual Hill radii; window caps in collisions-and-mergers.md sections 2.3 and 5.3 |
| Galactic centre | point-mass force diverges; Hernquist, Miyamoto-Nagai and NFW are finite | soften only the central black hole (eps = 1 pc) |
| Beyond the galaxy bounds | flung out | keep and flag |
| Mass change (merger) | mu, Hill radius and cover lists change | recompute mu and `max_mass`, rebuild the affected sector tables and cover lists |
| Skipped rows lose time | global `last_updated_at` | per-object epoch column (GEN.106) |
| Non-finite values | NaN, inf, e < 0 | `finite_domain` decorators and `ValueError` raises; test each |
| Ill-conditioned solves, trajectory fits, rotation drift | element growth; BFGS curvature turns non-positive; round-off leaves rotation matrices non-orthonormal | LU with partial pivoting, then rook pivoting or Householder QR (growth 1), then equilibrate and mpmath; Wolfe line search, Powell's damped update, reset to steepest descent; re-orthonormalise (SVD or quaternion) on every save |

## 8. Hill-radius warnings and the admin email (GEN.111)

- No library needed: standard `smtplib` with `email.message.EmailMessage`,
  `ssl.create_default_context()` for STARTTLS (587) or implicit TLS (465) and a
  timeout, matching the USR.3 plan. Send as a queued RQ job so a slow SMTP server
  never delays an update; log failures (never the password), retry with backoff;
  with SMTP unset the warning goes to the log only.
- "Inside another's Hill radius" is a state, true for every planet and moon and,
  in the field, for roughly a quarter to a third of stars [R: estimate]. New
  entries alone: `n pi r_H^2 v` = 0.1 x pi x 1.37^2 x 41 pc/Myr (40 km/s) = 24 per
  Myr per star, about 3.3e4 pair entries a day per 1e9 stars. Rules: exclude
  hierarchical members; trigger on entry only, and require a gravitational
  criterion (bound or capture candidate, or a straight-line deflection of about 1
  degree; collisions-and-mergers.md section 3.3; at 1 pc the escape speed from a
  solar mass is 0.3 km/s, so a 40 km/s flyby is never bound); one row per pair in a
  `hill_warnings` table (first seen, last seen, last emailed), re-sent after a
  cooldown (7 days) or on leaving and re-entering; one digest per run of up to 25
  pairs with the two IDs and positions, at most 5 emails a day; every event also
  goes to the activity log and the table.

## 9. Spin and tilt: decisions for GEN.104

GEN.104 is not in the code yet: there is no spin vector; the schema holds a
static `rotation_period_hours` on planets and moons, `spin` on black holes and
`spin_period_ms` on neutron stars. The cascade of
"Observational Kinetics for Rotational Vectors.md" is carried into
orbital-updates.md section 6 unchanged apart from these choices:

- Asteroid obliquity: the procedure's two Gaussians (10 and 170 degrees, sigma 8,
  probability 0.5 each) stand; drop the prose "+-90 degrees" and the conclusion's
  "perpendicular to orbital planes" (section 1, row 16).
- Small-body period: one distribution in settings (median 8.2 h from `mu = 2.1`,
  or 13 h), floor 2.2 h; `P_crit = sqrt(3 pi / (G rho))` gives the density-aware
  floor (2.33 h at 2,000 kg/m^3).
- Gyrochronology constants: take the Mamajek and Hillenbrand set above after the
  [R] check, or Barnes's own; test that a solar-age, B-V 0.65 star gives 25 to 27
  days.
- Tidal locking: confirm the coefficient (4/9 or 1/3) and recompute the worked
  example (section 1, row 19) before the "system age past locking time" test is
  written; the test should use the same function for the example and the code.
- Rotation axes: the axis is the orbit normal tilted by the obliquity at a random
  precession angle, as in the source; the pole's right ascension and declination
  are measured in the galactic frame, matching orbital-updates.md section 10.2.

## 10. What is carried from each uploaded document

| Source | Carried | Not carried |
|---|---|---|
| Orbital Update Full Algorithm.md | data structure, passes, `hill_occupants` and `r_hill_max`, dirty flags, atomic commit (orbital-updates.md section 4) | top-10 influencers, nightly year step, `ResolveNeighborSectors`, Phase B velocity (rows 12, 13) |
| Orbital Position and Vector Update Algorithms.md | the per-sector point-mass table (128 B a row; 551 rows is 70.5 KB against 1 ZB for a uniform grid and 32 MB for an octree), the supermassive registry, `max_child_mass` pruning, a temporary in-memory tree, "never store vector fields" | the `a_min` filter and octree depth formula (moot without stored grids), the cylindrical stencil, the 0.01 pc example (rows 10 to 12) |
| Orbital Position and Vector Mathematical Foundations.md | Brent bracket and tolerance, BFGS safeguards, LU growth and the QR fallback, mixed-precision refinement (section 7) | Lanczos, dense QR and SVD: no large sparse matrix problem, numpy supplies the 3 by 3 decompositions |
| Computational Astrodynamics.md | the potential (galactic-potential.md), frames, modified equinoctial elements, hyperbolic deflection, the micro-pass idea, phase wrapping, J2000 | the uncorrected universal equation and Stumpff forms, the year step, Morton keys, 1 pc softening of every mass |
| Observational Kinetics for Rotational Vectors.md | the cascade of orbital-updates.md section 6 | the inconsistencies in section 9 |

## Evidence notes

[C] results are in the research scripts (scratchpad, not in the repo)
`ell_accuracy.py`, `ell_speed2.py`, `hyp_accuracy.py`, `near_parabolic*.py`,
`univ_*.py`, `stumpff_test.py`, `verlet_*.py`, `rebound_test2.py`,
`dop_galaxy.py`, `hill_*.py`, `thresholds.py`, `precision.py`,
`geom_search_test.py`, `influence_radius_test.py`, `exit_bound.py`, and in the
reruns named in the evidence paragraph at the top.

[R] items to verify once paper access is allowed: the velocity dispersions
(Bensby et al. 2003); local stellar density 0.1 pc^-3 (RECONS); the Jacobi radius
formula (Binney and Tremaine 2008); satellite stability limits (Hamilton and Burns
1991; Domingos, Winter and Yokoyama 2006); YORP obliquity behaviour (Asteroids IV);
the tidal-locking coefficient (Gladman et al. 1996); Hermite not symplectic
(Makino and Aarseth 1992); Gooding and Odell (1988); that poliastro was archived;
the REBOUND integrator papers (Rein and Liu 2012; Rein and Tamayo 2015; Rein and
Spiegel 2015; Rein et al. 2019); Markley (1995) and Mikkola (1987), implemented
from memory and validated against mpmath, so the numbers stand even if a
coefficient differs; the gyrochronology constants' origin; the Kinetics reference
list; the FSF position on GPL and linking; the MariaDB write rate; the stellar
fraction inside a neighbour's Hill sphere.

## Sources

- PyPI JSON API, 2026-10-09: `pypi.org/pypi/{rebound,reboundx,poliastro,hapsira,astropy,galpy,gala,heyoka,pyorb,pykep,numba,scipy}/json` (versions, licences, `requires_python`, wheel tags).
- Run in the research sandbox: rebound 5.2.2, numpy 2.5.3, scipy 1.18.1, mpmath 1.4.1.
- Repo files read: `LICENSE.md`, `setup.py`, `requirements.lock`, `docs/TODO.md` (GEN.104 to GEN.111, GEN.115, ADM.36, USR.3), the uploaded documents of section 10, `cli/orbits.py`, `db/store.py`, `physics/kepler.py`, `physics/state_vectors.py`, `physics/position.py`, `galaxy/geometry.py`.
- No web search or fetch results were available.
