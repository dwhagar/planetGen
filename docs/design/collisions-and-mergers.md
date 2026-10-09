# Collisions and mergers: how often, what happens, and how to detect it

How often rogue planets and small bodies actually meet, what the physics says
happens when they do, how the orbital update should find a collision without
tunneling or rounding errors, and what the admin report records. It turns
Boss's rogue collision rule (GEN.110) into a table the code can follow, carries
the still-valid material of "Orbital Update Collision Detection.md" forward,
and lists where that document is wrong. The orbital update that calls this code
is described in [orbital-updates.md](orbital-updates.md); asteroid field
classes and the rendering plan are in
[nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md);
the rogue and debris densities are in
[interstellar-object-rates.md](interstellar-object-rates.md).

Informs: GEN.110, GEN.109, GEN.105, GEN.112, GEN.47, GEN.75, GEN.111, ADM.36

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result (URL in Sources), [C] computed in
the research sandbox, [R] recalled and unconfirmed. The research environment
could only read search-result text, not the papers themselves; every [R] is on
the Evidence notes list. The [C] values for giant radii, Hill radii,
catastrophic speeds, the terrestrial pair rate, fusion shares, the Roche table
and the section 2.1 failures were rerun while writing this document.

## Decisions already taken

- Rogue collision rule (Boss, 2026-10-03 05:38Z): "If rogue planets ever hit
  and at least one is terrestrial both are destroyed and they turn into an
  asteroid field. If two rogue gas giant collide we assume they become a
  bigger gas giant losing only a small % of their mass, then we have to find
  out if they would have sufficient mass for nuclear fusion to occur and then
  add a star to the map if that's yes. If a star does get added to the map
  that gets reported to admin as well and the star should get no planets and
  instead should be a planetary nebula. Anything like that."
- Influence rule (Boss, 2026-10-09): the "nearest 10 influencers of equal or
  greater mass" scheme is replaced. An object's influence radius is the Hill
  sphere of the largest nearby object; every point mass inside it, plus the
  galactic gradient, perturbs the object at an update. Section 3 applies it.
- Time step (Boss, 2026-10-07 17:11Z): each run advances by the real time
  since the last one, default 1 day = 1 day, with an option for a stated extra
  span in one go.
- Sectors (Boss, 2026-10-07 17:11Z): keep the ring, layer and slot sectors,
  with forces summed in galactic coordinates, "but verify the algorithm will
  work with our sector geometry."

## Findings in brief

- Rogue collisions are astronomically rare: one rogue Earth meets another about
  once per 4e22 years, and across a real-sized galaxy the first rogue-rogue
  collision of any kind comes about once per 2e9 years [C]. The handler is a
  safety net; expect none in a game run unless an admin aims bodies at each
  other (ADM.36).
- As built, nothing can collide: `rogue_planets` and `asteroid_fields` have no
  velocity columns, only a circular-rotation speed
  (`galactic_orbital_speed_kms`), so neighbours co-rotate at near-zero relative
  speed. A peculiar-velocity draw (real rogues: 20 to 35 km/s per axis) is a
  prerequisite for GEN.109 and GEN.110.
- Boss's terrestrial rule is close to the truth: an Earth-Earth hit at galactic
  speed is mostly a hit-and-run (63 to 75%) or a supercatastrophic disruption
  (22 to 24%), and a merger is 2% or less [C]. But a terrestrial body hitting a
  giant, brown dwarf or star is swallowed with no field, and two giants (at
  most 26 Mjup) cannot make a star, which needs 78.6 Mjup and so a brown dwarf
  in the pair.
- "Planetary nebula" is the wrong object: it is the envelope of a dying 1 to 8
  Msun star lit by a 30,000 K or hotter white dwarf [S]; a newly lit 0.075 Msun
  star is a 2,500 to 2,800 K red dwarf that ionises nothing.

## 1. What "Orbital Update Collision Detection.md" gives and what stays valid

The document replaces end-of-step distance checks with continuous collision
detection (CCD), a test of the whole straight path of a step. The same material
is in "Orbital Update Full Algorithm.md" (`CheckContinuousSphereOverlap`,
`ExecuteInelasticMerger`, `CheckPlanetaryHillCrossings`) and "Computational
Astrodynamics.md" (the section "Large Macro-Step Integration, Temporal
Baselines, and Continuous Collision Sweeping" and the encounter micro-pass).

**Kept: the geometry** (source section "Mathematical Foundation of Continuous
Collision Detection"). With dp0 = p1(0) - p2(0) and dv = v1 - v2 over tau in
[0, dt]: S(tau) = A tau^2 + B tau + C with A = |dv|^2, B = 2 dp0.dv, C = |dp0|^2;
the time of closest approach is tau_min = -B / 2A (negative: separating at the
start; above dt: after the step); contact when the minimum separation is at most
R1 + R2, entering at tau_impact = (-B - sqrt(B^2 - 4A(C - (R1+R2)^2))) / 2A.
Section 5.2 computes the same quantities in a form that survives floating point.

**Kept: the inelastic merger.** M_new = M1 + M2 and v_new = (M1 v1 + M2 v2) /
M_new; right for the cases that merge (section 5.1). Additions: the merged body
sits at the centre of mass at impact, (M1 p1(t) + M2 p2(t)) / M_new, not at
body 1's position (orbital-updates.md section 10.6), and the momentum check goes
in the admin report (section 7).

**Kept: Roche (tidal) disruption.** A smaller body of density rho_2 within
d_roche = alpha R_1 (rho_1 / rho_2)^(1/3) of a primary of radius R_1 and density
rho_1 is torn apart. The form is right; alpha depends on the body (section 2.3).

**Kept: collisions are tested inside the update**, not in a separate O(N^2)
pass. The candidate set the document names is wrong (section 2.4).

## 2. Corrections to the source document

### 2.1 The swept-sphere code loses precision exactly when it matters

`check_continuous_collision` is correct on paper and wrong in float64. For a
body starting far away, B^2 and 4AC are both near 1e42 and their difference,
which carries the miss distance, is near 1e23, below float64 rounding at that
size (about 1e26). Compared against exact 120-digit arithmetic, 3,000 random
geometries per row, 30 km/s relative speed, two Earth-size bodies (R_sum =
1.27e7 m) [C]:

| Step length | Start distance | Source code wrong | Closest-approach form wrong |
|---|---|---|---|
| 1 day | 2.6e9 m | 0 | 0 |
| 1 year | 9.5e11 m | 0 | 0 |
| 100 years | 9.5e13 m | 8 | 0 |
| 1,000 years | 9.5e14 m (6,300 AU) | 578 | 0 |
| 10,000 years | 9.5e15 m (0.3 pc) | 1,491 (every real hit missed) | 0 |
| two 1 km bodies, 1 day | 2.6e9 m | 0 | 0 |
| two 1 km bodies, 1 year | 9.5e11 m | 1,468 | 0 |

The relative impact-time error is up to 2.5e-8 for the source form and 1e-15
for the other; an independent rerun against exact rational arithmetic gave the
same pattern. A one-day step (Boss's default) is safe for planets, but the
"extra span in one go" option and small bodies are not. "Orbital Update Full
Algorithm.md" (`CheckContinuousSphereOverlap`) has the same cancellation.

Other defects in the same code:

- `C <= radius_sq` returns a hit at t = 0 even when the pair is separating, so
  fragments just made by a collision, or a moon and its planet, would collide
  forever. The caller must exclude pairs created in the same run.
- Each pair is tested twice (once per body in the loop) and the partner is
  deleted mid-loop. Sort the pair's ids so one pair gives one event, resolve
  events in time order, and never delete mid-loop.
- `A <= 1e-18` has no unit (in SI, 1 nm/s of relative speed); use `A == 0.0`.
  The text promises "Adaptive Time-to-Impact Bisection" but never implements it,
  and the closed form needs none.
- The focusing radius of orbital-updates.md section 10.5, R_eff = (R1 + R2)
  sqrt(1 + v_esc^2 / max(dv^2, floor)), is a capture cross-section. As a
  straight-line test it is valid only along an asymptote (speed at infinity);
  the contact time and point it implies are not physical, and the speed for the
  outcome is v_imp = sqrt(v_inf^2 + v_esc^2), not dv. Cap R_eff at the smaller
  body's Hill radius (two Earths at 1 m/s would give 1 AU).

### 2.2 The merged radius is wrong for giants

`new_radius = (R1^3 + R2^3)^(1/3)` is volume addition. It is fine for two
terrestrial bodies of equal density and wrong for gas giants, whose radius is
nearly independent of mass. `physics.planets.giant_radius_km` with
`GIANT_JOVIAN_MASS_RADIUS = (317.8, 11.21, -0.04)` gives, in equatorial Jupiter
radii (71,492 km) [C]:

| Mass | 1 Mjup | 2 Mjup | 13 Mjup | 78 Mjup |
|---|---|---|---|---|
| Radius (Rjup) | 1.00 | 0.97 | 0.90 | 0.84 |

Volume addition for two Jupiters gives 1.26 Rjup (about 90,000 km) against 0.97
(69,500 km), 30% too large. Chen and Kipping (2017) find the same Jovian slope,
-0.04 +/- 0.02, with brown dwarfs on the trend to about 80 Mjup [S:
https://iopscience.iop.org/article/10.3847/1538-4357/834/1/17; weakly
constrained at the brown-dwarf end]. Use `giant_radius_km(M_new)` for giants
and brown dwarfs.

### 2.3 The Roche factor 2.44 is for fluid bodies only

The source documents use d_roche = 2.44 R_1 (rho_1/rho_2)^(1/3) for every pair,
and orbital-updates.md section 10.3 repeats it. The coefficient is 1.26 for a
rigid body and 2.44 to 2.45 for a fluid, self-gravitating one [S:
https://academic.oup.com/mnras/article/509/2/2404/6408484, citing Davidsson
1999]. Rubble piles split near the fluid value; monoliths survive some way
inside the rigid one. Use 1.26 (up to about 2 for rubble) for rocky and icy
bodies and 2.44 for fluid ones (gas giants, stars, molten worlds). It matters
only where the primary is much denser than the body:

| Primary | Body | Rigid limit (1.26) | Contact distance | Verdict |
|---|---|---|---|---|
| Earth-like rogue | Earth-like rogue | 1.26 R | 2 R | Contact comes first; no tidal regime |
| Jupiter (1.33 g/cc) | Earth | 0.79 Rjup | 1.09 Rjup | Inside the giant; swallowed |
| Sun | Earth | 0.80 Rsun | 1 Rsun | Inside the star |
| White dwarf (0.6 Msun) | Earth | about 4.7e8 m | 1.5e7 m | Disrupted 30 times beyond contact [C] |
| Neutron star | Earth | about 4.7e8 m | 6.4e6 m | Disrupted well before contact [C] |
| 10 Msun black hole | Earth | about 1e9 m | horizon 3e4 m | Disrupted outside the horizon [C] |

A separate Roche event is needed only for a compact-remnant primary (the
project has isolated neutron stars and black holes), and even then a terrestrial
rogue meets one about once per 1e19 years [C].

Two implementation faults in the source procedures. In
`CheckPlanetaryHillCrossings` the Roche test sits inside the branch that runs
only on the step a body enters the Hill sphere, where the distance is near the
Hill radius and far outside d_roche, so it would almost never fire; the
micro-pass version tests every sub-step, which is right, and the Hill-crossing
version should test every pair already inside on every update. And the
procedure only deletes the body and adds its mass ("stage ring formation" is a
comment); for a compact remnant record a `tidal_disruption` event (section 5.1)
and add the mass.

### 2.4 The candidate set is wrong, and top-10 is superseded

The document tests collisions only against "candidate gravitational
influencers" (the top 10). Boss's 2026-10-09 decision replaces top-10 with the
influence radius (section 3), but even that set is the wrong collision list. A
partner must be found by a mass-blind distance query of radius (R1 + R2) + |dv|
dt. The old "as large or larger" restriction misses a small body running into a
large one when the small body is the one updating, and the Hill-sphere set is
too small for long steps (a body at 56 km/s travels farther than an Earth's
0.021 pc Hill radius in about 370 years).

### 2.5 Other points

- Positions must be differenced from sector-relative values, not absolute
  galactic coordinates (section 5.3).
- Brent (1973) and Ortega and Rheinboldt (1970), cited for the promised
  bisection, are root-finding references the closed form does not need.
- orbital-updates.md sections 4, 5, 10.3 and 10.5 still carry the uncorrected
  forms (top 10 influencers, volume-sum radius, fluid Roche factor for all
  bodies, R_eff as a contact point). This document is the authority for
  collisions until they are updated.

## 3. Influence radius, Hill spheres and warnings

### 3.1 The rule (Boss, 2026-10-09)

For an object at an update, take the Hill sphere of the largest nearby object
(its own included). The perturbers are every point mass inside that radius plus
the galactic gradient (GEN.115).

- A rogue has no host star, so its Hill sphere is taken against the galaxy. For
  a flat rotation curve, r_H = (G m / (2 Omega^2))^(1/3), Omega = 220 km/s / 8.2
  kpc = 26.8 km/s/kpc. This matches the project's
  `physics.orbits.calculate_hill_sphere(a, m, M_enc) = a (m / 3 M_enc)^(1/3)`
  within 12% when M_enc is the galaxy mass inside the orbit (1.26 pc against
  1.44 pc for the Sun) [C]. Use one function everywhere. With the GEN.115 curve
  (229 km/s at 8.128 kpc) the radii are about 4% smaller.
- The search covers this sector and the sectors touching it, which reaches at
  least 4 pc. A larger radius (a 50 Msun star has about 5.3 pc) is clipped to
  the stencil; the gradient term covers the rest.

Radius and expected number of other bodies inside it at local density (stars
0.14, all rogues 0.91, brown dwarfs 0.03 per pc^3) [C]:

| Largest object | r_H | Bodies inside, average |
|---|---|---|
| Earth | 0.021 pc (4,300 AU) | 4e-5 |
| Jupiter | 0.14 pc | 0.013 |
| 40 Mjup brown dwarf | 0.49 pc | 0.5 |
| 0.3 Msun star | 0.96 pc | 4 |
| Sun | 1.44 pc | 13.5 (1.75 stars, 11.7 rogues) |
| 10 Msun | 3.1 pc | 135 |

In empty space the set is empty and only the gradient acts; near a star it is a
dozen or so; in a cluster or the galactic centre it is unbounded and needs a
cap. Recommended: the 64 nearest by distance, with the 1 pc Plummer softening
GEN.115 already specifies, otherwise the cost is O(N^2).

### 3.2 How big the effect is

A star at 1 pc pulls with 1.4e-13 m/s^2, changing a body's velocity by 1e-8 m/s
in a day, a thousandth of the galactic centripetal acceleration (1.9e-10
m/s^2). A deflection of at least 1 degree at 56 km/s needs an impact parameter
b = G M_tot / (v^2 tan(0.5 degree)): 2.9e4 km for Earth + Earth (about contact),
0.062 AU for Jupiter + Jupiter, 32 AU for any rogue + Sun. A rogue passes the
Sun that closely once per 1.6e12 years [C], so flyby deflection of unbound
rogues is never seen in game time; the influence rule matters for bound
hierarchies (a planet, its star, a passing star) and for ADM.36 edits.

### 3.3 A Hill-sphere warning is rare or constant, depending on its definition

Another rogue entering a rogue Earth's 4,300 AU sphere happens 6.1e-8 times per
object per year. Any rogue entering a star's 1.4 pc sphere happens about 3.4e-4
times per star per year, with about 12 rogues inside every star's sphere at any
instant [C]. A body at 56 km/s crossing a Hill sphere is unbound (the escape
speed at the edge is 1.1 m/s for Earth, 7.6 m/s for Jupiter, 77 m/s for the
Sun). Left as written, a full real-size galaxy would raise about 200 warnings a
day (6.5e11 rogues x 6.1e-8 x 3.9 / 2 / 365), and 0.3 a day for 1e9 rogues [C].

So GEN.109 and GEN.111 should not warn on "inside the Hill radius" alone. Warn
(and mail the admin) only when the pair can interact: v_rel^2 < 2 G M_tot / d
(bound or capture candidate), or the straight-line pericentre gives a
deflection of at least about 1 degree, or the pair is in the same system. Log
plain Hill-sphere entries at debug level.

## 4. What happens in a collision

### 4.1 Impact regimes (Leinhardt and Stewart 2012; Asphaug et al. 2006)

Leinhardt and Stewart (2012, ApJ 745:79) classify outcomes as cratering,
merging, disruption, supercatastrophic disruption and hit-and-run, set by impact
speed, mass ratio and angle [S: https://iopscience.iop.org/article/10.1088/0004-637X/745/1/79].
Asphaug, Agnor and Williams (2006, Nature 439:155) showed hit-and-run is common
[S: https://www.nature.com/articles/nature04311]; Genda, Kokubo and Ida (2012)
give a critical speed in units of v_esc above which a collision is hit-and-run
[S: https://arxiv.org/pdf/1109.4330].

- Mutual escape speed at contact v_esc = sqrt(2 G (M1 + M2) / (R1 + R2)); impact
  speed from an unbound approach with speed at infinity v_inf:
  **v_imp = sqrt(v_inf^2 + v_esc^2)** [C].
- Specific impact energy **Q_R = 0.5 mu v_imp^2 / M_tot**, mu = M1 M2 / (M1 + M2).
- Catastrophic disruption threshold (largest remnant is half the total mass):
  **Q*_RD = [(gamma+1)^2 / (4 gamma)]^(2/(3 mu_bar) - 1) x Q*_RD,gamma=1**, with
  gamma = M_projectile / M_target and Q*_RD,gamma=1 = c* (4/5) pi rho_1 G R_C1^2,
  rho_1 = 1000 kg/m^3, R_C1 the radius of a body of the total mass at that
  density. Strengthless planets: c* = 1.9, mu_bar = 0.36; small bodies: c* = 5,
  mu_bar = 0.37 [S: https://arxiv.org/html/1106.6084 as quoted in search
  results; https://arxiv.org/pdf/2006.01881 for the formula].
- Largest remnant M_lr = M_tot (1 - Q_R / (2 Q*_RD)) below Q*_RD [R, derived
  from the definition: linear through (0, 1) and (Q*, 1/2)]; beyond about 1.8
  Q*_RD, M_lr = 0.1 M_tot (Q_R / (1.8 Q*_RD))^-1.5 [R, about 70% sure].
- Grazing: b_crit = R_target / (R_target + R_projectile), b the sine of the
  impact angle. Above it only a lens of the projectile overlaps, with
  interacting fraction (3 R_p l^2 - l^3) / (4 R_p^3), l = (R_t + R_p)(1 - b)
  [R]. For isotropic encounters P(b < x) = x^2, so for equal sizes (b_crit =
  0.5) **75% of hits are grazing**.

Head-on catastrophic speeds [C]:

| Pair | Q*_RD (J/kg) | v_cat (km/s) | v_esc (km/s) | v_inf needed | v_cat / v_esc |
|---|---|---|---|---|---|
| Earth + Earth | 6.4e7 | 22.6 | 11.2 | 19.7 | 2.0 |
| Earth + 0.1 Earth | 1.1e8 | 51.7 | 9.5 | 50.8 | 5.5 |
| Earth + 0.01 Earth | 6.4e8 | 362 | 10.0 | 362 | 36 |
| Moon-mass pair | 3.4e6 | 5.2 | 2.2 | 4.7 | 2.3 |
| Jupiter + Jupiter | 3.0e9 | 155 | 60 | 142 | 2.6 |

Protoplanets in a disc meet at a few km/s, below v_esc, so mergers dominate
there. Rogues meet at the galactic dispersion (mean v_inf about 56 km/s for 25
km/s per axis in each body), 2 to 5 times the catastrophic speed for comparable
terrestrial bodies. A small body hitting a big one is not disruptive even at 100
km/s; it erodes and merges. The coefficients for gas giants extrapolate scaling
laws calibrated on rocky and icy bodies.

Monte Carlo of isotropic encounters (Maxwellian v_inf, P(b) ~ b db, hit-and-run
when b > b_crit and v_imp > 1.5 v_esc [a design value, R], 35 km/s relative 1-D
dispersion, 20,000 trials per row) [C]:

| Pair | hit-and-run | merge | partial erosion | disruption | supercatastrophic | mean debris, fraction of M_tot |
|---|---|---|---|---|---|---|
| Earth + Earth | 74% | 2% | 0% | 2% | 22% | 0.55 |
| Earth + 0.1 Earth | 57% | 2% | 18% | 11% | 12% | 0.28 |
| Earth + 0.01 Earth | 37% | 28% | 36% | 0% | 0% | 0.01 |
| Jupiter + Earth (giant target) | 5% | 95% | 0% | 0% | 0% | 0.00 |
| Jupiter + Jupiter | 22% | 78% | 0% | 0% | 0% | 0.06 |

At a colder 22 km/s mean v_inf the Earth-Earth row is 63% hit-and-run, 22%
merge, 15% disruption or worse. These are design-grade shares, not predictions.

### 4.2 Gas giant mergers

- **Mass loss is small.** SPH simulations find that quick Jupiter-Jupiter mergers
  (pericentre below 2 Rjup) keep at least 97% of the mass, give fast-spinning
  puffy remnants, and can be treated as perfect inelastic collisions [S:
  https://academic.oup.com/mnras/article/501/2/1621/6027698, search snippets
  only]; a gas giant hit by an Earth-mass embryo kept about 90% of its envelope
  [S: https://arxiv.org/pdf/1410.6815 snippet]. The scaling-law formula gives 4
  to 12% mean loss across the galactic speed range, in line with Boss's "small %".
- **The remnant spins near breakup.** A grazing merger (b = 0.5) of two Jupiters
  at v_inf = 0 to 30 km/s spins at 1.0 to 1.1 times breakup (moment of inertia
  0.25 M R^2); at b = 0.8 and 50 km/s about 2 times [C]. Flag a merged giant
  "rapidly rotating, inflated, hot"; it may shed a disc (rings or moons, the
  existing `has_moons` field).
- **Radius** comes from the relation (section 2.2).

### 4.3 Fusion thresholds

| Threshold | Mjup | Msun | Source |
|---|---|---|---|
| Deuterium burning (50% burn) | 13.0 +/- 0.8 (range 11 to 16.3) | 0.0124 | [S] Spiegel, Burrows and Milsom 2011, ApJ 727:57 |
| Lithium burning | about 63 to 68 | 0.060 to 0.065 | [S] https://arxiv.org/abs/2110.11982 |
| Hydrogen burning, solar metallicity | 73 to 79 (75 cloudless; 78.5 with ATMO models) | 0.070 to 0.075 | [S] Chabrier and Baraffe 2000; Burrows 2001; https://arxiv.org/pdf/1312.1736 |
| Hydrogen burning, zero metallicity | about 96 | 0.092 | [S] Burrows 2001 as quoted |

1 Mjup = 9.54e-4 Msun [C], so 13 Mjup = 0.0124, 78.6 = 0.075, 80 = 0.0763 Msun.
The project already uses 13 to 80 Mjup for rogue brown dwarfs
(`ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER`) and 0.08 Msun as the lowest IMF break
(`IMF_BREAKS_SOL`, which also floors `star.py`'s initial mass). Recommended
"lights as a star" line: **0.075 Msun (78.6 Mjup)**, as a setting, flagged
"borderline" in the report within 3% of it. A merger star of 0.075 to 0.08 Msun
is below the generator's 0.08 floor, so either lower the floor for merger stars
or clamp them to 0.08.

A merger below the line becomes an ordinary **brown dwarf** (deuterium-burning
briefly above 13 Mjup): reclassify the row (mass_bin 'brown-dwarf', recompute
radius and temperature), do not make a star. Chance a merger crosses each line,
partners drawn from the generator's own classes (Jupiter bin: power law of slope
0.65 over 1 to 13 Mjup; brown dwarfs log-uniform over 13 to 80), no mass loss
[C]:

| Pair | above 13 Mjup | above 65 | above 75.4 | above 78.6 | above 83.8 |
|---|---|---|---|---|---|
| Jupiter-bin + Jupiter-bin | 8.2% | 0 | 0 | 0 | 0 |
| Jupiter-bin + brown dwarf | 100% | 14% | 5.8% | 3.5% | 0.6% |
| Brown dwarf + brown dwarf | 100% | 59% | 45% | 41% | 35% |
| Brown dwarf + brown dwarf, 5% mass loss | 100% | 54% | 40% | 36% | 29% |

Boss's "check if fusion occurs, add a star" therefore triggers only when a
rogue brown dwarf is involved; write the check for any two substellar bodies.

### 4.4 What a newly lit star looks like

A 0.075 to 0.08 Msun merger remnant is a very late M dwarf (T_eff about 2,500 to
2,800 K, L about 1e-4 to 1e-3 Lsun, radius about 0.1 Rsun [R]) and makes
essentially no ionising radiation. A transient is plausible but unobserved for
planets: the analogues are luminous red novae (V1309 Sco, about 1e45 erg [S])
and the planet-engulfment transient ZTF SLRN-2020 (about 6.5e41 erg [S:
https://www.nature.com/articles/s41586-023-05842-x]), against about 6e46 erg of
impact energy for a 40 + 40 Mjup merger [C], mostly heating the remnant. Honest
summary: a brief red infrared flare lasting weeks to months at most, or nothing
visible. The ejecta (a few percent of 80 Mjup at 100 to 400 km/s) would be a
cold, dusty, dark shell, about 1 cm^-3 at 1 ly after 1,000 years and below the
interstellar medium a few thousand years later [C].

## 5. Rewrite of Boss's rule, outcome table and detection

### 5.1 The rule and the outcome table

Kinds: T (rocky or icy rogue below `ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER`,
0.05 Mjup, about 16 Earth masses), G (gas giant, 0.05 to 13 Mjup, Saturn and
Jupiter bins), B (brown dwarf, 13 to 80 Mjup), S (star), X (neutron star or
black hole). Compute v_inf from the stored velocities, v_esc at contact, v_imp =
sqrt(v_inf^2 + v_esc^2), b from the swept geometry (impact parameter over R1 +
R2), and x = Q_R / Q*_RD (section 4.1, interacting mass for grazing hits).

| Kind 1 | Kind 2 | Condition | Outcome | Debris | Record |
|---|---|---|---|---|---|
| T | T | b < b_crit and x < 0.1 (or v_imp < 1.2 v_esc) | Merge into one T of summed mass; radius by volume | none | merger |
| T | T | b >= b_crit and v_imp <= 1.5 v_esc | Graze-and-merge, as above | none | merger |
| T | T | b >= b_crit, v_imp > 1.5 v_esc, x < 1 | Hit-and-run: both survive on slightly changed paths, each loses its interacting lens | field of the eroded mass if at least 1% of M_tot, else only logged | hit_and_run |
| T | T | 0.1 <= x < 1, b < b_crit | Partial erosion: M_lr = M_tot (1 - x/2) | field of M_tot - M_lr | erosion |
| T | T | x >= 1 | **Both destroyed**: one field of mass M_tot, with M_lr (0.5 down to about 0.1 M_tot) as a hero body | all of M_tot | disruption |
| T | G, B, S | any | T absorbed; the larger body gains the mass; for S a transient may be logged | none | absorbed |
| T | X | tidal limit (section 2.3) | Tidally disrupted before contact, accreted | none stored | tidal_disruption |
| G or B | G or B | any | Merge, lose f = clamp(Q_R / (2 Q*_RD), 0.01, 0.5); M_new = (1 - f) M_tot; radius from relation; fast-spin flag; gas lost, no field | none | giant_merger |
| (merged) | | M_new < 13 Mjup | stays a gas giant (re-pick the bin) | | |
| (merged) | | 13 <= M_new < 78.6 Mjup | becomes a brown dwarf (mass_bin 'brown-dwarf') | | brown_dwarf_formed |
| (merged) | | M_new >= 78.6 Mjup | **new star**: no planets, age 0, "merger remnant" flag, admin report | | star_created |
| S or X | S or X | any | outside Boss's scope; default: merge conserving momentum, log, admin report | | stellar_merger |

If Boss prefers the literal rule, replace the T + T rows with "any contact: both
destroyed, one field of M_tot" (setting `ALWAYS_DESTROY_TERRESTRIAL`).

Options for the nebula clause:

| Option | What | Consequence |
|---|---|---|
| A (recommended) | No nebula. Star flagged "merger remnant"; event log and admin report carry the story; optional flavour line "a faint infrared flare was recorded" | Physically honest; no new object type; no conflict with `PLANETARY_NEBULA_CENTRAL_STAR_TYPES` (O3VII to B0VII, 30,000 K and up) |
| B | A small, dark, short-lived "merger ejecta shell" in the dark nebula family (radius 0.1 to a few ly, removed after about 1e4 years) | Needs a seeded `NebulaShape` (GEN.75) per event and a removal job; adds a nebula type no class table describes |
| C | Keep the name "planetary nebula" as a label | Wrong; classes H to L in the nebula classes document require a hot white dwarf; not recommended |

**The asteroid field an event produces.**

- Mass: M_tot minus any surviving bodies. `asteroid_fields` has no mass column
  today; add `mass_kg`, `formed_at`, `origin_event_id` and a velocity (the
  centre-of-mass velocity at impact; the update needs one to move the field).
- Position: the centre of mass at impact, (M1 p1(t_hit) + M2 p2(t_hit)) / M_tot,
  advanced to now.
- Composition: stony by default; carbonaceous or icy if the parents were
  water-rich (`ROGUE_WATER_RICH_CHANCE`); mixed if they differ. Class letter U,
  collisional family.
- Radius: starts at about R1 + R2 (1e4 km) and grows as r0 + sigma_disp t with
  sigma_disp 1 to 3 km/s. It reaches the smallest field, 0.001 ly (63 AU), after
  100 to 300 years and 1 ly after 1e5 to 3e5 years [C]. The size digit is
  clamped to 1 to 4 (`ASTEROID_FIELD_SIZE_DIGIT_RANGE`,
  `ASTEROID_FIELD_RADIUS_RANGE_LY` from 0.001 ly), so store the real radius
  (display "<63 AU") and recompute `radius_ly` from `formed_at` at each update.
- Members (rendering only, GEN.112): cumulative N(>D) ~ D^-2.5 with explicit hero
  bodies (differential slope 3.5 aged, 2.2 to 2.7 for a fresh disruption [S:
  https://arxiv.org/pdf/0911.3937 snippet]); the largest fragment is M_lr.
- Life span: free fields disperse in 1e6 to 1e7 years (Raymond et al. 2020, in
  interstellar-object-rates.md). The generator makes no free fields
  (`PHENOMENON_DENSITY_PC3["asteroid-field"]` is 0), so event-made and hand-made
  fields are the only ones. Let galactic shear stretch the radius, which the
  update does anyway, or retire them in the maintenance job.

**Recording an event.** A new table `collision_events` holds one row per event
and the products' ids. Destroyed bodies are not deleted outright: set
`destroyed_by_event_id` and hide them (or move them to a history table) so the
admin can inspect the inputs and replay the outcome from the stored seed.
Outcome draws use
`draw.Stream(f"{galaxy_seed}:collision:{uid_a}:{uid_b}:{t_hit_unix}")`, the same
pattern as `phenomenon_scatter`.

### 5.2 The detection routine (closest-approach form)

Positions are SI, differenced from sector-relative values; R = R1 + R2; the
window is [0, dt].

```python
def sweep(dp, dv, R, dt):
    """dp, dv: relative position and velocity at the window start
    (SI, differenced from sector-relative values, never from absolute
    galactic coordinates); R = R1 + R2; window [0, dt].
    Returns (hit, t_hit)."""
    d2 = dp[0]**2 + dp[1]**2 + dp[2]**2
    R2 = R * R
    if d2 <= R2:
        return True, 0.0
    A = dv[0]**2 + dv[1]**2 + dv[2]**2
    if A == 0.0:
        return False, dt
    t_min = -(dp[0]*dv[0] + dp[1]*dv[1] + dp[2]*dv[2]) / A
    if t_min < 0.0:                     # separating
        return False, dt
    mx, my, mz = (dp[0] + dv[0]*t_min, dp[1] + dv[1]*t_min, dp[2] + dv[2]*t_min)
    m2 = mx*mx + my*my + mz*mz          # miss distance squared, vector form
    if m2 > R2:
        return False, dt
    t_hit = t_min - math.sqrt((R2 - m2) / A)
    return (t_hit <= dt), (max(t_hit, 0.0) if t_hit <= dt else dt)
```

It finds the time of closest approach and the miss distance as a vector first,
then the entry time, so nothing large is subtracted from something nearly equal.
Tested [C]: head-on hit and time; arrival after the window; hit at the window
end; separating outside; already overlapping at t0 (a hit at 0, so the caller
excludes pairs created in the same run); tangent graze b = R (hit); miss by 1e-7
relative (no hit); zero relative speed outside; zero dt; moving partner. NaN
input is rejected before the call. All pass, and 20,000 random small-scale cases
against dense time sampling (4,000 samples per window) show 0 mismatches in
hit/miss and in impact time. Large-scale accuracy is from the exact-arithmetic
comparison in section 2.1.

### 5.3 Pair windows, curved paths, precision

- **Windows.** Bodies update asynchronously (`next_update_due`), each storing
  (position, velocity, epoch) valid until its next update. When X updates over
  [e_X_prev, e_X_now], test it against every other body Y's stored segment
  extended to e_X_now, over **[max(e_X_prev, e_Y), e_X_now]**. Y's earlier path
  was covered when Y updated, so the union has no gap and each pair is tested
  once per update of either body. Resolve events in time order, drop pairs that
  include a body already consumed in the same run, and never delete mid-loop.
- **Relative frame.** The galactic acceleration (1.9e-10 m/s^2) is common to
  both bodies and cancels; only the tidal part matters. Against an RK4
  orbit in a flat-rotation-curve potential (two Earths, 30 km/s, miss distance
  5,000 km) the straight relative path is off by 14 to 26 km after 1,000 years
  and 22,000 km after 10,000 years (more than R_sum = 12,700 km) [C]. The fit is
  about 0.03 Omega^2 v_rel dt^3, so cap the window at
  tau_max = (0.01 R_sum / (0.03 Omega^2 v_rel))^(1/3): about 1,800 years for two Earths, 100 years for two 1 km bodies. Cut longer spans
  into windows, each re-propagated by the Verlet step. The one-day default is
  far inside all of these.
- **Sector reach.** A pair must not move more than one sector (4 pc) in a
  window: 7e4 years at 56 km/s, 3,900 years at 1,000 km/s. Cap windows by the
  smaller of that and tau_max.
- **Precision.** 8 kpc is 2.5e20 m, where float64 spacing is 3.3e4 m (28 km in
  parsec doubles). `rogue_planets.center_x/y/z_pc` are absolute doubles, so
  rogue positions are good to about 14 km at 8 kpc and 28 km at 16 kpc: fine for
  planets (R_sum 1.3e4 km), useless for asteroids. Difference from
  sector-relative offsets (the star systems table has local `position_*_mpc`
  doubles; rogues would need the same) or from exact integers on a milliparsec
  grid (1 mpc = 206 AU).

### 5.4 Broad phase

Reach in a window is |dv_max| dt + R_sum: about 2.6e6 km (8.4e-8 pc) for one
day at 30 km/s, a sphere holding about 2e-21 other rogues at local density
[C]. The design is a cheap index, not pair loops: per sector and its touching
neighbours (the ring, layer and slot address maths planned for GEN.109), a
sorted array along one axis or a k-d tree over the due bodies' segment bounding
boxes inflated by R_sum. Skip a pair when |dp| - |dv| dt > R_sum. A sector holds
about 58 rogues and 9 stars at local density, so a sorted sweep is O(n log n)
with n below 100, and the narrow phase almost never runs. The same neighbour
search supplies the perturbing point masses under the influence rule, but
collision candidates use their own mass-blind query (section 2.4).

## 6. Rates

For one body of species 1 moving through species 2 of density n (Maxwellian
relative speeds, 1-D relative dispersion s_rel = sqrt(s1^2 + s2^2)):

**Gamma = n pi b^2 [ <v> + (2 G M_tot / b) <1/v> ]**, with b = R1 + R2, <v> = sqrt(8/pi) s_rel, <1/v> = sqrt(2/pi) / s_rel.

The second term is gravitational focusing, sigma = pi b^2 (1 + v_esc^2 / v^2).
Check against the published focusing-dominated stellar collision time, T_coll
about 7e14 yr (R/Rsun)^-1 (M/Msun)^-1 (n/pc^-3)^-1 (sigma/km s^-1) [S:
https://arxiv.org/html/1302.2549]: 1.2e17 yr for n = 0.14 and sigma = 25 km/s;
this formula gives 1.26e17 yr [C]. Inputs are the project's densities
(`tuning.PHENOMENON_DENSITY_PC3`, `tuning.ROGUE_PLANET_MASS_BINS`): terrestrial
0.78, sub-Neptune 0.10, Saturn-class 0.024, Jupiter-bin 0.0039, brown dwarf 0.03,
neutron star 7e-4, star 0.14 per pc^3, and 1e12 pc^-3 for comets and debris
(`INTERSTELLAR_DEBRIS_DENSITY_PC3`). The thin disc has (sigma_U, sigma_V,
sigma_W) about (33, 28, 23) km/s overall [S: https://arxiv.org/pdf/1710.08479;
https://www.aanda.org/articles/aa/full_html/2024/08/aa49445-24/aa49445-24.html,
search snippets]; rogues are ejected at a few km/s [R], so they inherit it.
Default 25 km/s per axis per body (mean relative speed 56 km/s).

| Pair | n_target (pc^-3) | v_esc,mutual (km/s) | Rate per object (per yr) | Mean wait (yr) | Galaxy-wide, 1e11 stars, per day |
|---|---|---|---|---|---|
| terrestrial rogue + terrestrial rogue | 0.78 | 11.2 | 2.5e-23 | 4e22 | 7.6e-14 |
| terrestrial rogue + giant rogue | 0.028 | 58 | 7.3e-23 | 1.4e22 | 4.4e-13 |
| giant + giant | 0.028 | 60 | 3.6e-23 | 2.8e22 | 2.8e-14 |
| brown dwarf (40 Mjup) + brown dwarf | 0.03 | 380 | 6.6e-21 | 1.5e20 | 7.4e-13 |
| Jupiter-bin + brown dwarf | 0.03 | 272 | 3.4e-21 | 2.9e20 | 1.0e-13 |
| rogue (Earth-like) + star | 0.14 | 615 | 2.0e-18 | 5e17 | 1.4e-8 (all 6.5 rogues per star) |
| terrestrial rogue + neutron star | 7e-4 | 7,600 | 1.3e-22 | 8e21 | 7.5e-13 |
| star + star (reference) | 0.14 | 618 | 7.9e-18 | 1.3e17 | n/a |
| terrestrial rogue + 1 km comet (harmless hit) | 1e12 | 11 | 8e-12 | 1.2e11 | 0.05 |

The galaxy-wide column is N1 x per-object rate x 3.9 (the density weighting of
an exponential disc, exp(R0/L)/6 with L = 2.6 kpc, R0 = 8.2 kpc, no bulge), with
a factor 1/2 for identical species [C]. The five rogue-rogue lines sum to about
1.4e-12 per day, one event per 2e9 years for a real-size galaxy; rates scale as
N^2, so a generated galaxy of 1e9 rogues is smaller by a further 3e-6. Comet
hits are craters, not events to model, which is why comets must not be rows.
Dense places do not change the picture: a terrestrial rogue pair waits about
1.3e17 years in a globular cluster core and 1e14 to 1e15 years at the galactic
centre [S/R: https://arxiv.org/html/1302.2549;
https://iopscience.iop.org/article/10.3847/2041-8213/ad251f; the central density
is recalled]. Belt and moon collisions (main belt: 2.85 +/- 0.66 e-18 km^-2
yr^-1 at 5.81 km/s [S: Farinella and Davis 1992]) are the orbit code's concern,
not this handler's.

## 7. What the admin report contains

One row per event in `collision_events` (also emailed per GEN.111 when SMTP is
configured and the outcome is in the notify set), plus a count line in the
orbital-update summary (GEN.107). It records:

- the event: id, run id, impact time t_hit, time detected, outcome code (the
  Record column of section 5.1), regime name and seed string;
- each participant: id, name, kind, class, mass, radius, composition summary,
  galactic velocity, sector and galactic coordinates at impact, and whether it
  was admin-edited (ADM.36);
- the numbers: separation at start and impact, b and b_crit, v_inf, v_esc,
  v_imp, Q_R, Q*_RD, x, interacting-mass fraction, the pair's catastrophic
  velocity, the window used and any cap, and whether a focusing radius was used;
- a conservation check: mass and momentum before and after to a stated
  tolerance, impact kinetic energy and the share assumed deposited as heat;
- the products: ids, masses, radii; for a field, class, mass, radius and growth
  law; for a giant merger, mass lost, new radius and spin-to-breakup ratio; for
  a fusion check, the merged mass, the thresholds used and "borderline" within
  3%. Star creation always goes to the admin;
- context and flags: the containing nebula or remnant (`NebulaShape.contains`),
  nearest stars, the number of influence-radius bodies, sectors marked dirty,
  suppressed duplicates, a pair already overlapping at window start, and
  `numeric_warning` when the position resolution (14 to 28 km) exceeds 1% of
  R_sum.

## Evidence notes

Items marked [R] or partly confirmed; check when arXiv, ADS and Wikipedia are
reachable (the research environment read search-result text only):

- Leinhardt and Stewart: the linear remnant branch is derived from the
  definition of Q*; the second branch, the lens formula and b_crit are recalled
  (70 to 80% sure of the constants). The 1.5 v_esc hit-and-run threshold is a
  design value; Genda et al.'s formula was not retrieved.
- The Monte Carlo shares depend on a simplified model, an Earth mass-radius law
  R ~ M^0.27 and isotropic encounters; giant behaviour extrapolates rocky
  scaling. They are order-of-magnitude.
- Merger-star properties (T_eff, luminosity, radius), ejecta speeds, the galactic
  centre density (1e6 to 1e7 stars per pc^3 inside 0.1 pc), open cluster inputs,
  thick disc (60, 40, 35) and bulge (about 110 km/s) dispersions are recalled.
- The galaxy-wide factor 3.9 ignores the bulge and thick disc.
- The Jacobi radius (G m / 2 Omega^2)^(1/3) is the standard tidal radius in a
  logarithmic potential; no source was found.
- The giant radii in section 2.2 are in equatorial Jupiter radii (71,492 km);
  the research draft's 1.02/0.97/0.92/0.86 used a smaller Jupiter radius.

## Sources

Seen in search results (search text only):

- Leinhardt and Stewart 2012, ApJ 745:79: https://iopscience.iop.org/article/10.1088/0004-637X/745/1/79 ; https://arxiv.org/html/1106.6084 ; https://arxiv.org/pdf/2006.01881 (formula restated). Stewart and Leinhardt 2012: https://arxiv.org/pdf/1109.4588 (title only).
- Asphaug, Agnor and Williams 2006: https://www.nature.com/articles/nature04311 ; Genda, Kokubo and Ida 2012: https://arxiv.org/pdf/1109.4330
- Spiegel, Burrows and Milsom 2011: https://www.osti.gov/etdeweb/biblio/21567589 ; hydrogen limit: https://arxiv.org/pdf/1312.1736 ; https://www.researchgate.net/publication/367543291 ; https://arxiv.org/pdf/astro-ph/9902015 ; lithium limit: https://arxiv.org/abs/2110.11982
- Giant mergers (SPH): https://academic.oup.com/mnras/article/501/2/1621/6027698 ; https://arxiv.org/pdf/1410.6815
- Chen and Kipping 2017: https://iopscience.iop.org/article/10.3847/1538-4357/834/1/17 ; https://arxiv.org/pdf/2311.12593
- ZTF SLRN-2020: https://www.nature.com/articles/s41586-023-05842-x ; V1309 Sco: https://arxiv.org/pdf/1311.6522
- Planetary nebulae: https://arxiv.org/pdf/1411.2365 ; https://arxiv.org/pdf/1002.1525
- Roche coefficients: https://academic.oup.com/mnras/article/509/2/2404/6408484 ; https://arxiv.org/html/2110.07601
- Main belt: https://sciencedirect.com/science/article/abs/pii/001910359290060K ; https://arxiv.org/pdf/2403.03248
- Velocity dispersions: https://arxiv.org/pdf/1710.08479 ; https://www.aanda.org/articles/aa/full_html/2024/08/aa49445-24/aa49445-24.html ; https://arxiv.org/html/2412.07089v1
- Stellar encounter rates: https://arxiv.org/html/1302.2549 ; nuclear star cluster: https://iopscience.iop.org/article/10.3847/2041-8213/ad251f ; https://arxiv.org/html/0810.0204
- Free-floating planets: https://arxiv.org/html/2503.11597 ; https://iopscience.iop.org/article/10.3847/1538-3881/ace688/pdf
- Fragment size distribution: https://arxiv.org/pdf/0911.3937

Repo files read: "Orbital Update Collision Detection.md" (in full), "Orbital
Update Full Algorithm.md", "Computational Astrodynamics.md", orbital-updates.md,
interstellar-object-rates.md, nebula-and-asteroid-field-classes.md,
`docs/TODO.md` (GEN.105 to GEN.112, ADM.36), `src/planetgen/tuning.py`,
`physics/planets.py`, `physics/orbits.py`, `physics/constants.py`,
`galaxy/galactic_orbit.py`, `generation/phenomena/rogue.py`,
`generation/phenomena/asteroid_field.py`, `db/models.py`.
