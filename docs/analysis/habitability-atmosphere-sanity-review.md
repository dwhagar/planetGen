# Habitability / Atmosphere Sanity Review

**Date:** 2026-09-06
**Scope:** Manual, domain-informed sanity review of `planetPhysics.calculate_surface_gravity`
and `planetPhysics.calculate_atmospheric_conditions` output, specifically to answer the
question TODO.md raises: now that the atmospheric-pressure formula has been fixed
(barometric-formula rewrite, `planetPhysics.py` ~L391-405), are the disabled Class M/P
forced clamps (gravity ~L342-347, pressure/temperature ~L407-419) actually safe to leave
removed?

**Method:** Built `Planet` objects directly (bypassing `StarSystem`'s orbit-placement
loop, the same lightweight pattern `tests/test_planets.py` uses) against one fixed G2V
host star (mass ≈0.89 M☉, luminosity ≈0.67 L☉, habitable zone ≈0.78–1.13 AU — a
plausible, if not dead-center-solar, Sun-analog). For each of the 23 `PLANET_CLASSES`
entries with at least one valid zone, generated 300 planets (moons disabled) and recorded
`gravity`, `surface_temperature`, `atmospheric_pressure`, `density`, `atm_density`,
`atm_molar_density`, `scale_height`, `mass`, `radius`. Terrestrial ecosphere-only classes
used a fixed distance at the habitable-zone midpoint for the main batch (isolates
non-distance sources of variance); a second, supplementary run for Class M and P instead
drew a random distance uniformly across the whole zone width per sample, to approximate
what real system generation would actually produce. Class R was skipped — it declares
`h/e/c` all `False` in `PLANET_CLASSES`, so it has no valid zone under normal generation
(it can only be force-assigned via `zone_override`, out of scope here).

Scripts used are not part of the repo (scratch analysis only); the underlying formulas
referenced below are in `stellarObjects/planetPhysics.py` and
`stellarObjects/physical_constants.py`.

## Class definitions confirmed (program_constants.PLANET_CLASSES)

| Class | Radius range (km) | Type | Valid zone(s) | Atmosphere |
|---|---|---|---|---|
| M | 5000–10000 | t | e only | "a mix of oxygen, nitrogen, and argon" |
| P | 5000–10000 | t | e only | "a mix of oxygen, nitrogen, and argon (thinning with age)" |

M and P have **identical radius ranges, identical density ranges** (both terrestrial,
`PLANET_DENSITY["t"]` = 3.93–5.51 g/cm³), and **identical atmosphere-parameter ranges**
(`ATMOSPHERE_DENSITY["t"]` = 0.02–1.2 kg/m³, `ATMOSPHERIC_MOLAR_DENSITY["t"]` =
0.02897–0.04347 kg/mol) — nothing about their generation differs except their text
description and life-chemical list. Every other terrestrial atmosphere-bearing class
draws from the same two `"t"` ranges too (Class N is the sole special case, hard-coded to
`atm_density = 65`, `atm_molar_density = max`).

## Headline finding: the removed Class P clamp was carrying all of P's "cold world" identity, and it's gone

Both the fixed-distance and randomized-distance runs show **Class P is not colder than
Class M**:

| Run | M mean temp (K) | P mean temp (K) | M stdev | P stdev |
|---|---|---|---|---|
| Fixed distance (zone midpoint) | 276.35 | 279.12 (P *warmer*) | 20.64 | 20.23 |
| Randomized distance (full HZ) | 281.77 | 276.75 (P colder, but by <1 stdev) | 24.87 | 25.10 |

The sign of the M-vs-P difference **flips between the two runs**, and in both cases the
gap (2.8 K and 5.0 K respectively) is small relative to the ~20-25 K standard deviation —
i.e. statistically indistinguishable from noise, not a systematic "P runs colder" effect.
This is expected once you trace the formula: P has no distance bias, no albedo bias, no
atmosphere-parameter bias relative to M or any other terrestrial ecosphere class. The
*only* thing that ever made P reliably cold was the disabled clamp
(`elif planet.planet_class == "P" and planet.surface_temperature >= 283: ...`, ~L416-422).
With it gone, **"a cold, glaciated world" is currently indistinguishable from "a
terrestrial Earth-like world" in every generated physical quantity.** This is squarely
the class-differentiation failure mode the task asked me to check for.

**Verdict: not safe to leave removed for Class P.** Recommend restoring some form of
cold-bias for P — either the original clamp, or (better long-term) giving P its own
higher-albedo range or an outer-HZ-biased distance preference so the *unclamped* physics
naturally trends cold instead of being force-overridden after the fact.

## Class M: pressure never actually reaches Earth-like values

| Field | Fixed-distance run | Randomized-distance run | Real Earth |
|---|---|---|---|
| `gravity` (g) | min 0.59, max 1.54, mean 1.00, sd 0.23 | min 0.57, max 1.54, mean 1.02, sd 0.22 | 1.0 |
| `surface_temperature` (K) | min 233.4, max 318.8, mean 276.4, sd 20.6 | min 218.5, max 341.3, mean 281.8, sd 24.9 | ~288 |
| `atmospheric_pressure` (Pa) | min 1358, max 80947, mean 33681, sd 19398 | min 1572, max 90507, mean 32404, sd 20101 | 101325 |

- **Gravity** clusters right on 1.00g on average, which is good — but the *spread*
  (0.57-1.54g) is a direct, mechanical consequence of the class's own `radius_range`
  (5000-10000 km) crossed with `PLANET_DENSITY["t"]` (3.93-5.51 g/cm³): gravity is
  proportional to `density × radius`, and the class's own declared ranges alone produce
  a ~2.7x spread between the smallest/least-dense and largest/densest instance. This is a
  property of the class's declared radius/density envelope, not of the atmospheric-pressure
  fix under review — restoring the gravity clamp is a separate judgment call about whether
  "Earth-like" should mean "tightly Earth-mass" or "plausible super-Earth/sub-Earth
  terrestrial", independent of the pressure/temperature question.
- **Pressure** is the real problem: across 600 total samples (both runs combined),
  **zero** reached even 0.9 atm (91,193 Pa), let alone Earth's actual 101,325 Pa. The
  mean sits at roughly **1/3 atm**. This directly contradicts the TODO's premise that
  "the corrected atmospheric-pressure formula ... lands close to realistic ranges on its
  own" — for Class M specifically, it does not; it systematically undershoots.
  - Root cause: `atmospheric_pressure = atm_density * R * T / atm_molar_density` (the
    `scale_height`/gravity terms cancel algebraically, see below), and reproducing
    Earth's real 101,325 Pa needs `atm_density` near/above the sampled range's *ceiling*
    (1.2 kg/m³, itself only just above Earth's real 1.225 kg/m³) **simultaneously** with
    `atm_molar_density` near the range's *floor* (0.029 kg/mol, Earth's real air molar
    mass) — i.e. Earth sits in a corner of the two independently-uniformly-sampled
    parameters' joint space, not their center. Most draws land well short of that corner,
    hence the ~1/3-atm mean.
- **Temperature**: mean 276-282K is a plausible few degrees below Earth's 288K, but the
  spread is wide enough that a non-trivial share of "Earth-like" worlds come out well
  below freezing (36% of the randomized-distance M sample was below 0°C, with one draw as
  low as 218K / -55°C). A single named "Earth-like" class regularly generating worlds
  colder than Antarctica's average is a sanity failure regardless of the clamp-removal
  question.

**Verdict: not safe to leave removed for Class M either**, specifically for pressure (never
reaches realistic Earth values) and, more marginally, temperature (long cold tail). Gravity
alone looks fine on its own (mean ≈1g) but inherits a wide spread from the class's radius/
density definition that's an independent question from the pressure formula fix.

## Likely root cause behind both M's and P's issues: the greenhouse-factor formula looks physically inverted

```python
greenhouse_factor = abs((planet.atm_molar_density - physical_constants.CO2_BASE_MOLAR_DENSITY)
                         / physical_constants.CO2_BASE_MOLAR_DENSITY
                         * program_constants.CO2_MAX_GREENHOUSE_FACTOR)
```

This computes "how far is this atmosphere's average molar mass from pure CO2's own molar
mass" and uses *that distance* as the greenhouse boost — meaning an atmosphere whose molar
density happens to sit **near CO2's own value gets ~zero greenhouse boost**, while one that's
**far from CO2 (in either direction)** gets a large boost. That's backwards from reality
(CO2-heavy atmospheres are strong greenhouse absorbers; drifting away from CO2 isn't what
produces warming). Two concrete symptoms observed in the data:

- **Class N** ("a hot world with a dense, reducing atmosphere") is hard-coded to
  `atm_density = 65` (correctly Venus-like) and `atm_molar_density = max` — which happens to
  sit almost exactly at `CO2_BASE_MOLAR_DENSITY` (0.04345 vs 0.04347). Its greenhouse_factor
  is therefore ≈0, and its measured surface temperature (mean 240K, range 231.7-249.9K) comes
  out **colder than Class M** despite N carrying a ~30-atm CO2-dominated atmosphere (mean
  pressure 2.98 MPa) that in reality would produce an extreme, Venus-like greenhouse effect. A
  class explicitly designed to be "hot" with a "dense" atmosphere is the coldest terrestrial
  class sampled.
- **Gas giants (I, J, S, T, U)** draw `atm_molar_density` from `ATMOSPHERIC_MOLAR_DENSITY["g"]`
  (0.00226-0.00416 kg/mol), which is always *very* far from CO2's molar mass, so their
  greenhouse_factor sits pinned near its ceiling (`CO2_MAX_GREENHOUSE_FACTOR = 5`) regardless
  of the actual draw. H2/He atmospheres are not meaningfully greenhouse-active in reality, so
  applying a near-maximum CO2-style greenhouse boost to every gas giant is another instance of
  the same inversion.
- Because `atm_molar_density` is drawn independently at random for every class (it is *not*
  derived from the class's own described atmosphere composition text at all, except for the
  N special case), this formula effectively injects large, physically ungrounded temperature
  noise into every terrestrial and gas-giant class alike — which is very likely the dominant
  reason the disabled M/P clamps existed in the first place, and why removing them produces
  wide, undifferentiated temperature distributions across nearly every class.

This is a substantive greenhouse-model design question (what should stand in for "how much
CO2/greenhouse gas is actually present," since the codebase doesn't currently track that as
its own quantity) rather than a one-line fix, so **I did not change it** — flagging for the
user's judgment rather than a clamp-restoration decision.

## Other terrestrial/atmosphere-bearing classes (type "t", atmosphere ≠ None)

Sampled A, B, E, F, G, H, K, L, N, O, Q, V, W, Y (X, C, D have `atmosphere: None` and were
sampled only for gravity/temperature/mass/radius, not atmospheric fields).

| Class | Zone | Gravity mean (sd) | Temp mean K (sd) | Pressure mean Pa (sd) |
|---|---|---|---|---|
| A | h | 0.38 (0.17) | 433.8 (30.7) | 51,350 (31,541) |
| B | h | 0.38 (0.18) | 434.9 (31.7) | 55,974 (32,349) |
| E | e | 1.00 (0.21) | 278.2 (20.7) | 34,418 (19,907) |
| F | e | 1.02 (0.21) | 278.7 (20.2) | 34,203 (19,809) |
| G | e | 1.01 (0.22) | 278.2 (20.7) | 34,418 (19,971) |
| H | e | 1.00 (0.21) | 278.7 (19.6) | 34,359 (19,672) |
| K | e | 0.68 (0.20) | 278.0 (20.1) | 33,933 (19,355) |
| L | e | 0.84 (0.13) | 278.9 (19.7) | 34,588 (19,774) |
| M | e | 1.00 (0.23) | 276.4 (20.6) | 33,681 (19,398) |
| N | e | 1.02 (0.22) | 240.1 (5.3) | 2,983,190 (66,105) |
| O | e | 1.02 (0.22) | 280.1 (19.9) | 33,360 (19,466) |
| P | e | 1.01 (0.21) | 279.1 (20.2) | 34,879 (20,058) |
| Q | e | 0.66 (0.23) | 278.7 (19.4) | 34,155 (19,969) |
| V | e | 1.67 (0.25) | 279.1 (19.3) | 33,775 (19,743) |
| W | e | 0.73 (0.38) | 276.3 (20.5) | 34,351 (19,445) |
| Y | h | 0.84 (0.13) | 435.5 (30.4) | 52,121 (30,888) |

**Finding: nearly every ecosphere-zone terrestrial class (E, F, G, H, K, L, M, O, P, Q, V, W)
produces a statistically indistinguishable atmospheric-pressure and surface-temperature
distribution.** Compare, e.g., K ("an adaptable world with a *thin* atmosphere") to M
("Earth-like"): K's pressure (mean 33,933 Pa) is not meaningfully thinner than M's (mean
33,681 Pa) — if anything they're identical within noise. Q ("variable atmosphere, thin to
dense") and V ("Super-Earth", the only class with a visibly distinct gravity mean, ≈1.67g
from its larger radius range 10000-15000km) both land on the same pressure/temperature
distribution as everyone else. This is because `atm_density`/`atm_molar_density` are drawn
from the same `"t"`-keyed range for every terrestrial class regardless of the individual
class's flavor text — the class descriptions ("thin", "dense", "variable", "toxic") are
currently cosmetic and not reflected in the actual generated physics for anything but the
one N special case. This is a pre-existing characteristic of the generation model (not
something the pressure-formula fix changed), but it's directly relevant to "do classes
differentiate from each other when they should," so it's included here.

A, B, and Y (hot-zone-only classes) show a visibly different, hotter distribution — but
that's purely because they're placed at a much closer/hotter distance (zone 'h'), not
because of any atmosphere-parameter difference from the ecosphere classes.

## Gas giants (type "g": I, J, S, T, U)

| Class | Radius range (km) | Real analog | Gravity mean (range) | Density mean g/cm³ (range) | Pressure mean Pa (range) |
|---|---|---|---|---|---|
| I | 15,000–50,000 | Neptune/Uranus | 0.33 (0.017–1.20) | 0.36 (0.026–0.91) | 428,199 (228,039–754,212) |
| J | 25,000–250,000 | Jupiter-and-beyond | 1.42 (0.06–5.99) | 0.36 (0.035–0.93) | 433,488 (222,441–736,523) |
| S | 250,000–50,000,000 | "supergiant" (no real analog) | 262.5 (1.97–1267.6) | 0.37 (0.026–0.95) | 428,982 (215,351–765,181) |
| T | 250,000–25,000,000 | "gas dwarf" (no real analog) | 134.9 (1.35–553.6) | 0.37 (0.027–0.97) | 417,907 (232,454–717,192) |
| U | 25,000,000–60,000,000 | "ultragiant" (no real analog) | 413.0 (25.5–1355.4) | 0.34 (0.027–0.96) | 418,861 (218,257–743,364) |

Real-world reference: Jupiter surface (1-bar level) gravity ≈2.53g, bulk density ≈1.33
g/cm³, reference pressure by definition ≈1 atm (101,325 Pa); Saturn gravity ≈1.06g, bulk
density ≈0.69 g/cm³; Neptune gravity ≈1.14g, bulk density ≈1.64 g/cm³.

Two clear findings here, one of them a demonstrable, mechanical bug in how gas-giant
density is blended, and one a mathematically provable formula property:

### 1. Gas-giant density blending (planetPhysics.py L308-311) produces implausibly "fluffy" planets

```python
if planet.body_type == 'g':
    core_to_atmosphere_ratio = random.uniform(*program_constants.GAS_GIANT_CORE_ATMOSPHERE_RATIO)
    planet.density = planet.density * core_to_atmosphere_ratio + (1 - core_to_atmosphere_ratio) * (planet.atm_density / 1000)
```

`GAS_GIANT_CORE_ATMOSPHERE_RATIO` (0.03-0.6) is documented as "the average ratio of a gas
giant's core mass to its **total mass**" — a *mass* fraction. This code uses it directly as
a weight in a simple arithmetic blend of two *densities* (bulk-rock-like core density ~0.69-
1.64 g/cm³, and atmosphere density rescaled to ~0.0007-0.0013 g/cm³). Blending densities by
mass-fraction like this is not physically meaningful (density is intensive; combining two
components' densities correctly requires their volume fractions, roughly a harmonic-mean-like
combination, not a mass-fraction-weighted arithmetic mean of densities). Because the
atmosphere's density contribution is ~1000x smaller than the core's, and the core's mass
fraction is frequently the low end of its range (down to 3%), the blend frequently produces
overall planet densities as low as **0.026-0.035 g/cm³** (min values seen across all five gas
giant classes) — meaningfully fluffier than any known giant planet or even known "puffy"
hot-Jupiter exoplanets, and it's this artificially low density that in turn produces the
very low gravity values seen for I/J (as low as 0.017g and 0.064g respectively) despite
still having Neptune-to-Jupiter-scale radii. This looks like a genuine modeling bug, but
fixing it correctly requires deciding what `GAS_GIANT_CORE_ATMOSPHERE_RATIO` should really
mean (mass fraction vs. volume fraction) and re-deriving the blend accordingly — a design
decision, not a one-line change, so left for the user's judgment rather than fixed here.

### 2. Atmospheric pressure is mathematically independent of surface gravity — for every class, not just gas giants

Tracing the formula:

```
scale_height_m        = R * T / (atm_molar_density * g_ms2)
atmospheric_pressure   = atm_density * g_ms2 * scale_height_m
                        = atm_density * g_ms2 * [R * T / (atm_molar_density * g_ms2)]
                        = atm_density * R * T / atm_molar_density        # g_ms2 cancels exactly
```

Verified numerically (holding `atm_density`, `atm_molar_density`, `T` fixed and varying
gravity from 0.2g to 50g): the resulting pressure is **bit-identical** at every gravity
value. A supplementary correlation check across 300 randomly-generated Class M and P
planets each gives `corr(gravity, pressure) ≈ 0.02` and `-0.14` respectively — statistical
noise, not a relationship. This explains why classes I, J, S, T, and U — spanning a
~1,200,000km radius range and a >1000x gravity range (0.06g to 1355g) between them — all
produce essentially the **same** atmospheric-pressure distribution (all cluster around
400,000-430,000 Pa / ~4 atm, regardless of class). Real atmospheric retention is strongly
gravity-dependent (it's a large part of why Mars, at 0.38g, has ~0.6% of Earth's pressure
despite comparable early volatile inventories); this model's `atm_density` (the free,
independently-drawn parameter that stands in for "how much atmosphere was retained") has
no coupling to gravity/escape velocity at all, so this decoupling is real and applies to
every class, terrestrial and gas giant alike — it isn't a Class-M/P-specific issue, and the
module's own docstring ("estimates the atmospheric pressure based on the atmospheric mass
and planet's gravity") no longer matches what the code actually does once `scale_height` is
substituted in. This is a modeling-completeness question (should `atm_density` be derived
from a gravity/escape-velocity-aware retention model rather than sampled independently) —
worth the user's attention, but re-deriving `atm_density` from first principles is a design
decision, so it's flagged here rather than changed.

## Fix applied

`stellarObjects/physical_constants.py`'s `ATMOSPHERE_DENSITY["t"]` comment previously read
"Terrestrial: Range from Mars to Venus" for the range `(0.02, 1.2)` kg/m³. Real Venus surface
air density is ~65 kg/m³ — roughly 50x this range's own ceiling — while real Earth sea-level
air density (~1.225 kg/m³) matches the range's ceiling almost exactly. The comment was
unambiguously wrong (not a judgment call) and made the terrestrial atmosphere-density range
look far more permissive than it actually is, which is directly relevant to why Class M's
pressure undershoots real Earth values (see above) — so it was corrected in place; no
behavior changed, this is a comment-only fix.

## Summary of recommendations

1. **Do not leave the Class P clamp removed as-is.** Without it, P is statistically
   indistinguishable from M (and from every other ecosphere terrestrial class) in
   temperature — its entire "cold, glaciated" identity currently depends on that disabled
   clamp. Either restore it, or give P its own parameters (higher albedo range, an
   outer-HZ-biased placement, or a class-specific greenhouse penalty) so the physics
   produces cold worlds without a post-hoc override.
2. **Do not leave the Class M pressure clamp removed as-is.** Across 600 samples spanning
   two different distance-sampling strategies, zero reached even 0.9 atm; the mean sits at
   ~1/3 atm. The "fixed formula lands close to realistic ranges on its own" premise does
   not hold for M's pressure. Either restore the pressure clamp, or narrow/re-center
   `ATMOSPHERE_DENSITY["t"]`/`ATMOSPHERIC_MOLAR_DENSITY["t"]` specifically for M so its
   own distribution centers nearer Earth's actual values.
3. **Class M's gravity clamp is a separable, lower-priority question.** Its mean already
   centers on 1.00g; the wide spread (0.57-1.54g) comes from the class's own declared
   radius/density envelope, not from anything the pressure-formula fix touched. Whether
   that spread is "too wide for Earth-like" is a judgment call about the class's radius/
   density range design, independent of the atmosphere-formula question.
4. **Investigate the greenhouse_factor formula** (`planetPhysics.py` ~L402) — it currently
   rewards atmospheres for being *far* from CO2's own molar density rather than for
   *containing* more CO2/greenhouse gas, which produces backwards results for at least
   Class N (a "dense," CO2-heavy world coming out colder than "Earth-like" M) and
   arguably for every gas giant (near-max greenhouse boost applied uniformly to
   non-greenhouse-active H2/He atmospheres). This is very likely the dominant source of
   the wide, undifferentiated temperature spread seen across nearly all classes, and is
   probably the real reason the M/P clamps existed in the first place.
5. **Most ecosphere terrestrial classes (E, F, G, H, K, L, M, N-partially, O, P, Q, V, W)
   are statistically indistinguishable from each other in pressure and temperature** since
   they all draw from the same `"t"`-keyed atmosphere-parameter ranges. Whether that's
   acceptable (they're meant to differ mainly in flavor text/life chemistry, not raw
   physics) or should be tightened per-class is a design call for the user; flagged here
   as a "classes don't differentiate" observation per the task brief.
6. **Gas-giant density blending (I/J/S/T/U) produces implausibly low densities** (as low
   as 0.026 g/cm³) from mixing a mass-fraction ratio into a simple arithmetic density
   blend; this is very likely a genuine bug, but the correct fix depends on clarifying
   what `GAS_GIANT_CORE_ATMOSPHERE_RATIO` is meant to represent, so it's left to the user.
7. **Atmospheric pressure is mathematically decoupled from gravity for every class** (the
   gravity term cancels exactly in the current formula) — real atmospheric retention is
   strongly gravity-dependent, and this is worth a deliberate design decision about
   whether/how to reintroduce that coupling (e.g. deriving `atm_density` from an
   escape-velocity-aware retention model) rather than sampling it independently.

No `planetPhysics.py` clamp logic was touched, per the task's constraints — all six
substantive findings above are flagged for the user's decision. Only the one
unambiguous, behavior-preserving documentation fix (item under "Fix applied") was made.
