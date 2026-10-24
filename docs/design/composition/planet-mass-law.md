# The planet mass law without a class (GEN.650)

Status: built as a pure law (`physics/planet_mass.py`, tests in `tests/test_planet_mass.py`), not wired into generation. The class cut-over (GEN.646) calls it. Written by the condition-based classes lane with Habitability lane I (the bulk contract). Times EDT.

## Decided

- Boss, 6:39 PM EDT 2026-10-23, pipeline card 1, option 1: the planet mass is the isolation mass at the formation distance and the disc's solids, with a scatter and the giant core-accretion rule, tuned against the observed occurrence by the GEN.503 harness (builders GEN.509 and GEN.510). The isolation law is checked above a few Earth masses. The class mass ranges of bulk-composition.md decision 28 are not kept as a hidden stand-in beyond the spine's class-first epoch.

## The law (`planet_mass.planet_mass(distance_au, star_mass_kg, snow_line_au, fe_h, stream)`)

1. **Solids.** The isolation mass at the formation distance, `formation.isolation_mass_kg`, from the minimum-mass nebula scaled by the star's mass and 10^[Fe/H] (`disk_surface_density_scale`) and boosted by 30/7 beyond the star's own snow line (Lissauer 1993; Kokubo & Ida 2000; Hayashi 1981), times `PLANET_DISC_SOLIDS_ENHANCEMENT` (5, stand-in: Kepler's inner planets need several times the minimum-mass nebula, Chiang & Laughlin 2013, recalled). The isolation formula is built and unchanged; at 1 times the nebula the inner planets come out at lunar masses.
2. **Assembly.** Giant impacts merge a feeding zone's oligarchs, and `GIANT_IMPACT_SURVIVAL_FRACTION` (0.4) of them survive as planets. By mass conservation a planet is on average 1 / 0.4 = 2.5 isolation masses, so the assembly factor is not a new knob: it is the same fraction that already sets how many planets a star gets (`estimate_num_objects`), which keeps count and mass consistent. The spread is log-normal, `PLANET_ASSEMBLY_SCATTER_DEX` = 0.8 dex (stand-in).
3. **Runaway.** A core of `GIANT_CORE_CRITICAL_MASS_EARTH` = 10 Earth masses or more, beyond the snow line, becomes a gas giant with probability `GIANT_GAS_SURVIVAL_SHARE` (0.5, stand-in) times the occurrence factor for the star's [Fe/H] (`metallicity.giant_occurrence_factor`, the existing 10^(2 [Fe/H]) law, capped at 1). The giant's mass follows dN/dlogM ~ M^-0.31 (Cumming et al. 2008, the slope already in `GIANT_MASS_FUNCTION_SLOPE`) between `GIANT_RUNAWAY_MIN_MASS_EARTH` = 30 and the lighter of 13 Jupiter masses and 1% of the star. A core that finds no gas stays a solid or ice-giant core of its assembled mass (the envelope share of such a body is the composition lane's draw, GEN.428).
4. **Cap.** No planet's solids exceed `PLANET_MAX_STAR_MASS_RATIO` (1%, the existing `GIANT_MAX_STAR_MASS_RATIO`) of its star.

Draw order on the body's stream: the assembly scatter, then (only for a core past the critical mass beyond the line) the gas-survival draw, then the runaway mass. The same stream gives the same planet.

## Not in this law

- **Moons** are separate (GEN.466, parent and the circumplanetary ice line).
- **The formation distance.** The law takes it as an argument. Today a planet's slot distance is its formation distance. The migrant draw (a giant formed beyond the line now sitting close in, the "born beyond the snow line or migrant" logistic of GEN.510) supplies a different formation distance when it lands; until then there are no hot giants, which is the main difference from today's mix (below).
- **Radius, density, envelope.** Derived from the mass and composition by the spine (GEN.428, GEN.429), not drawn here.
- **Belts and small bodies** keep their own laws (GEN.442).

## What changes against today's class-first masses

Measured with 400 generated systems (3,205 planets, seed 11) by the sweep in the lane's report `reports/planet-mass-law-sweep-2026-10-24.md`: today the mass is a class draw with almost no link to the distance. The law is physical, so it differs where the class draw was not:

- Giants (30 Earth masses or more): 13.7% of planets today, 8.4% under the law. Inside the snow line today 5 to 14% of planets are giants; under the law under 1% (no migrants yet). Beyond the line 14 to 22% today, 9 to 21% under the law, rising with distance.
- Median mass: about 0.05 to 0.13 Earth masses at every distance today; the law gives 0.004 to 0.02 inside the snow line and 0.2 to 1.3 beyond, rising with distance (dim red-dwarf discs are light), and a wider spread (p10 to p90 spans four to five orders of magnitude against three today).
- This is explained, not hidden: the law is the intended physics and the numbers are stand-ins. The three knobs that move the giant share and the median are `PLANET_DISC_SOLIDS_ENHANCEMENT`, `PLANET_ASSEMBLY_SCATTER_DEX` and `GIANT_GAS_SURVIVAL_SHARE`; the migrant draw (GEN.510) restores close-in giants. Tuning them against the observed occurrence is the harness's job (GEN.503, GEN.509, GEN.510); the stand-ins were chosen by hand so that the giant share is the same order as today's (8.4% against 13.7%) and the median rises with distance.

## Stand-ins (labelled)

| Constant | Value | Source | State |
|---|---|---|---|
| `PLANET_DISC_SOLIDS_ENHANCEMENT` | 5 | Chiang & Laughlin 2013 (several times the MMSN), recalled | stand-in, tune with GEN.509 |
| `PLANET_ASSEMBLY_SCATTER_DEX` | 0.8 | planetGen default | stand-in, tune with GEN.509 |
| `GIANT_GAS_SURVIVAL_SHARE` | 0.5 | planetGen default, chosen against the sweep's giant share | stand-in, tune with GEN.510 |
| `GIANT_CORE_CRITICAL_MASS_EARTH` | 10 | Pollack et al. 1996; Rafikov 2006 give 5 to 20 | recalled, not looked up this session |
| `GIANT_RUNAWAY_MIN_MASS_EARTH` | 30 | planetGen default (class J's floor today) | stand-in |
| assembly factor | 1 / 0.4 | mass conservation against `GIANT_IMPACT_SURVIVAL_FRACTION` (Chambers 2001, recalled) | derived |

Nothing here is fixed from the Solar System: every value is a distribution parameter, and the Solar System is not an input.
