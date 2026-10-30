# The planet mass law without a class (GEN.650)

Status: built as a pure law (`physics/planet_mass.py`, tests in `tests/test_planet_mass.py`) and wired into generation by the planet architecture (GEN.728, wired on merge train 4): `planet_architecture` builds a whole chain of these cores (`oligarch_core_kg`) and thins it, and a planet asked for alone (HABITABLE_WORLD, FORCE_CLASS, SLOTS) draws one. Written by the condition-based classes lane with Habitability lane I (the bulk contract), reworked by the Planets 1 lane. Times EDT.

## Decided

- Boss, 6:39 PM EDT 2026-10-23, pipeline card 1, option 1: the planet mass is the isolation mass at the formation distance and the disc's solids, with a scatter and the giant core-accretion rule, tuned against the observed occurrence by the GEN.503 harness (builders GEN.509 and GEN.510). The isolation law is checked above a few Earth masses. The class mass ranges of bulk-composition.md decision 28 are not kept as a hidden stand-in beyond the spine's class-first epoch.

## The law (`planet_mass.planet_mass(distance_au, star_mass_kg, snow_line_au, fe_h, stream)`)

1. **Solids.** The isolation mass at the formation distance, `formation.isolation_mass_kg`, for a feeding zone `OLIGARCH_SEPARATION_MEDIAN` (15.07) mutual Hill radii wide (drawn per system, `OLIGARCH_SEPARATION_SIGMA_LN`), from the minimum-mass nebula scaled by the star's mass and 10^[Fe/H] (`disk_surface_density_scale`), the disc's own scatter (`DISC_SOLIDS_SCATTER_DEX`, 0.7 dex, Andrews et al. 2013) and boosted by `condensate_shares.WATER_ICE_SOLID_BOOST` (2.17, the water ice of Lodders 2003 Table 11 over the rock; GEN.950 replaced Hayashi's 30/7) beyond the star's own snow line (Lissauer 1993; Kokubo & Ida 2000), times `PLANET_DISC_SOLIDS_ENHANCEMENT` (0.4647, stand-in PA1: the solids against the minimum-mass nebula, tuned to the cold-giant occurrence), and inside 1 AU times `(distance / 1 AU) ^ -PLANET_INNER_DISC_STEEPENING` (0.75, stand-in PM2). Beyond `PLANET_OUTER_GROWTH_AU` (11.36 AU) the slow outer disc has built less: the mass falls as the distance to the power `PLANET_OUTER_GROWTH_EXPONENT` (2.889, stand-in PA1). The isolation formula is built and unchanged.
2. **Assembly.** Giant impacts merge a feeding zone's oligarchs; that is the architecture's thinning (`planet_architecture.sculpt`, Pu and Wu 2015), not a constant survival share, so the old `GIANT_IMPACT_SURVIVAL_FRACTION` (0.4) and `PLANET_ASSEMBLY_SCATTER_DEX` are gone. A planet asked for alone is one oligarch of its feeding zone.
3. **Runaway.** A core of `GIANT_CORE_CRITICAL_MASS_EARTH` = 10 Earth masses or more, beyond the snow line, becomes a gas giant with probability `GIANT_GAS_SURVIVAL_SHARE` (0.766, stand-in PA2) times `10^(GIANT_GAS_SURVIVAL_METALLICITY_EXPONENT [Fe/H])` (0.6498, PA2) times the occurrence factor for the star's [Fe/H] (`metallicity.giant_occurrence_factor`, the existing 10^(2 [Fe/H]) law, capped at 1), times `(GIANT_GAS_SURVIVAL_CUT_AU / distance)^2` beyond `GIANT_GAS_SURVIVAL_CUT_AU` (9.022 AU, PA2) (the outer cores grow too slowly to reach the critical mass before the gas is gone). The giant's mass follows dN/dlogM ~ M^-0.31 (Cumming et al. 2008, the slope already in `GIANT_MASS_FUNCTION_SLOPE`) between `GIANT_RUNAWAY_MIN_MASS_EARTH` = 30 and the lighter of 13 Jupiter masses and 1% of the star. A core that finds no gas stays a solid or ice-giant core, at most `GIANT_NO_GAS_CORE_MAX_EARTH` = 20 Earth masses (the envelope share of such a body is the composition lane's draw, GEN.428).
4. **Cap.** No planet's solids exceed `PLANET_MAX_STAR_MASS_RATIO` (1%, the existing `GIANT_MAX_STAR_MASS_RATIO`) of its star.

Draw order on the body's stream when a planet is asked for alone: the disc scatter and the spacing (`draw_disc_scatter`, `draw_separation`), then (only for a core past the critical mass beyond the line) the gas-survival draw, then the runaway mass. The same stream gives the same planet.

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

## Tuned against the occurrence of planets (2026-10-28, GEN.711, Aux 16's check)

Aux 16's distribution check (`reports/planet-class-distribution-check-2026-10-25.md`) found the class draw far off the observed occurrence: about 50 times too many hot Jupiters around M dwarfs, about 10 times too many giants around M dwarfs, a missing 2 to 20 Earth-mass population. The first version of this law (enhancement 5, survival 0.5, no cut) fixed the hot Jupiters but gave 1.2 giants of 100 Earth masses or more per Sun-like star and far too few super-Earths around M dwarfs. Measured on 400 generated systems' real slot distances (seed 5; 317 M dwarfs, 28 Sun-like stars), planets per star:

| Per star | First version | Tuned | Observed (recalled, rounded) |
|---|---|---|---|
| Giants (100 Earth masses or more) 0.3 to 10 AU, Sun-like | 0.51 | 0.16 | 0.10 to 0.15 (Fulton et al. 2021) |
| Giants, any orbit, M dwarf | 0.19 | 0.04 | 0.03 to 0.06 (Johnson et al. 2010; Bonfils et al. 2013) |
| 1 to 20 Earth masses inside 0.45 AU, Sun-like | 0.52 | 0.88 | 0.5 to 0.9 (Petigura et al. 2018; Zhu et al. 2018) |
| 1 to 20 Earth masses inside 0.45 AU, M dwarf | 0.24 | 1.30 | 1 to 2.5 (Dressing & Charbonneau 2015) |
| Earth-size (0.5 to 2 Earth masses) at 0.7 to 1.5 AU, Sun-like | 0.15 | 0.15 | 0.2 to 0.6 (Bryson et al. 2021) |

The rates are quoted from memory and need a check against the papers before they are cited elsewhere. The eta-Earth row stays at the low edge of the observed band; the checks of the tuned values are the tests `test_planet_mass.py::test_a_sun_like_star_has_about_a_tenth_of_a_giant_between_a_third_and_ten_au` and the sweep script `src/tests/planet_mass_sweep_cli.py`.

## Stand-ins (labelled)

| Constant | Value | Source | State |
|---|---|---|---|
| `PLANET_DISC_SOLIDS_ENHANCEMENT` | 5 | Chiang & Laughlin 2013 (several times the MMSN), recalled | stand-in, tuned with the steepening (registry PM1) |
| `PLANET_INNER_DISC_STEEPENING` | 0.75 | Chiang & Laughlin 2013 (the extrasolar nebula falls as about the -1.6 to -2 power), recalled; tuned | stand-in (registry PM2) |
| `PLANET_ASSEMBLY_SCATTER_DEX` | 0.8 | planetGen default | stand-in (registry PM3) |
| `GIANT_GAS_SURVIVAL_SHARE` | 0.12 | tuned to Fulton et al. 2021's giant occurrence, recalled | stand-in (registry PM6) |
| `GIANT_GAS_SURVIVAL_CUT_AU` | 5 | tuned to the rarity of giants beyond 10 AU in direct imaging | stand-in (registry PM7) |
| `GIANT_NO_GAS_CORE_MAX_EARTH` | 20 | Uranus and Neptune are 14.5 and 17.1 Earth masses | stand-in (registry PM8) |
| `GIANT_CORE_CRITICAL_MASS_EARTH` | 10 | Pollack et al. 1996; Rafikov 2006 give 5 to 20 | recalled, not looked up this session |
| `GIANT_RUNAWAY_MIN_MASS_EARTH` | 30 | planetGen default (class J's floor today) | stand-in |
| assembly factor | 1 / 0.4 | mass conservation against `GIANT_IMPACT_SURVIVAL_FRACTION` (Chambers 2001, recalled) | derived |

Nothing here is fixed from the Solar System: every value is a distribution parameter, and the Solar System is not an input.
