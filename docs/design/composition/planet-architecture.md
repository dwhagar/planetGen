# The planet architecture (GEN.728, GEN.729, GEN.730, GEN.828, GEN.829, GEN.920, TEST.395)

Status: wired into generation (`StarSystem._plan_architecture`, merge train 4, after GEN.950 and the star epoch GEN.932). The pure modules are `physics/planet_architecture.py`, `physics/planet_occurrence.py`, `physics/kepler_multiplicity.py` and `physics/data/planet_occurrence.json`. The wiring replaced `estimate_num_objects`, the uniform spacing ratio (`PLANET_SPACING_RATIO_RANGE`), the first-slot habitable-zone fraction and the constant survival (`GIANT_IMPACT_SURVIVAL_FRACTION` 0.4). Written by the Planets 1 lane. Times EDT.

## Decided

- Boss, 8:16 PM EDT 2026-10-29: build all of 8.1; the physics found in research wins; secondary sources are accepted by default; stand-ins are fine where they are indexed and registered (GEN.725, section 6.36 of the gap-filling plan).
- Lock-free first: the observed bands, the pure model, the harness and the tests merge on their own. The wiring changes generated output and waits for its train. At the wiring the architecture's own copy of the disc and gas numbers (the `ARCH_*` constants) became the mass law's (`PLANET_DISC_SOLIDS_ENHANCEMENT`, `PLANET_OUTER_GROWTH_*`, `GIANT_GAS_SURVIVAL_*`), and the sweep (L3) was run again on GEN.950's 182 K snow line and ice boost of 2.17.

## What it replaces

| Today | Fault found | Replaced by |
|---|---|---|
| A uniform planet count (`estimate_num_objects`), about 5 per Sun-like star and 8 per M dwarf | 60 percent of Sun-like stars had a planet of 1 Earth mass inside 1 AU; the surveys find 30 percent | A Bernoulli draw for a Kepler-like inner system at the observed share, then a chain |
| A uniform spacing ratio 1.3 to 1.8 (GEN.729) | Spacing in mutual Hill radii depends on the masses; the observed mode is near 20 with 93 percent of pairs over 10 | Oligarch spacing in mutual Hill radii, thinned by instability |
| A constant survival 0.4 of the oligarchs (GEN.730) | No time or spacing dependence | Pu and Wu 2015 survival law with event-driven collisions |
| No hot Jupiters; 0.02 giants per Sun-like star at 2 to 8 AU (GEN.828, GEN.920) | Observed 0.6 to 1.2 percent hot; 14 percent cold | Giant gas outcome with metallicity, inward migration, hot-Jupiter migration |
| No 2 to 20 Earth-mass population (GEN.829) | Missing | The inner chain's masses |

GEN.722 (radius valley) is already built on main (`radius_valley.py`, rows RV1 to RV4); it only needs closing.

## Sources (all in `physics/data/planet_occurrence.json`, each band READ or 2ND-HAND)

Zhu et al. 2018 (arXiv:1802.09526): 30 percent of Sun-like stars have Kepler-like planets, mean about 3, the transit multiplicity {1:432, 2:99, 3:42, 4:11, 5:3, 6:2}, the inclination law. Yang et al. 2020 (2002.02840): about 75 percent of stars under 5000 K. Zhu and Dong 2021 (2103.02127, review). Pu and Wu 2015 (1502.05449): survival law, eq. 11. Weiss et al. 2018 (1706.06204): peas in a pod, spacing mode near 20. Mulders et al. 2015 (1406.7356): inner edge, mass ratio. Johnson et al. 2010 (1005.3084): giants rise with mass and with 10^(1.2 [Fe/H]). Fulton et al. 2021 (2105.11584) and Cumming et al. 2008 (0803.3357): giant occurrence peaks at 1 to 3 AU. Petigura et al. 2018 (1712.04042): hot Jupiters go as 10^(3.4 [Fe/H]). Gan et al. 2023 (2210.08313): 0.27 percent of early M dwarfs. Wright et al. 2012 (1205.2273), Bryson et al. 2021 (2010.14812), Dressing and Charbonneau 2015 (1501.01623), Andrews et al. 2013 (1305.5262), Pascucci et al. 2016 (1608.03621).

## The model (`planet_architecture.build`)

1. **Disc scatter.** The disc's solids scatter 0.7 dex (Andrews 2013).
2. **Inner system.** Whether the star has one: a Bernoulli draw at the share anchors (0.85 at 0.15 Msun, 0.75 at 0.6, 0.30 at 1.0 and 1.4, log-linear in mass; Zhu 2018, Yang 2020). It starts with `1 + Poisson` oligarchs (the mean itself log-normal), neighbours alike in mass (peas in a pod), the first at a period near 5 days (Mulders 2015), each next at a spacing in mutual Hill radii.
3. **Outer chain.** From 1.3 times 400 days outward, at the isolation mass of the feeding zone with a steep taper beyond `PLANET_OUTER_GROWTH_AU`.
4. **Sculpting.** Pu and Wu's law, `x = 0.5 K - 0.6 - 0.35 log10 tau`, survival `1 - exp(-x^2)`, zero below `K < 1.2 + 0.7 log tau`. The earliest unstable pair merges (angular momentum conserved, mass retention 0.9) until every pair outlasts 1 Gyr.
5. **Floors.** Bodies lighter than 0.1 Earth masses (inside the snow line) or 0.60 (beyond) are debris. No chain starts, and no hot Jupiter settles, closer than the distance where a blackbody would be at the dust sublimation temperature (`DUST_SUBLIMATION_K`, 1500 K; `sublimation_au`), which keeps the luminous stars' planets out of the vapour of their own rock.
6. **Gas.** `planet_mass.gas_outcome` per core with an extra metallicity exponent.
7. **Migration.** Hot Jupiters at a share 0.0092 times 10^(2.2 [Fe/H]) of giants (period 1 to 10 d), clearing everything inside; other giants move in by e^-|N(0, 0.72)|, scattering planets within 8 mutual Hill radii. Giants move in on their own, so when two land within 8 mutual Hill radii the heavier throws the lighter out (`_scatter_crowded_giants`).

A whole-disc oligarch chain alone gave about three times too many planets (tested), which is why the inner system is its own draw.

## Calibration (re-run on the wiring with GEN.950, n = 1000 per host mass, seed 99; observed in brackets)

| Host | Inner share | Inner multiplicity | Giants 2 to 8 AU | Hot Jupiters |
|---|---|---|---|---|
| 0.3 Msun | 0.78 (0.75) | 3.0 (2.8) | 0.02 | 0.4% (0.27%) |
| 0.6 Msun | 0.75 | 3.1 | 0.07 | 0.5% |
| 1.0 Msun | 0.33 (0.30 +-0.03) | 3.2 (3.0 +-0.3) | 0.11 (0.141 +0.020 -0.018) | 1.0% (0.6 to 1.2%) |
| 1.3 Msun | 0.32 (0.35) | 2.8 (1.8) | 0.15 | 0.9% |

Hot Jupiters are compared on a mixture of metallicities (N(-0.03, 0.2) dex, the survey's), not at solar [Fe/H], because the rate rises with 10^(2.2 [Fe/H]). Spacing: median 14.6 mutual Hill radii, 93 percent of pairs over 10 (observed mode near 20, 93 percent). Giants rise 7 times between [Fe/H] -0.4 and +0.4; the 10^(1.2 [Fe/H]) law of Johnson 2010 gives 9 over that range (tweak L3). A Sun-like star holds 6.8 planets on average, 0.11 of them giants between 2 and 8 AU: GEN.950's lower solids beyond the snow line took the cold giants from 0.17 to 0.11 at the same constants, and the re-tune brought back most of it.

## Known weaknesses (tracked, not hidden)

- Giants inside 2.5 AU and the exact metallicity slope are still loose. Giants between 2 and 8 AU are 0.11 per Sun-like star against 0.141 observed, the best the search found with the 182 K line. Several tuned constants sit near the edge of their search range (`PLANET_OUTER_GROWTH_AU` 11.4); a disc-evolution model (GEN.509, GEN.510) replaces them.
- Hot-Jupiter statistics need n of several thousand to test.
- The inner masses pile up near the inner core cap (`PLANET_INNER_CORE_MAX_EARTH` 12); warm super-Earths of 1 Earth mass or more are about 0.06 per star against 0.20 observed.
- An F-type host has a multiplicity (2.9) above the observed 1.8 of Yang 2020; the multiplicity has no stellar-mass dependence yet.
- The transit multiplicity has more doubles than Zhu 2018 (0.22 against 0.17), because the harness uses its inclination law and the model's singles are intrinsically fewer.

## Stand-ins

Section 6.36 of `docs/design/research/gap-filling-plan.md`, rows PA1 to PA12, each with tier, tweak, falsifier test and what replaces it. The listing `stand-in-values.md` is regenerated from it.

## Wiring (train 4)

`StarSystem._plan_architecture` builds the `planet_architecture.Plan` for the star (the architecture cut to the planets the system holds, then further oligarchs for a system asked for more); `_generate_planets` takes the planet of slot `i` from it: its mass and formation distance, at the slot's distance (never inside the previous planet's orbit limit). `PLANETS` False gives none, `NUM_ORBITS` exactly that many, `MAX_PLANETS` True the most the disc holds (an inner system and a starting count two widths above the mean), False the fewest. SLOTS, FORCE_CLASS and HABITABLE_WORLD requests keep the formation law (`planet_mass.planet_mass`). A hot Jupiter keeps its formation distance for the composition and its final distance for the orbit. Planets inside the orbit floor of the star (`_orbit_floor_au`) are left out. A belt takes the place of a planet in its slot only where the disc built a body under `BELT_REPLACES_BELOW_EARTH` (1 Earth mass; PA13), and the architecture's own spacing is kept (the old `previous planet + 5 Hill radii` rule pushed planets outward and a census of 300 Sun-like systems showed it halving the giants); `validate_system` still enforces the stable spacing of the whole. Census after the fixes: the planets of the plan and of the finished systems agree (7.0 and 6.7 per Sun-like star; giants 0.19 in both; inner planets 59 of 59). The removed constants are `ARCH_*`, `PLANET_ASSEMBLY_SCATTER_DEX`, `GIANT_IMPACT_SURVIVAL_FRACTION`, `PLANET_SPACING_RATIO_RANGE`, `FIRST_PLANET_HZ_FRACTION_RANGE`, `INITIAL_PLANET_DISTANCE_FACTOR` and `PLANET_OCCURRENCE_REFERENCE_OBJECTS`; registry rows PM1, PM3, PM4, PM6 and PM7 are retired into PA1 and PA2. The change alters the plan stage: `--affects plan`, a new generator epoch, and the life, air and class shares move with it.

## Tests

`test_planet_occurrence.py` (reader), `test_kepler_multiplicity.py` (harness), `test_planet_architecture.py` (law, model and the statistical falsifiers), `test_prevalence_occurrence.py` with `tests/audit/prevalence_occurrence.py` (TEST.395: the census against sourced bands widened by the sampling error of n systems; a miss carries its GEN item).
