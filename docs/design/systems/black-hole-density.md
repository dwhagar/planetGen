# Black hole density: metallicity-dependent (GEN.567, design for GEN.563)

Status: built (GEN.563, see Built below); the plan below is what was decided. Boss answered the seven cards BH1 to BH7 (= Cards 1 to 7) at 12:33 to 12:35 AM EDT on 2026-10-18 (as the design lane reported; the Decisions table below), six on the recommended option and card 5 on the cap, which overrides the widest verified range. Written at 12:13 AM EDT, 2026-10-18 (as the design lane's container clock reported) by a read-only design thread (no code). Inputs: `docs/design/systems/black-hole-formation-sources.md` (GEN.560, verified; the handoff copy is `/mnt/project-files/handoffs/gen560-black-hole-formation-2026-10-18.md`), the rate sites read on main at 0c27fd8, and the numbers computed below (own computation, reproducible with the generator's own IMF and `population_densities`). Items: GEN.567 (this note), GEN.563 (the code).

## Built (GEN.563)

Code: `physics/black_hole_formation.py` (the IMF share, the formed share of a window, the failed share, the cap and the draws), `galaxy/remnant_density.py` (a place's factor in its scalar, array, profile and component forms, the galaxy's mean), the two shares as `Structure.black_hole_solar_share` and `black_hole_metal_free_share` (`galaxy/structure_draw.py`, not stored columns, so the shape carries them to every consumer and no seed has to be passed around), and the consumers of section 4. Where the build departs from the plan above:

- **The draws live in the structure.** The two shares are drawn from `unit_seed(seed, "structure", "black-hole-solar")` and `"black-hole-metal-free"`, after the neutron stars' unbound share they are capped by (`ejecta_field.unbound_fraction`). A galaxy built with no seed has the reference draws, `tuning.REFERENCE_BLACK_HOLE_DRAWS` (the middle of the metal-free range and 0.55 of the solar one, so that the reference structure holds the 1e8 black holes of the older estimates after the lighter halo of GEN.1154; row BR1 of the registry). The planned `f_cap` of 12.1% to 24.7% and the order f_sun first, f_low then `max(f_low, f_sun)` are as decided.
- **The bound for the object-first reach** is the metal-free share (no population's factor can exceed it) rather than a tighter bound at each group's lowest mean: it costs a few rejected candidates in the halo, where the factor is near that share anyway, and cannot be violated (T6 checks the lowest and highest draws).
- **The galaxy's total** (the ejecta, T5 and T9) is `remnant_density.mean_factor`: the place factor averaged over the density on a grid of 72 log radii, 24 polar angles and 6 azimuths, to a thousand disk scale lengths. The clusters' and satellites' stars (a percent) are not in it, as `galaxy_systems` leaves them out.
- **Streams.** Every class draws from `{prefix}:{kind}:{subtype}` (`phenomenon_scatter.class_stream`): per shell, per reach and far run, and per component, so the neutron star, planetary nebula and remnant rows no longer move when the black hole counts do (T2).
- **A standalone sector** (no galaxy) rolls the black holes at the failed share at solar metallicity of the middle draw, 9.3%, times `base_rate_per_star()`: 5.8e-4 per star against the old 5.0e-4.
- **Numbers at the solar circle**, Milky Way-like shape, middle draw: 7.7e-4 black holes per star, 1.5 times the old 5.0e-4; 150 sampled galaxies range over 0.44 to 2.84 times, the mean 1.46 (T9 asserts a looser range).

Stand-ins, all tunable in `tuning`: the uniform shape across each verified range (only the range is known), the linear-in-Z form between the anchors and its hold above solar, the generic cluster and satellite peaks in the scatter's weights, and the 100 Myr completeness window averaged over 64 ages and 64 progenitor masses.

## Decisions (Boss's taps on the seven cards, 12:33 to 12:35 AM EDT, 2026-10-18)

| Card | Answer | Tapped |
|---|---|---|
| 1 range and draw | Widest range, two ordered draws per galaxy (recommended) | 12:33 AM EDT |
| 2 shape | Linear in Z, held at the solar value above solar (recommended) | 12:33 AM EDT |
| 3 birth [Fe/H] | No, density only; GEN.553 draws it in 8.4 (recommended) | 12:33 AM EDT |
| 4 intermediate-mass holes | Hold today's rate, 1.5e-9 per star (recommended) | 12:34 AM EDT |
| 5 neutron star plus black hole budget | **Cap f so the sum fits (not the recommended option; see the cap below)** | 12:34 AM EDT |
| 6 mass law | Leave it for GEN.553 and GEN.554 in 8.4 (recommended) | 12:34 AM EDT |
| 7 mass islands | Record them for GEN.554 (recommended) | 12:35 AM EDT |

**The cap (card 5).** Per galaxy, `f_cap = 1 - NS_formed / IMF_SHARE`, where `NS_formed` = the neutron star rate per star (5.6e-4 / 0.14 = 4.0e-3) times (1 + the galaxy's neutron star unbound share, `ejecta_field.unbound_fraction("neutron-star", seed)`, drawn from the seed in 0.2 to 0.4). So `f_cap` is 12.1% to 24.7% depending on that share. The draws of section 2 become: **GEN.554 (slice 3):** the unbound share is now each galaxy's own (`natal-kick.md`, 0 to 74%), so the cap is the Milky Way-like galaxy's, the same for every galaxy (`structure_draw.draw_black_hole_shares`): a per-galaxy cap would be 0 where more than 59% of the neutron stars are ejected. The change of the cap's reading is on Boss's card.

```
f_sun ~ U(0, min(0.45, f_cap))
f_low  = max( min( U(0.15, 0.32), f_cap ), f_sun )
```

so neutron stars plus black holes formed never exceed the stars of 8 to 150 Msun in any galaxy (the sum reaches exactly 100% only when f_low = f_cap). Order dependency: the black hole draw reads the neutron star unbound share, so that draw comes first; both are pure functions of the galaxy seed. A later change of the neutron star rate (GEN.554) moves the cap with it; nothing else changes. Effect over 300,000 sampled draws (own computation, solar-circle mix, default galaxy of 2.953e11 systems): holes range from 0.5 to 3.1 times today's 1.5e8, mean 1.5 times (2.2e8), instead of 0.6 to 5.6 times and mean 3.1 times uncapped. **Consequence for the sources:** the cap overrides the widest verified range: solar values above 25% (Sukhbold's 26 to 45%) and metal-free values above the cap (Ebinger's up to 32%) are never reached; label the cap "derived" and the truncated ranges "verified, cut by the budget" in the registry. At the cap every massive star becomes a neutron star or a black hole, leaving none for other fates; GEN.554 may want a margin.

Boss's decision (11:51 PM EDT, 2026-10-17, GEN.563 card): the black hole density is **metallicity-dependent**, replacing `PHENOMENON_DENSITY_PC3["black-hole"] = 7e-5` per pc^3. Neutron stars are untouched.

## 1. What the density becomes

At a place, the expected black holes per star (a "star" is a system the density counts):

```
r_BH(place) = IMF_SHARE x RETAIN x sum over populations p [ w_p(place) x C_p x F_p(place) ]
```

| Term | Meaning | Value | Source state |
|---|---|---|---|
| IMF_SHARE | stars born with 8 to 150 Msun, per star born, from the generator's own IMF (`tuning.IMF_BREAKS_SOL`, `IMF_SLOPES`: the Chabrier 2003 function of `galaxy/stellar-mass-function.md`) | 0.7127% (0.6372% on the Kroupa function before GEN.927) | derived (own computation; reproduces GEN.560's 0.637%) |
| w_p | the population's share of the stars at the place: `galaxy.density.population_densities` (young, intermediate, old, bulge, halo, cluster, satellite and its two young windows, nuclear) | from the galaxy | existing code |
| F_p | the mean of f(Z) over population p's [Fe/H] distribution at the place (section 3) | 0 to 45% | verified range, see 2 |
| C_p | share of the population's window in which its remnants have formed: 1 for every window that starts after 100 Myr (an 8 Msun star lives 45 Myr), **0.79** for the young window (0 to 0.1 Gyr, uniform birth, IMF-weighted lifetimes from `stellar_evolution.main_sequence_lifetime_gy` times `remnant_age_fraction`); also divide by (1 - IMF_SHARE) = 1.0064 because the generator's star counts leave out stars already collapsed (`evolve_star` returns None for them) | 0.79 young, 1.0064 others | derived (own computation) |
| RETAIN | black holes that stay bound: 1 - `BLACK_HOLE_UNBOUND_FRACTION` | 0.98 | verified through `compact-remnant-regions.md` (Sweeney et al., arXiv 2210.04241, "2% of black holes leave"); the primary paper was not opened here |

Why sum over populations: the three disk windows (young, intermediate, old) are the partition of the disk's star formation history (`disk_history`), so the stars born in all windows are the stars that made the remnants at the place; the bulge, halo, cluster and satellite windows are the old ones. Their birth [Fe/H] distributions are the GEN.158/159 ones (`physics/metallicity.py`, `tuning.METALLICITY_BY_POPULATION`, `METALLICITY_HALO`, and the cluster and satellite peaks of `ComponentField.fe_h_peaks`).

The existing GEN.132 place factors stay (vertical thickening `REMNANT_SCALE_HEIGHT_RATIO["black-hole"]` = 2.5, capped at 4; core excess `BLACK_HOLE_CORE_EXCESS` = 1 over 2 kpc). They move black holes after birth (kicks), and the core excess has the opposite sign to the metallicity effect (the inner galaxy is metal rich, so lower f), so they do not double count. Nothing else in `remnant_distribution.py` changes.

## 2. The anchors and the draw (cards 1 and 2)

Sourced anchors (all from GEN.560, section 3 and the table in 4; the share is of core-collapse events, stars of 8 to 150 Msun, Salpeter IMF):

| [Fe/H] | O'Connor and Ott 2011 | Ebinger et al. 2020 (two calibrations) | Sukhbold et al. 2016 (four calibrations) | Widest range |
|---|---|---|---|---|
| solar (0) | 0 to 7% | 5 to 8% and 16 to 21% | 26 to 45% | **0 to 45%** |
| 1e-4 solar or zero (below the generator's -3.5 clamp) | ~15% | 18 to 20% and 27 to 32% | none (solar only) | **15 to 32%** |
| above solar | none | none | none | none |

Observations (Neustadt 2021: 0.16, 90% interval 0.04 to 0.39; Byrne and Fraser 2022: upper limits 0.23 to 0.61) exclude nothing in these windows except the last few points of 45% against Neustadt's 90% upper edge.

Design (Card 1 and Card 5 decided; the draws below are the uncapped ranges, the Decisions block shows them with the cap applied): two numbers per galaxy, drawn once from the galaxy seed (the precedent is `ejecta_field.unbound_fraction`, a seeded draw per galaxy; the galaxy laws table of GEN.368 is where they are listed), each uniform over its widest verified range because only the range is known:

```
f_sun ~ U(0, 0.45)       (solar-metallicity failed-supernova share)
f_low ~ U(0.15, 0.32)    (metal-free share), then f_low := max(f_low, f_sun)
f(Z)  = f_low + (f_sun - f_low) x min(10^[Fe/H], 1)     # Z measured in solar units; held at f_sun above solar
```

The `max` keeps "more at low metallicity" true in every galaxy: no source gives a lower share at low Z than at solar (GEN.560 section 4), and without it a galaxy with f_sun = 40% and f_low = 20% would have fewer holes in the halo than in the disk. Cost: where f_sun > f_low draw the metal-free value is f_sun itself, above 32% at most up to 45%, a value no low-Z source gives but a verified solar source implies.

The shape between the anchors is a model choice (no source gives a third point). Linear in Z (not in [Fe/H]) is chosen because Ebinger's own grid is flat between Z = 0 and Z = 1e-4 solar (18 to 20% and 27 to 32%): the share saturates toward low Z, which linear-in-Z reproduces and linear-in-[Fe/H] does not (it would keep rising to the clamp). Label: bracket (anchors verified, shape a stand-in, tunable). Above solar nothing is sourced, so f is held at f_sun (stand-in).

Expected values at the Sun with these laws (own computation, `population_densities` on a Milky Way-like shape of 2.6 kpc by 300 pc, bulge 1.58 kpc; the real default galaxy will differ in detail):

| Place | f_sun 22.5%, f_low 29% (middle draw) | f_sun 5%, f_low 15% (low draw) | f_sun = f_low = 40% (high draw) | Today (7e-5 / 0.14 = 5.0e-4) |
|---|---|---|---|---|
| solar circle | F 0.247, 1.57e-3 per star (x3.1) | 0.085, 5.4e-4 (x1.1) | 0.40, 2.55e-3 (x5.1) | 5.0e-4 |
| inner disk 3 kpc | 0.234, 1.49e-3 (x3.0) | 0.065, 4.1e-4 (x0.8) | 0.40, 2.55e-3 | |
| bulge | 0.239, 1.52e-3 (x3.0) | 0.071, 4.5e-4 (x0.9) | 0.40, 2.55e-3 | |
| thick disk, 1 kpc up | 0.250, 1.59e-3 (x3.2) | 0.089, 5.6e-4 (x1.1) | 0.40, 2.55e-3 | |
| inner halo, 15 kpc | 0.286, 1.82e-3 (x3.6) | 0.145, 9.2e-4 (x1.9) | 0.40, 2.55e-3 | |
| outer halo, 40 kpc | 0.288, 1.84e-3 (x3.7) | 0.149, 9.5e-4 (x1.9) | 0.40, 2.55e-3 | |

(Per star = IMF_SHARE x F, before RETAIN, C_p and the GEN.132 place factors.) Two findings Boss should see before approving:
- **The level matters far more than the metallicity pattern.** The disk, bulge and thick disk differ by only 5 to 8% from each other at a given draw; the halo is 15 to 20% richer at the middle draw (up to 2.3x at the low draw). The draw moves the whole galaxy: in the default galaxy (2.953e11 systems) the holes range from 9e7 to 8.3e8 over 200,000 sampled draws (0.6 to 5.6 times today's 1.5e8), with a mean of 4.5e8 (3.1 times). Today's 5.0e-4 per star (f = 7.8%) is near the low end of the range. The floor is not zero because the metal-free value is at least 15%, so even a galaxy with f_sun = 0 still makes holes in its metal-poor old stars and halo.
- **The old 7e-5 stays inside the new range** (as GEN.560 found), so it stays a test target (section 7), not an error.

### Card 1: the range and the draw
- **Widest range, two ordered draws per galaxy (recommended).** f_sun uniform over 0 to 45%, f_low over 15 to 32%, kept in order. Consequence: follows the source rule as the coordinator stated it (widest verified range), nothing hidden; galaxies range from 0.6 to 5.6 times today's black hole count, the mean is 3.1 times (4.5e8 in the default galaxy).
- **Newest only (Ebinger 2020).** f_sun uniform over 5 to 21%, f_low over 18 to 32%, same ordering. Consequence: reads "when two verified sources are in doubt, the newest wins" and drops Sukhbold and O'Connor and Ott; the mean falls to 2.1 times today's (3.1e8 in the default galaxy) and the range to 1.2 to 3.1 times.
- **One fixed middle curve (22.5% and 29%).** Consequence: simplest to test and reproduce, but it is a point where the sources give a range, against the standing rule.

### Card 2: the shape between the anchors
- **Linear in Z, held at the solar value above solar (recommended).** Consequence: reproduces Ebinger's flat low-Z plateau; a bulge or inner disk above solar gets f_sun.
- **Linear in [Fe/H].** Consequence: simplest in the log metallicity the code stores, but keeps rising down to -3.5 where Ebinger's grid has already saturated; at the middle draw the disk gets about 7% less and the halo about 11% less than the recommended form.

## 3. Birth [Fe/H] of a remnant (card 3)

A remnant is born at the end of its progenitor's few Myr, in the population the place's stars were born into, so its birth [Fe/H] is the population's at that place (GEN.158/159), not the present gas. The density does not need a per-object value: F_p is the closed-form mean of f over the population's distribution, `metallicity.expected_factor(function, peaks)` (five-point Gauss-Hermite on each normal, already used by the GEN.236 planet cull in `galaxy/profile.py:290-292`), so it stays a pure function of the place. The linear-in-Z form has a kink at solar; the five-point rule is within 2% of a fine integral there (own check over means -1.0 to +0.6 and scatters 0.1 to 0.4).

Two places use it, and they must use the same function: the scatter's longitude-bin weights (arrays over bins, using `metallicity.mean_fe_h` per population; the cluster and satellite columns use the generic `tuning.METALLICITY_CLUSTER_PEAKS` and `METALLICITY_SATELLITE_PEAKS` as constants because those populations are a percent of the systems and the per-component peaks need a per-cell lookup that has no array form), and the sector profile's `kind_factors_at` (exact peaks from `ComponentField.fe_h_peaks`).

### Card 3: does a black hole row carry its own [Fe/H]?
- **No, the density only (recommended).** Consequence: no new column and no migration here; the per-object birth [Fe/H] is drawn in GEN.553 (8.4, "origin, birth, state and [Fe/H]") weighted by f(Z) times the population's [Fe/H] density, so metal-poor tails are over-represented consistently with this note.
- **Yes, store it on the scatter rows in 8.1.** Consequence: one more column on the scatter table (an Alembic revision, in the same wiped rescatter), black hole objects show a [Fe/H] in 8.1, and GEN.553 reads it instead of drawing.

## 4. The three consumers and the exact sites (as read on 0c27fd8)

| Site | Line(s) | Today | GEN.563 change |
|---|---|---|---|
| `tuning.py` | 3437 `"black-hole": 7e-5`; 3610-3614 `phenomenon_rate_per_star`; 3554/3559 unbound fractions | one constant | the constant and its comment go; `phenomenon_rate_per_star("black-hole")` is gone (the kind has no single rate), a new function returns the place rate; the old 7e-5 is kept only as a test constant |
| `phenomenon_scatter.class_rate_per_star` | 158-169 | `phenomenon_rate_per_star(kind)` split by the intermediate-mass chance | black hole classes take the place rate in the weights instead of a constant per star (see 5) |
| `phenomenon_scatter._scattered_rates` | 229-239 | `(kind, subtype, rate)` | the black hole rate becomes a weight, so rate = 1 x IMF_SHARE x RETAIN x galaxy draws, and the per-bin `sum w_p C_p F_p` joins the weights below |
| `phenomenon_scatter._band_rows` | 292-369, weights at 324-336 | `base x vertical x radial` | `base` becomes `sum_p population_densities_array[p] x C_p x F_p` for the black hole table (neutron stars keep `base`) |
| `phenomenon_scatter.below_cut_draws` | 415-443, mean at 436-438 | `class_rate_per_star x expected_stars x below x placement_factor` | multiply by the place factor `sum w_p C_p F_p` at the sector centre (the profile's black hole factor, one call) |
| `phenomenon_scatter.shell_expected`, `_expectation_factor` | 446-477 | remnant factors only | add the place factor, so the progress bar matches |
| `phenomenon_reach.class_fractions`, `class_bounds`, `_weights` | 43-67, 70-86, 89-102 | `base x remnant factors`, bound by `extent_majorants` | the object-first bound needs `F_p` at its largest in the group: f is non-increasing in [Fe/H], and the disk mean falls with radius, so bound each population by F at the lowest mean the group reaches (the same trick `class_bounds` plays for the radial factor); `_weights` returns the new weight; the existing `check` flag asserts weight <= bound |
| `galaxy/profile.kind_factors_at` | 256-300 (`"black-hole": 1.0` at 296) | 1 everywhere | `place factor / reference`, where the run_sector rolled path and the below-cut draw read it |
| `run_sector.generate_sector_phenomena` | 223 | `phenomenon_rate_per_star(kind) x star_count x kind_factors[kind]` | the black hole kind uses the new rate (standalone sectors keep factor 1 through the solar circle's rate) |
| `galaxy/ejecta_field.unbound_count` | 59-62 | `phenomenon_rate_per_star(kind) x galaxy_systems x unbound_fraction` | the galaxy total formed is the integral of the place rate over the galaxy (a cached function like `galaxy_systems`), times `BLACK_HOLE_UNBOUND_FRACTION`; bound plus ejected equals formed (today bound = rate and ejected = rate x 2%, a 2% double count) |

About 10 functions in 5 files (`tuning`, `phenomenon_scatter`, `phenomenon_reach`, `profile`, `ejecta_field`, plus the rolled path in `run_sector`). The neutron star classes read none of this: `NEUTRON_STAR` weights, rates and `row_density_divisor` stay as they are.

**Stream hazard.** Every shell draws all classes from one stream in `SCATTER_CLASSES` order, black hole first (`phenomenon_scatter.py:63-68` says a new class goes last so earlier ones keep their draws). Changing the black hole counts therefore moves the realised neutron star, planetary nebula and supernova remnant rows even though their laws are unchanged. GEN.563 should give each class its own stream (`f"{seed}:phenomena:shell:{shell}:{kind}:{subtype}"`): every class's rows move once, in the one wiped rescatter this is going into anyway, and the neutron star rows are then bit-identical across any later black hole change. Test T2 below needs it.

**What GEN.563 must also not break:** the intermediate-mass hole chance (card 5), the companion chances in `multistar.MEMBER_STATE_BY_PLACED` (174-186; per placed object, with their own sources, left as they are), and the scatter mass cut (14 Msun, card 6).

## 5. Intermediate-mass black holes and the neutron star budget (cards 4 and 5)

**Intermediate-mass holes** are drawn as `BLACK_HOLE_INTERMEDIATE_MASS_CHANCE` = 3e-6 of the black hole rate (`class_rate_per_star`, GEN.259: "about 300 in a galaxy of 1e8 black holes"). They are not remnants of failed supernovae, so f(Z) should not scale them: at the middle draw they would go from about 440 to about 1,200 in the default galaxy.

### Card 4: the intermediate-mass hole count
- **Hold their rate at today's (recommended).** Their rate per star becomes 3e-6 x 5.0e-4 = 1.5e-9 (the old black hole rate), independent of f(Z). Consequence: GEN.259's "hundreds" stays true; they are no longer a fraction of the stellar black holes.
- **Keep the 3e-6 share.** Consequence: they scale with the black hole count, 0 to about 2,500 across galaxies, about 1,200 at the middle draw.

**Neutron star budget.** Neutron stars are untouched, and the code reads their rate as the formed count: bound at 4.0e-3 per star and ejected on top at 20 to 40% (`ejecta_field.unbound_count`), i.e. 0.48 to 0.56% of stars, against 0.637% of stars born with 8 to 150 Msun. Formed neutron stars plus formed black holes cannot exceed the massive stars, which leaves room for f of only 25% (ejection 20%) down to 12% (ejection 40%); over the draws, the sum overshoots the massive stars in 47% of galaxies at 20% ejection, 65% at 30% and 85% at 40% (widest range; with the newest-only range 0%, 37% and 92%), by up to 0.21% of stars at the top (45% x 0.637% + 0.56% = 0.85%, 1.3 times the budget). The inconsistency exists in the current constants (neutron stars were tuned to the observed bound count, not to the mass function); the black hole range just makes it visible. GEN.554 ("remnant counts follow from the mass function and the star formation history") is the item that reconciles them.

### Card 5: the budget (decided: cap f; the recommended option below was not the one chosen)
- **Leave it, register the gap, add a report test (recommended).** The gap registry (section 6 of `gap-filling-plan.md`) lists "neutron star plus black hole budget exceeds the massive-star share in about half to two thirds of galaxy draws; reconciled in GEN.554"; a test prints the overshoot and fails only if it grows. Consequence: follows both of Boss's decisions (widest range, neutron stars untouched); no galaxy is blocked.
- **Cap f so the sum fits.** f_cap = 1 - NS_formed / IMF_SHARE (37% with no ejection, 12 to 25% with ejection). Consequence: the sum can never exceed the massive stars, but it cuts the sourced upper range and ties the black hole density to the neutron star one, against the widest-range rule.

## 6. The mass law and the mass islands (cards 6 and 7)

Today a stellar-mass black hole's mass is uniform over 5 to 20 Msun (`BLACK_HOLE_MASS_RANGE_SOLAR`, mean 12.5). It is read at: `phenomenon_scatter.py:154` (the law, and `share_above`, the share above the 14 Msun cut: `(20 - c)/15` = 40%), `compact_remnant.py:316` (the draw), `multistar.py:476` (the companion draw, with progenitors 20 to 100 Msun, and `DIRECT_COLLAPSE_PROGENITOR_SOL` 40), `place_terms.py:99` (uniform moments for place terms), `phenomena_plausibility.py:388`, `web/lib/classref.py:314`, plus tests (`test_place_terms.py:90` asserts a mean between 9 and 14.2; `test_phenomena.py:64`, `test_phenomenon_scatter.py:680-720`). That is 6 source sites and 3 test files.

The verified masses (GEN.560 section 4, point 3): IMF-weighted mean 7.7 to 9.2 Msun if the hole keeps only the helium core, 13.2 to 13.8 if it keeps the whole star (Sukhbold et al.; the authors prefer the helium core), lower fallback holes at 3.2 to 7.3, gravitational mass 0.9 of baryonic, heavier at low Z. The present mean, 12.5, lies between the two readings, so nothing in the density change forces the mass law to move.

### Card 6: the mass law in 8.1
- **Leave it for GEN.553 and GEN.554 in 8.4 (recommended).** Consequence: the GEN.301 rescatter keeps today's masses; 8.4 changes the mass law once, together with the progenitor and the initial-final relation, which rewrites the black hole rows' masses and the number above the cut a second time (a mass-only change to rows that are already stored; GEN.564 refill covers the sector fill, but the scatter rows above the cut need their share and masses redrawn).
- **Change it in 8.1.** Consequence: one rescatter for density and masses together, but 8.4 then changes it again; a mean near 8 Msun puts almost every hole below the 14 Msun cut (the scatter keeps only a thin tail and the sector fill draws the rest), so the cut, `share_above`, `place_terms` and the progress estimates all move with it; the form of the new law is a further card.

### Card 7: the mass islands
The islands (which progenitor masses fail: from about 20 Msun, mostly failure above 35 Msun at low Z, per Sukhbold, O'Connor and Ott, Maltsev et al. 2025) change the progenitor and the mass of each hole, not the share f: the density above already integrates over the whole IMF. They matter where the generator reads a progenitor: the direct-collapse flag (`DIRECT_COLLAPSE_PROGENITOR_SOL` = 40, which GEN.560 says is not verified at solar Z where 23 to 25 Msun also fail) and `_remnant` (`multistar.py:466`, progenitors 20 to 100).
- **Not in 8.1, record them for GEN.554 (recommended).** Consequence: no new rules in 8.1; GEN.554 gets the islands as its mass-law input, with a test that the island windows integrated over the IMF land inside this note's f range.
- **Add a progenitor window in 8.1.** Consequence: black hole progenitors are drawn from the island windows (carbon-oxygen core rule of Maltsev et al.: 6.6 to 7.2 and above 13 Msun at solar, shifting down with Z), 3 more functions touched in `multistar.py` and a metallicity input into the companion draw; a larger change than the density, with only abstract-level support for the windows' IMF mapping.

## 7. Tests (what GEN.563 must show)

- **T1 (density against the derived share).** At three places with pinned draws (f_sun, f_low): the disk at the solar circle, the bulge, the halo, the scatter's expected black hole count and the sector fill's mean match `IMF_SHARE x RETAIN x sum w_p C_p F_p x place factors` within Poisson error; `shell_expected` agrees.
- **T2 (neutron stars unchanged).** With the black hole rate changed between two values, the neutron star rows (and planetary nebula and supernova remnant rows) are bit-identical (needs the per-class streams of section 4); the neutron star weights, rate and `row_density_divisor` are untouched.
- **T3 (monotone).** For every draw in a seed sweep, F at [Fe/H] -3.5 is at least F at solar, and the halo's rate is at least the disk's.
- **T4 (the old constant is a target that shows the gap).** `f = 7.8%` reproduces 7e-5 per pc^3 at the reference density 0.14 (per star 5.0e-4); the new middle is reported as a multiple of it; the test asserts 7e-5 lies inside the new range and prints the ratio at the three places.
- **T5 (ejecta).** Bound plus ejected equals the galaxy total formed; ejected = 2% of it.
- **T6 (bound).** The `check` flag of `_group_rows` raises on no group in a galaxy with the lowest and highest draws.
- **T7 (draw).** Over 1,000 seeds f_sun stays in 0 to min(45%, f_cap), f_low >= f_sun, f_low <= max(f_cap, f_sun), f_cap lies in 12.1% to 24.7%, and a seed reproduces its pair.
- **T8 (intermediate-mass holes, if card 4 recommended).** Their count in the default galaxy stays at about 440.
- **T9 (budget cap).** Over the seed sweep, neutron stars formed plus black holes formed is at most the 8 to 150 Msun share of stars in every galaxy, and equals it only where f_low = f_cap; the galaxy's black hole count stays between 0.5 and 3.1 times the old constant.

## 8. What goes into the one wiped GEN.301 rescatter

- The black hole scatter rows (stellar-mass and intermediate-mass classes): counts change by the draw (0.5 to 3.1 times today's with the budget cap), positions follow `sum w_p C_p F_p`. Every scatter class's rows are redrawn once because of the per-class streams.
- Below-cut black holes in every sector's fill (stage `BELOW_CUT`) and the rolled black holes of standalone sectors: the sector stream's draw count changes; **generator epoch bump** (the next free number in `src/tests/golden/epochs.json`; one pending epoch rides until a release stamps it) and the golden digest changes. No schema change unless card 3's second option is chosen (then one Alembic revision, next free at merge, with the others).
- Ejecta (`SectorProfile.ejecta` mean counts for black holes) change with the galaxy total.
- Black hole masses: **unchanged** under the recommended card 6; a second pass in 8.4 if it is left.
- Unchanged: neutron stars (after the stream split), planetary nebulae, supernova remnants, hypervelocity stars, the nucleus.

## 9. Gap registry stand-ins (labels, all tunable, replaced when verified)

| Number | State |
|---|---|
| IMF_SHARE 0.713% | derived from the generator's IMF |
| f_sun range 0 to 45%, f_low range 15 to 32% | verified ranges; the uniform shape across them is a stand-in (only the range is sourced) |
| f_low >= f_sun ordering, f_low up to 45% | bracket (consistency with a verified solar value; no low-Z source above 32%) |
| linear in Z between the anchors, held above solar | stand-in (anchors verified, shape and the above-solar hold not) |
| anchor at the -3.5 clamp for Ebinger's 1e-4 solar ([Fe/H] about -4) | derived (a small shift) |
| C_young 0.79, 1/(1 - IMF_SHARE) | derived |
| RETAIN 0.98 | verified through a secondary note, primary not opened |
| budget cap f_cap = 1 - NS_formed / IMF_SHARE (12.1% to 24.7%) | derived; it cuts the verified ranges (Boss, card 5); GEN.554 revisits it when the neutron star rate follows the mass function |
| generic cluster and satellite peaks in the scatter's array path | stand-in (the sector profile uses the exact component peaks) |

Not verified here (carry to GEN.563 or GEN.554): the primary paper behind the 2% retention (Sweeney et al.); Luo 2025, Zapartas 2021, Patton and Sukhbold 2020 are abstract-level in GEN.560; Maltsev 2025 has no IMF map in the text read; whether Sukhbold's four calibrations and Ebinger's two are independent enough to treat their ranges as additive evidence (they are not independent physics, GEN.560 section 3).
