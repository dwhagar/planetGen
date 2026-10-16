# Where each star type and phenomenon belongs: what is built and what is left (GEN.103)

Boss (2026-10-07 11:47Z) asked for the tiered structure of where stars of different types sit in the galaxy, and for the
kinds of phenomena with similar restrictions, added to the density and type draws. This note is the review of
2026-10-16 turned into a table: which placement rules the code already has, and which three population-dependent rules
for the type draws are still missing, each with a way to calculate it and what it needs from Boss.

Nothing here is a fixed catalogue. Every number a rule uses comes from the sector profile (`galaxy/profile.py`), which is
a function of the galaxy shape, the cell's address and the galaxy seed. A rule below either reads a profile field or
names the physical law that turns profile fields into the number. Where a law needs a constant from the literature the
note says so, marks it for Boss's research, and the code reads it from `tuning.py` with its source named.

Informs: GEN.103, GEN.269, GEN.270, GEN.271, GEN.272, GEN.273, GEN.275.
Status: design, 2026-10-16. Evidence tags: [S] read from the repository, [R] recalled and unconfirmed, [C] computed here.

## 1. What decides where a thing goes

A cell's profile carries nine populations (`profile.populations`, each with a density, a share of the cell's stars, an age
window, a metallicity peak, a velocity dispersion and a rotation rule) and a `kind_factors` map. Two kinds of rule sit
on top of them, and the first is complete.

| Layer | Question it answers | Where it lives | Status |
|---|---|---|---|
| Density | How many systems does the cell hold? | `galaxy/density.relative_density`; `profile.expected_systems` | built |
| Population mix | Which populations are they from? | `galaxy/density.population_densities`; `profile.mix`, `profile.populations` | built (GEN.234) |
| Age and metallicity | How old, how metal rich? | `tuning.STELLAR_POPULATION_AGE_RANGES_GY`; `profile._fe_h`; `bright_stars.py` reads `populations[...].fe_h` | built (GEN.158, GEN.234) |
| Per-kind rate factors | How much of a phenomenon does the place scale? | `profile.kind_factors_at`; `run_sector.py` multiplies `phenomenon_rate_per_star` | built (GEN.236) |
| Scatter placement | Where do the single objects go? | `phenomenon_scatter.py`, one weight per band and longitude bin | built (GEN.278) |
| Type draws by population | Does the place change how many stars a system holds and how far apart they sit? | `generation/multistar.py` | **not built** (section 3) |

## 2. Built: the placement table

Each row is a kind the registry (`generation/phenomena/registry.py`) already places, with the rule that decides where. The
"Follows" column is what the draw is weighted by; no row has a position or a count written into it.

| Kind (registry gate) | Follows | Why | Where it is in the code |
|---|---|---|---|
| star-system (stellar density) | total density; the population draws the age and [Fe/H] | the field stars are the populations | `density.population_densities`, `star_population.py` |
| bright-star (stellar density) | density, population by population | the bright tail of each population's mass function | `bright_stars.py` (pre-placed, per mass cell) |
| black-hole, neutron-star (stellar density) | stars, factor 1 | born where massive stars were, then kicked into the disk's thickness | `kind_factors` 1; `phenomenon_scatter` band density |
| neutron-star, magnetar (young stars) | the young population and its arm contrast, none in the halo or a cluster | too young to have drifted | `SCATTER_CLASS_GATES`, `phenomenon_scatter.py` |
| nebula (hot stars) | the cold gas fraction around the cell's own hot stars | emission and reflection nebulae need a hot star and gas | `kind_factors["nebula"]` |
| molecular-cloud (nebula field) | cold gas fraction times the young share against the solar circle's | the places that form stars | `kind_factors["molecular-cloud"]` |
| supernova-remnant (stellar density) | the gas fraction, none in the halo or in clusters, kept inside a superbubble | the shell shines only where it meets gas | `kind_factors["supernova-remnant"]` |
| planetary-nebula, brown-dwarf (stellar density) | stars, factor 1 | every population makes both | `kind_factors` 1 |
| rogue-planet, comet (stellar density) | the planet occurrence at the cell's [Fe/H], averaged over its populations | planets favour metal-rich stars | `metallicity.planet_occurrence_factor` |
| hypervelocity-star, runaway-star | the central black hole; the stars' own speeds | not a population effect | `wanderer` pass; `flag_fast_stars` |
| quasar | the nucleus | one per galaxy | `nucleus` pass |

The cluster and satellite populations (GEN.160, GEN.237) are placed by the component field, not by a gate: their density
is the components' systems in the cell's cube.

## 3. Left: three rules for the type draws

`multistar-and-compact-systems.md` section 6.6 names three effects of the population on a system's stars. None is in
`multistar.py` yet: the star count, the pairing and the axis limit read the star's mass and nothing about the place.
Each rule below is a function of fields the profile already holds. The constants are the only things to look up, and
they sit in `tuning.py` with a source.

### 3.1 Close pairs rise as metallicity falls

Observed: the fraction of solar-type stars with a companion inside about 10 AU rises steeply toward low [Fe/H]
(Moe, Badenes and Kratter 2019 [R], with the Raghavan 2010 solar-neighbourhood value as the high-metallicity end [R]).
Rule: the chance that a drawn pair is a close one (an axis under 10 AU) is a function of the star's [Fe/H], which the
cell already draws per population (`PopulationProfile.fe_h`, peak and scatter). The thick disk, bulge and halo carry
low peaks, so they get more close pairs without any place being named.

How: a monotone interpolation through the observed points, so the rule reads the star's own [Fe/H] and no population
name. **Needs from Boss:** the observed close-binary fraction against [Fe/H] (one table of 4 to 6 points).

### 3.2 Massive multiples sit in young regions

Observed: high-mass stars are far more often in multiples than low-mass ones, and the most massive sit in clusters and
associations. The mass dependence is already in `BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS` and
`STAR_COUNT_PROBABILITIES`; the place dependence is not. The `prevalence` argument of `draw_member_tree` already scales
the chance of five or more stars, so the place only has to set it.

Rule: `prevalence` is the young share against the solar circle's, clamped to 0 to some cap, where the young share is the profile's
`populations["young"].share` over the reference `profile._reference_young_share(shape)` that `kind_factors_at` already
uses. That makes the knob 1 where a cell looks like the solar circle, above 1 in arms and associations and 0 in the halo. The cap on the rise is the one constant (a design default for Boss, since the five-or-more column is a tail).
The five-or-more columns are the ones that scale, which is how the knob works today. 
### 3.3 Wide pairs survive less where it is dense

Wide pairs are pulled apart by the Galactic tide and by passing stars. The first limit is a calculation: a pair of total
mass M in a medium of mass density rho has a tidal (Jacobi) radius of about (M / (3 rho))^(1/3). The cell gives rho:
its systems per cubic parsec times `tuning.MEAN_SYSTEM_MASS_SOLAR`, plus the gas and dark matter the profile holds
(`gas.mass_solar / volume_pc3`, `dark_matter.density_solar_pc3`) [C].

| Place (today's default shape, stars only) | rho (solar masses per pc^3) | Tidal radius for a 2 solar mass pair | Against the 10,000 AU cap |
|---|---|---|---|
| 8.2 kpc out, in the plane [C] | 0.048 | 2.4 pc, 500,000 AU | cap does not bind |
| bulge, 500 pc from the centre [C] | 6.9 | 0.46 pc, 95,000 AU | cap does not bind |
| bulge, 100 pc from the centre [C] | 77 | 0.21 pc, 42,000 AU | cap does not bind |
| 10 pc from the centre [C] | 660 | 0.10 pc, 21,000 AU | cap does not bind |
| 2 pc from the centre [C] | 61,000 | 0.022 pc, 4,600 AU | cap shrinks the axis |
| a globular cluster core [R] | 10^3 to 10^5 | 0.09 to 0.009 pc, 19,000 to 1,900 AU | cap shrinks the axis |

(The density is the model's `relative_density` times `E / edge^3` times `tuning.MEAN_SYSTEM_MASS_SOLAR`, 0.4; it leaves out
the gas and dark matter, which would shorten the radii a little. The radii follow from the formula. The
cluster row is an observed range, to be replaced by the cluster's own stored density when the component carries one.)

So the existing 50 to 10,000 AU range (`WIDE_BINARY_SEPARATION_*`) is right almost everywhere and the rule only has to
cut it in dense cores. Rule: the widest axis a wide pair may have is the smaller of the 10,000 AU cap and the tidal
radius of the cell's density, with no floor above the 50 AU lower edge. The second limit, disruption by passing stars,
scales with the local stellar density times the velocity dispersion the profile already holds; its constant is the
one thing to look up. **Needs from Boss:** the encounter disruption limit (a hard-soft boundary form, with the
constant) from the binary-survival literature; the tidal term can go in without it.

## 4. How the rules reach the draw

`run_sector.py` already passes `profile.kind_factors` into the fill. The three rules need the profile's populations and
density in the system draw, so the system config gets one more input, a small `PlaceTerms` record built once per cell
from the profile (the [Fe/H] per population, the young-share prevalence, the density in solar masses per cubic parsec),
read by `multistar.draw_member_tree` (section 3.2 through `prevalence`), the pair-pick (3.1) and `wide_cap` (3.3). It is
a function of the cell, so the same seed gives the same stars wherever a run lands (the GEN.39 rule). It stays out of
`tuning.py`'s stellar evolution constants, so Core 3's evolution work is not touched.

## 5. What closes GEN.103

- 3.2 and the tidal term of 3.3 need no new science and can be built first.
- 3.1 waits for the [Fe/H] table, and the encounter term of 3.3 for its constant (both asked of Boss).
- GEN.269 and GEN.270 to GEN.275 are the multi-star items that own the code (Core 3), so this note is their input: the
  build goes there, and GEN.103 closes when they have taken the three rules in, or the rules are written as items of
  their own.
