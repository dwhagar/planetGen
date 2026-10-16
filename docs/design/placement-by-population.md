# Where each star type and phenomenon belongs: what is built and what is left (GEN.103)

Boss (2026-10-07 11:47Z) asked for the tiered structure of where stars of different types sit in the galaxy, and for the
kinds of phenomena with similar restrictions, added to the density and type draws. This note is the review of
2026-10-16 turned into a table: which placement rules the code already has, and which three population-dependent rules
for the type draws are still missing, each with a way to calculate it and, for two of them, the constants Boss supplied on 2026-10-16.

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

Moe, Kratter and Badenes (2019, ApJ 875, 61) put the intrinsic fraction of solar-type stars with a companion inside 10 AU
at **53% +-12% at [Fe/H] -3.0, 40% +-6% at -1.0, 24% +-4% at -0.2 and 10% +-3% at +0.5** (the values in the paper's abstract, as
the search summary quoted them) [R]. Boss's research (2026-10-16) gave the same trend and the same source, with a rounder
"50 to 60% at -1.0 or lower, 10 to 15% at +0.3 or higher"; the paper's own point at -1.0 is 40%, so this note uses the
paper's four points, not the round figures.

Rule: the chance that a drawn pair is a close one (an axis under 10 AU) scales with the star's [Fe/H], which the cell
already draws per population (`PopulationProfile.fe_h`, peak and scatter). The existing multiplicity tables
(`BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS`, `STAR_COUNT_PROBABILITIES`) keep their mass dependence; the metallicity
enters as a factor `F(fe_h) / F(0)`, where `F` is a straight line between the four points above (in [Fe/H]), held flat
past the ends. `F(0)` is about 20%, so a star at -1.0 gets about twice the close pairs of a solar one and a star at +0.5 about
half. The thick disk, bulge and halo carry low peaks, so they get more close pairs without any place being named.

What the sources say about the limits of the rule, and what could not be confirmed:

- Scope: the points are for solar-type primaries (about 0.6 to 1.2 solar masses). Applying the factor to every mass is a
  design choice, not an observation; Moe et al. argue the OB-star close fraction does not change with metallicity, so the factor
  should fade toward 1 for primaries above about 2 solar masses. [R]
- Wide pairs: Boss's text says the dependence vanishes beyond about 200 AU (El-Badry and Rix 2019, Gaia DR2). A later
  Gaia paper found a strong, non-monotonic wide-binary dependence on metallicity (arXiv 2010.02920, Hwang et al. 2021)
  [R, found by search, not read]. The rule therefore changes the close fraction only and leaves the wide fraction as it is
  until Boss decides which to follow.
- Mazzola et al. (2020, arXiv 2007.09059) find the close fraction tracks alpha abundance more closely than iron; the
  generator has no alpha abundance, so [Fe/H] is the proxy [R].
- Citation fixes: the El-Badry and Rix paper is dated 2019 in the search results, not 2018 as in the pasted list.

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

Two limits cut the widest axis a pair can keep. Both are functions of the cell's density, and the second also of the stars' age.

**The tidal limit.** A pair of total mass M in a medium of mass density rho has a tidal (Jacobi) radius of about
(M / (3 rho))^(1/3). The cell gives rho: its systems per cubic parsec times `tuning.MEAN_SYSTEM_MASS_SOLAR`, plus the gas and dark
matter the profile holds (`gas.mass_solar / volume_pc3`, `dark_matter.density_solar_pc3`). For the solar neighbourhood
(a total dynamical density of about 0.1 solar masses per pc^3 [R]) and a 2 solar mass pair this gives 1.9 pc, which agrees
with the 1.7 pc Jiang and Tremaine (2010, MNRAS 401, 977; arXiv 0907.2952) quote from the Oort constants for a pair of solar-mass
stars and with Boss's 1.5 to 2 pc [R]. The Oort form (M / (2 A (A - B)))^(1/3) needs the galaxy's rotation, which a
cell does not carry, so the local-density form is the one the code can use.

| Place (today's default shape, stars only) | rho (solar masses per pc^3) | Tidal radius for a 2 solar mass pair | Against the 10,000 AU cap |
|---|---|---|---|
| 8.2 kpc out, in the plane [C] | 0.048 | 2.4 pc, 500,000 AU | cap does not bind |
| bulge, 500 pc from the centre [C] | 6.9 | 0.46 pc, 95,000 AU | cap does not bind |
| bulge, 100 pc from the centre [C] | 77 | 0.21 pc, 42,000 AU | cap does not bind |
| 10 pc from the centre [C] | 660 | 0.10 pc, 21,000 AU | cap does not bind |
| 2 pc from the centre [C] | 61,000 | 0.022 pc, 4,600 AU | cap shrinks the axis |
| a globular cluster core [R] | 10^3 to 10^5 | 0.09 to 0.009 pc, 19,000 to 1,900 AU | cap shrinks the axis |

(The density is the model's `relative_density` times `E / edge^3` times `tuning.MEAN_SYSTEM_MASS_SOLAR`, 0.4; it leaves out
the gas and dark matter, which would shorten the radii a little. The cluster row is an observed range, to be replaced by the
cluster's own stored density when the component carries one.)

**The encounter limit.** Passing stars and clouds unbind soft pairs (Heggie 1975; Bahcall, Hut and Tremaine 1985; Jiang
and Tremaine 2010). Boss's text gives about 2 x 10^4 to 5 x 10^4 AU for the solar neighbourhood over 10 Gyr, and Jiang and Tremaine's
50% survival over 10 Gyr falls near 30,000 AU (a secondary source's reading [R]); Gaia catalogues truncate near the same
scale, well inside the Jacobi radius.

The formula in Boss's pasted text could not be used as written. As typed, a_disrupt = (G M_tot^2 / (8 pi G^2 n M_*^2 ln(Lambda) T))^(1/2)
does not come out as a length (its units are mass times time), and it holds no velocity dispersion, which the diffusive
limit needs. The source papers could not be opened from this environment (arxiv.org is not reachable), so the constant of the
exact expression is unconfirmed. What the physics fixes is the scaling: the disruption time of a pair of semi-major axis a is
t = C M_tot sigma / (G rho <M_p> a ln Lambda), with sigma the velocity dispersion of the field, rho its density, <M_p> the
mean perturber mass and C a number near 10 [R, my reconstruction]. Solved for a at the stars' age T it is dimensionally a length:

    a_disrupt = a_ref * (M_tot / M_ref) * (sigma / sigma_ref) * (rho_ref / rho) * (T_ref / T)

with `a_ref` the 10 Gyr limit in the reference place. The code takes `rho` and `sigma` from the cell's profile (the density of
section 3.3 and the dispersion of its populations, `tuning.STELLAR_POPULATION_DISPERSION_KMS`) and `T` from the star's age,
so a young star in a dense cell keeps its wide pair and an old one does not, and nothing is written in by place.
`a_ref`, `rho_ref` and `sigma_ref` are one calibration: about 3 x 10^4 AU at 0.08 solar masses per pc^3 and 30 km/s (Boss's
research; Jiang and Tremaine), stored in `tuning.py` with this source named. Once the constants C and ln Lambda of the
impulse integral are confirmed, `a_ref` is computed instead of stored, and `rho_ref` goes away.

The widest axis a wide pair may have is the smaller of the tidal radius and `a_disrupt`. This takes the place of the fixed
10,000 AU cap (`WIDE_BINARY_SEPARATION_MAX_AU`), as section 6.6 of
[multistar-and-compact-systems.md](multistar-and-compact-systems.md) now says, so a solar-neighbourhood pair may be as wide as
about 30,000 AU and a pair in a dense core much narrower. Section 6.6 takes the tidal radius from the galactic potential
at the cell, the same one the orbital update uses; the local-density form above is the check on it and the fallback where
the potential is not at hand.

## 4. How the rules reach the draw

`run_sector.py` already passes `profile.kind_factors` into the fill. The three rules need the profile's populations and
density in the system draw, so the system config gets one more input, a small `PlaceTerms` record built once per cell
from the profile (the [Fe/H] per population, the young-share prevalence, the density in solar masses per cubic parsec),
read by `multistar.draw_member_tree` (section 3.2 through `prevalence`), the pair-pick (3.1) and `wide_cap` (3.3). It is
a function of the cell, so the same seed gives the same stars wherever a run lands (the GEN.39 rule). It stays out of
`tuning.py`'s stellar evolution constants, so Core 3's evolution work is not touched.

## 5. What closes GEN.103

- 3.1 and 3.3 have their constants (section 3.1's four points; section 3.3's one calibration, with the open question of the
  exact impulse constant). 3.2 needs only a design cap on the prevalence.
- Boss's pasted references are unverified here. Confirmed by search summaries: Moe, Kratter and Badenes 2019 and its four
  points; Jiang and Tremaine 2010 (arXiv 0907.2952), the 1.7 pc Jacobi radius and the survival near 30,000 AU; El-Badry and Rix 2019
  (the paper is 2019, not 2018). Not confirmed: Badenes et al. 2018's title and numbers, Bahcall, Hut and Tremaine 1985, Heggie 1975,
  El-Badry, Rix and Heintz 2021's quoted truncation, Kratter and Lodato 2016, Bate 2018, and Stegmann et al. 2024 (its title as
  pasted, "Close encounters of wide binaries induced by the Galactic tide", was not seen).
- GEN.269 and GEN.270 to GEN.275 are the multi-star items that own the code (Core 3), so this note is their input: the
  build goes there, and GEN.103 closes when they have taken the three rules in, or the rules are written as items of
  their own.
